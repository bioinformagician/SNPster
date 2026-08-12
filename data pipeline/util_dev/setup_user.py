import argparse
import os
import secrets
import string
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))

from db_handler import DbHandler  # noqa: E402
from db_config import USERNAME, PASSWORD, HOST, PORT  # noqa: E402


REPORT_LIBRARY_DIR = Path(__file__).resolve().parents[1] / "pgs_libraries"


def _report_names_from_library(report_library_dir: Path = REPORT_LIBRARY_DIR) -> list[str]:
	"""Return report template names defined by the local pgs_libraries files."""
	if not report_library_dir.exists():
		raise FileNotFoundError(f"Report library folder not found: {report_library_dir}")

	report_names = []
	for report_file in sorted(report_library_dir.glob("*.txt")):
		if "all" in report_file.stem:
			continue
		report_names.append(report_file.stem)

	if not report_names:
		raise ValueError(f"No report templates found in {report_library_dir}")

	return report_names


def _report_names_from_shop(db_handler: DbHandler) -> list[str]:
	"""Return report names currently present in snpster_users.pgs_reports_shop."""
	rows = db_handler.execute_query(
		"""
			SELECT DISTINCT report_name
			FROM snpster_users.pgs_reports_shop
			ORDER BY report_name ASC;
		"""
	)
	return [str(row[0]) for row in rows if row and row[0]]


def _random_password_hash(length: int = 64) -> str:
	alphabet = string.ascii_letters + string.digits
	return "".join(secrets.choice(alphabet) for _ in range(length))


def _random_email_for_username(username: str) -> str:
	suffix = secrets.token_hex(4)
	return f"{username}.{suffix}@example.com"


def create_or_update_user(db_handler: DbHandler, username: str) -> None:
	"""Create or update user_information with synthetic email/password_hash."""
	email = _random_email_for_username(username)
	password_hash = _random_password_hash()

	upsert_user_query = """
		INSERT INTO snpster_users.user_information (user_id, email, password_hash)
		VALUES (%s, %s, %s)
		ON CONFLICT (user_id)
		DO UPDATE SET
			email = EXCLUDED.email,
			password_hash = EXCLUDED.password_hash;
	"""
	db_handler.execute_query(upsert_user_query, (username, email, password_hash))


def add_user_files(db_handler: DbHandler, username: str, file_paths: list[str]) -> list[int]:
	"""Ensure user_files rows exist for the given user/file paths and return file_ids."""
	if not file_paths:
		raise ValueError("At least one file path must be provided.")

	normalized_paths = [os.path.abspath(path) for path in file_paths]
	file_ids: list[int] = []

	for full_path in normalized_paths:
		insert_query = """
			INSERT INTO snpster_users.user_files (user_id, genefile_location)
			VALUES (%s, %s)
			ON CONFLICT DO NOTHING;
		"""
		db_handler.execute_query(insert_query, (username, full_path))

		select_query = """
			SELECT file_id
			FROM snpster_users.user_files
			WHERE user_id = %s AND genefile_location = %s
			ORDER BY file_id ASC
			LIMIT 1;
		"""
		rows = db_handler.execute_query(select_query, (username, full_path))
		if not rows:
			raise RuntimeError(f"Could not resolve file_id for path: {full_path}")

		file_ids.append(int(rows[0][0]))

	return file_ids


def create_imputation_job(db_handler: DbHandler, username: str, file_ids: list[int]) -> int:
	"""Create or consolidate to one queued imputation job and attach all provided file_ids."""
	if not file_ids:
		raise ValueError("At least one file_id is required to create an imputation job.")

	existing_job_query = """
		SELECT DISTINCT ij.imputation_id
		FROM snpster_users.imputation_jobs ij
		JOIN snpster_users.imputation_job_parameters ijp
			ON ij.imputation_id = ijp.imputation_id
		WHERE ij.user_id = %s
		  AND ij.imputation_status = 'queued'
		  AND ij.started_at IS NULL
		  AND ij.completed_at IS NULL
		  AND ijp.file_id = ANY(%s)
		ORDER BY ij.imputation_id ASC;
	"""
	existing_rows = db_handler.execute_query(existing_job_query, (username, file_ids))

	if existing_rows:
		imputation_id = int(existing_rows[0][0])
		redundant_ids = [int(row[0]) for row in existing_rows[1:]]
	else:
		insert_job_query = """
			INSERT INTO snpster_users.imputation_jobs (user_id, imputation_status)
			VALUES (%s, 'queued')
			RETURNING imputation_id;
		"""
		job_rows = db_handler.execute_query(insert_job_query, (username,))
		if not job_rows:
			raise RuntimeError("Failed to create imputation job.")
		imputation_id = int(job_rows[0][0])
		redundant_ids = []

	insert_param_query = """
		INSERT INTO snpster_users.imputation_job_parameters (imputation_id, file_id)
		VALUES (%s, %s)
		ON CONFLICT DO NOTHING;
	"""
	for file_id in file_ids:
		db_handler.execute_query(insert_param_query, (imputation_id, int(file_id)))

	if redundant_ids:
		delete_query = """
			DELETE FROM snpster_users.imputation_jobs
			WHERE imputation_id = ANY(%s);
		"""
		db_handler.execute_query(delete_query, (redundant_ids,))

	return imputation_id


def create_prsc_jobs_for_all_reports(
	db_handler: DbHandler,
	imputation_id: int,
	report_names: list[str] | None = None,
) -> list[tuple[int, str]]:
	"""
	Create one queued PRSC job per report and subscribe each job to all pgs_ids
	for that report from snpster_users.pgs_reports_shop.
	Returns list of (prsc_id, report_name).
	"""
	resolved_from_shop = report_names is None
	if resolved_from_shop:
		report_names = _report_names_from_shop(db_handler)
		if not report_names:
			raise ValueError(
				"No report_name rows found in snpster_users.pgs_reports_shop. "
				"Run setup_raw_data.py setup_pgs_reports() first."
			)

	created_jobs: list[tuple[int, str]] = []

	for report_name in report_names:
		pgs_rows = db_handler.execute_query(
			"""
				SELECT pgs_id
				FROM snpster_users.pgs_reports_shop
				WHERE report_name = %s
				ORDER BY pgs_id ASC;
			""",
			(report_name,),
		)
		if not pgs_rows:
			if resolved_from_shop:
				# Defensive guard; DISTINCT report_name query should make this rare.
				print(f"Skipping report with no shop entries: {report_name}")
				continue
			raise ValueError(f"No PGS IDs found in snpster_users.pgs_reports_shop for report {report_name}.")

		insert_prsc_query = """
			INSERT INTO snpster_users.prsc_jobs (imputation_id, prsc_status)
			VALUES (%s, 'queued')
			RETURNING prsc_id;
		"""
		prsc_rows = db_handler.execute_query(insert_prsc_query, (int(imputation_id),))
		if not prsc_rows:
			raise RuntimeError(f"Failed to create PRSC job for report {report_name}.")

		prsc_id = int(prsc_rows[0][0])

		insert_params_query = """
			INSERT INTO snpster_users.prsc_job_parameters (prsc_id, pgs_id)
			VALUES (%s, %s)
			ON CONFLICT DO NOTHING;
		"""
		for (pgs_id,) in pgs_rows:
			db_handler.execute_query(insert_params_query, (prsc_id, pgs_id))
		created_jobs.append((prsc_id, report_name))

	return created_jobs


def setup_user_with_jobs(db_handler: DbHandler, username: str, file_paths: list[str]) -> dict:
	"""
	End-to-end helper:
	1) upsert user
	2) ensure user_files exist for provided full paths
	3) create one imputation job linked to all files
	4) create one queued PRSC job per report type with all panel pgs_ids
	"""
	create_or_update_user(db_handler, username)
	file_ids = add_user_files(db_handler, username, file_paths)
	imputation_id = create_imputation_job(db_handler, username, file_ids)
	prsc_jobs = create_prsc_jobs_for_all_reports(db_handler, imputation_id)

	return {
		"username": username,
		"file_ids": file_ids,
		"imputation_id": imputation_id,
		"prsc_jobs": prsc_jobs,
	}


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Create/update a user, attach files, create imputation job, and queue all report PRSC jobs."
	)
	parser.add_argument("--username", required=True, type=str, help="User ID to set up.")
	parser.add_argument(
		"--file_paths",
		required=True,
		nargs="+",
		type=str,
		help="One or more full file paths to register in snpster_users.user_files.",
	)
	return parser.parse_args()


#frede:
#/srv/raw/frederik_tolberg_MyHeritage_raw_dna_data.zip
#/srv/raw/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt

#bullino:
#/srv/raw/sara_christensen_raw_dna_data.zip

if __name__ == "__main__":
	#args = parse_args()
	db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
	db_handler.connect()
	try:
		#result = setup_user_with_jobs(db_handler, args.username, args.file_paths)
		#print("Setup completed:")
		#print(f"username={result['username']}")
		#print(f"file_ids={result['file_ids']}")
		#print(f"imputation_id={result['imputation_id']}")
		#print("prsc_jobs=[")
		#for prsc_id, report_name in result["prsc_jobs"]:
		#	print(f"  (prsc_id={prsc_id}, report_name={report_name})")
		#print("]")
		prsc_jobs = create_prsc_jobs_for_all_reports(db_handler, 309) #bullino 309, frede all 307, 
  
	finally:
		db_handler.close()
