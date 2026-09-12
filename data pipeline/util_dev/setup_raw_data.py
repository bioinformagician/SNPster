import os
import shutil
from pathlib import Path
import sys
import re
import gzip
import csv
from urllib.parse import urlparse
from urllib.request import urlopen
import random as rd
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))
print(sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module")))
from db_handler import DbHandler
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH


RAW_DATA_DIR = "/home/frederik/snpster_project/zipped"
TARGET_DIR = "/srv/raw"


def validate_scoring_file_effect_weight(scorefile_path: str) -> tuple[bool, str]:
    """Return (is_valid, reason) after checking numeric score values."""
    if not os.path.exists(scorefile_path):
        return False, "file not found"

    opener = gzip.open if scorefile_path.endswith(".gz") else open
    try:
        with opener(scorefile_path, "rt", encoding="utf-8", errors="replace") as handle:
            # PGS files contain metadata lines prefixed with '#'; the actual table
            # header is the first non-comment, non-empty line.
            header_fields = None
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                header_fields = [part.strip() for part in raw_line.rstrip("\n\r").split("\t")]
                break

            if not header_fields:
                return False, "missing table header row"

            columns_by_name = {
                name.lstrip("#").strip(): name
                for name in header_fields
                if name
            }
            effect_weight_col = columns_by_name.get("effect_weight")

            if effect_weight_col is None:
                return False, "missing effect_weight column"

            numeric_columns = {
                name: columns_by_name[name]
                for name in ("effect_weight", "OR", "HR")
                if name in columns_by_name
            }
            reader = csv.DictReader(handle, delimiter="\t", fieldnames=header_fields)

            for row_nr, row in enumerate(reader, start=1):
                if not any(value and str(value).strip() for value in row.values()):
                    continue

                for column_name, column in numeric_columns.items():
                    raw = row.get(column)
                    value = "" if raw is None else str(raw).strip()
                    if value == "":
                        return False, f"empty {column_name} at data row {row_nr}"
                    try:
                        float(value)
                    except ValueError:
                        return False, f"non-numeric {column_name} '{value}' at data row {row_nr}"
    except Exception as exc:
        return False, f"read error: {exc}"

    return True, ""




def transfer_files(source_dir:str, target_dir:str, db_handler:DbHandler) -> None:
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)



    for _, filename in enumerate(os.listdir(source_dir)):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            if os.path.exists(target_file) and os.path.getsize(target_file) > 0:
                print(f"File already exists, skipping copy: {target_file}")
            else:
                shutil.copy2(source_file, target_file)
                print(f"Copied: {source_file} to {target_file}")

            user_id = _
            email = f"{_}@example.com"
            password_hash = "hashed_password"  # Placeholder, replace with actual hash if needed
            genefile_location = target_file

            upsert_user_query = """
            INSERT INTO snpster_users.user_information (user_id, email, password_hash)
            VALUES (%s, %s, %s)
            ON CONFLICT (user_id)
            DO UPDATE SET
                email = EXCLUDED.email,
                password_hash = EXCLUDED.password_hash;
            """
            db_handler.execute_query(upsert_user_query, (str(user_id), email, password_hash))

            insert_file_query = """
            INSERT INTO snpster_users.user_files (user_id, genefile_location)
            SELECT %s, %s
            WHERE NOT EXISTS (
                SELECT 1
                FROM snpster_users.user_files
                WHERE user_id = %s AND genefile_location = %s
            );
            """
            db_handler.execute_query(insert_file_query, (str(user_id), genefile_location, str(user_id), genefile_location))
            print(f"Upserted user and ensured user_file row for {target_file}")

            # Trigger on user_files inserts creates a queued imputation job and imputation_job_parameters row.


def setup_pgs_reports(
    categorized_scores_file: str = "/home/frederik/github_projects/SNPster/data pipeline/pgs_libraries/pgs_disease_relevant_categorized.csv",
):
    """Download and register disease-relevant PGS scores from one categorized CSV."""

    scoring_target_dir = "/srv/scoring_files"
    required_columns = {"pgs_id", "category"}

    try:
        categorized_scores = pd.read_csv(categorized_scores_file, dtype=str)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Categorized PGS file not found: {categorized_scores_file}") from exc

    missing_columns = required_columns.difference(categorized_scores.columns)
    if missing_columns:
        raise ValueError(
            f"Categorized PGS file is missing required columns: {sorted(missing_columns)}"
        )

    categorized_scores = categorized_scores[["pgs_id", "category"]].dropna()
    categorized_scores["pgs_id"] = categorized_scores["pgs_id"].str.strip()
    categorized_scores["category"] = categorized_scores["category"].str.strip()
    categorized_scores = categorized_scores[
        (categorized_scores["pgs_id"] != "") & (categorized_scores["category"] != "")
    ].drop_duplicates()

    if categorized_scores.empty:
        raise ValueError(f"No valid pgs_id/category rows found in {categorized_scores_file}")

    invalid_ids = categorized_scores.loc[
        ~categorized_scores["pgs_id"].str.fullmatch(r"PGS\d{6}"), "pgs_id"
    ].unique()
    if len(invalid_ids):
        raise ValueError(f"Invalid PGS IDs in {categorized_scores_file}: {', '.join(invalid_ids[:10])}")

    os.makedirs(scoring_target_dir, exist_ok=True)

    def extract_pgs_ids(content: str) -> list:
        # Handles both comma-separated files and free text lines containing PGS IDs.
        return sorted(set(re.findall(r"PGS\d{6}", content)))

    def download_scoring_file(ftp_link: str, pgs_id: str) -> str:
        parsed = urlparse(ftp_link)
        basename = os.path.basename(parsed.path) or f"{pgs_id}.txt.gz"
        local_name = f"{pgs_id}_{basename}"
        local_path = os.path.join(scoring_target_dir, local_name)

        # Reuse any existing non-empty file for this PGS ID across reruns.
        existing_files = sorted(
            f for f in os.listdir(scoring_target_dir)
            if f.startswith(f"{pgs_id}_")
        )

        # Prefer a cached file that matches the exact expected (GRCh38 harmonized) filename.
        expected_cached_path = os.path.join(scoring_target_dir, local_name)
        if os.path.exists(expected_cached_path) and os.path.getsize(expected_cached_path) > 0:
            print(f"Scoring file already exists for {pgs_id}: {expected_cached_path}")
            return expected_cached_path

        if existing_files:
            print(
                f"Found existing cached files for {pgs_id} that do not match expected target file "
                f"{local_name}. Downloading expected file instead."
            )

        if not os.path.exists(local_path):
            print(f"Downloading scoring file for {pgs_id} from {ftp_link}")
            with urlopen(ftp_link) as response, open(local_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)
        else:
            if os.path.getsize(local_path) > 0:
                print(f"Scoring file already exists for {pgs_id}: {local_path}")
            else:
                print(f"Scoring file exists but is empty for {pgs_id}. Re-downloading: {local_path}")
                with urlopen(ftp_link) as response, open(local_path, "wb") as out_file:
                    shutil.copyfileobj(response, out_file)

        return local_path

    score_files: dict[str, str] = {}
    for pgs_id in categorized_scores["pgs_id"].unique():
        ftp_query = """
        SELECT ftp_link
        FROM data_libraries.pgscatalog_data
        WHERE pgs_id = %s
          AND ftp_link IS NOT NULL
          AND ftp_link <> ''
        LIMIT 1;
        """
        ftp_rows = db_handler.execute_query(ftp_query, (pgs_id,))
        if not ftp_rows:
            print(f"No catalog FTP link found for {pgs_id}; it will not be registered.")
            continue

        ftp_link = ftp_rows[0][0]
        if "_hmPOS_GRCh38" not in ftp_link:
            ftp_link = (
                f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/"
                f"ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh38.txt.gz"
            )

        try:
            local_scoring_file = download_scoring_file(ftp_link, pgs_id)
        except Exception as exc:
            print(f"Failed to download {pgs_id} from {ftp_link}: {exc}")
            continue

        is_valid, reason = validate_scoring_file_effect_weight(local_scoring_file)
        if not is_valid:
            print(f"Skipping invalid score file for {pgs_id}: {local_scoring_file} ({reason})")
            continue
        score_files[pgs_id] = local_scoring_file

    db_handler.execute_query("DELETE FROM snpster_users.pgs_reports_shop;")

    insert_query = """
    INSERT INTO snpster_users.pgs_reports_shop (pgs_id, report_name, scoring_file_path)
    VALUES (%s, %s, %s);
    """
    registered = 0
    for row in categorized_scores.itertuples(index=False):
        scoring_file = score_files.get(row.pgs_id)
        if scoring_file is None:
            continue
        db_handler.execute_query(insert_query, (row.pgs_id, row.category, scoring_file))
        registered += 1

    print(
        f"Rebuilt snpster_users.pgs_reports_shop from {categorized_scores_file}: "
        f"{registered} category mappings for {len(score_files)} downloaded scores."
    )

def setup_prsc_jobs():
    

    #create an insert to have jobs available
    
    
    imputation_id_query = """SELECT imputation_id
                    FROM snpster_users.imputation_jobs
                    WHERE imputation_status = 'completed';"""
    
    imputation_id_rows = db_handler.execute_query(imputation_id_query)

    imputation_ids = [row[0] for row in imputation_id_rows]
    print(f"Using existing completed imputation IDs for PRSC job setup: {imputation_ids}")

    prsc_status = "queued"

    for imputation_id in imputation_ids:
        insert_query = """
        INSERT INTO snpster_users.prsc_jobs (imputation_id, prsc_status)
        SELECT %s, %s
        WHERE NOT EXISTS (
            SELECT 1
            FROM snpster_users.prsc_jobs
            WHERE imputation_id = %s
        );
        """

        db_handler.execute_query(insert_query, (imputation_id, prsc_status, imputation_id))
        print(f"Inserted prsc job for imputation_id {imputation_id} with status {prsc_status}")
    

    all_panels_query = """
        INSERT INTO snpster_users.prsc_job_parameters (prsc_id, pgs_id)
        SELECT DISTINCT
            pj.prsc_id,
            prs.pgs_id
        FROM snpster_users.prsc_jobs AS pj
        JOIN snpster_users.pgs_reports_shop AS prs
            ON TRUE
        WHERE pj.imputation_id = ANY(%s)
          AND NOT EXISTS (
              SELECT 1
              FROM snpster_users.prsc_job_parameters AS existing
              WHERE existing.prsc_id = pj.prsc_id
                AND existing.pgs_id = prs.pgs_id
          )
        ON CONFLICT DO NOTHING;
        """
    
    db_handler.execute_query(all_panels_query, (imputation_ids,))
    print("Populated prsc_job_parameters table with PGS IDs for all panels.")
    
    
    
def update_ftp_links_to_grch38():
    #example link https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000001/ScoringFiles/Harmonized/PGS000001_hmPOS_GRCh38.txt.gz
    
    """CREATE TABLE data_libraries.pgscatalog_data (
        pgs_id varchar(100) PRIMARY KEY,
        pgs_name VARCHAR(255),
        reported_trait VARCHAR(255),
        mapped_trait_efo_label VARCHAR(255),
        efo_id VARCHAR(255),
        pgs_development_method VARCHAR(255),
        pgs_development_details TEXT,
        original_genome_build VARCHAR(20),
        number_of_variants INTEGER,
        number_of_interaction_terms INTEGER,
        type_of_variant_weight TEXT,
        pgp_id varchar(100),
        publication_pmid int,
        publication_doi VARCHAR(255),
        score_and_results_match_original_publication BOOLEAN,
        ancestry_distribution_source_of_variant_associations_gwas VARCHAR(255),
        ancestry_distribution_score_development_training VARCHAR(255),
        ancestry_distribution_pgs_evaluation VARCHAR(255),
        ftp_link VARCHAR(255),
        release_date DATE,
        license_terms_of_use TEXT
    );"""
    
    #get all rows in table
    rows = db_handler.execute_query("SELECT pgs_id FROM data_libraries.pgscatalog_data;")
    if not rows:
        print("No PGS IDs found in pgscatalog_data table.")
        return
    pgs_ids = [row[0] for row in rows]

    for pgs_id in pgs_ids:
        ftp_link = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh38.txt.gz"
        
        #check if ftp_link is valid
        try:
            with urlopen(ftp_link) as response:
                if response.status == 200:
                    print(f"FTP link is valid for {pgs_id}: {ftp_link}")
                else:
                    print(f"FTP link returned status {response.status} for {pgs_id}: {ftp_link}")
                    continue
        except Exception as exc:
            print(f"Error accessing FTP link for {pgs_id}: {ftp_link} - {exc}")
            continue
        
        update_query = """
        UPDATE data_libraries.pgscatalog_data
        SET ftp_link = %s
        WHERE pgs_id = %s;
        """
        db_handler.execute_query(update_query, (ftp_link, pgs_id))
        print(f"Updated ftp_link for {pgs_id} to {ftp_link}")


if __name__ == "__main__":
    with DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST) as db_handler:
        update_ftp_links_to_grch38()
        transfer_files(RAW_DATA_DIR, TARGET_DIR, db_handler)
        setup_pgs_reports()
        setup_prsc_jobs()
