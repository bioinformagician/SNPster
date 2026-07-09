import os
import shutil
from pathlib import Path
import sys
import re
from urllib.parse import urlparse
from urllib.request import urlopen
import random as rd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))
print(sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module")))
from db_handler import DbHandler
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH


RAW_DATA_DIR = "/home/frederik/snpster_project/zipped"
TARGET_DIR = "/srv/raw"



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


def setup_pgs_reports():


    # read data:
    report_library_folder = "/home/frederik/github_projects/SNPster/data pipeline/pgs_libraries"
    scoring_target_dir = "/srv/scoring_files"
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
        if existing_files:
            existing_path = os.path.join(scoring_target_dir, existing_files[0])
            if os.path.getsize(existing_path) > 0:
                print(f"Scoring file already exists for {pgs_id}: {existing_path}")
                return existing_path
            print(f"Existing scoring file for {pgs_id} is empty. Re-downloading: {existing_path}")

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

    for filename in os.listdir(report_library_folder):

        if "all" in filename:
            continue  # skip the all file, which is just a combined file of all reports for testing purposes

        if not filename.endswith(".txt"):
            continue

        report_name = filename.replace(".txt", "")
        file_path = os.path.join(report_library_folder, filename)

        with open(file_path, "r", encoding="utf-8") as f:
            report_content = f.read()

        pgs_ids = extract_pgs_ids(report_content)
        if not pgs_ids:
            print(f"No PGS IDs found in report file: {file_path}")
            continue

        for pgs_id in pgs_ids:
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
                print(f"No ftp_link found for {pgs_id}, skipping download and DB update.")
                continue

            ftp_link = ftp_rows[0][0]

            try:
                local_scoring_file = download_scoring_file(ftp_link, pgs_id)
            except Exception as exc:
                print(f"Failed to download scoring file for {pgs_id} from {ftp_link}: {exc}")
                continue

            insert_query = """
            INSERT INTO snpster_users.pgs_reports_shop (pgs_id, report_name, scoring_file_path)
            VALUES (%s, %s, %s)
            ON CONFLICT DO NOTHING;
            """
            db_handler.execute_query(insert_query, (pgs_id, report_name, local_scoring_file))

            update_query = """
            UPDATE snpster_users.pgs_reports_shop
            SET scoring_file_path = %s
            WHERE pgs_id = %s AND report_name = %s;
            """
            db_handler.execute_query(update_query, (local_scoring_file, pgs_id, report_name))

            print(
                f"Report '{report_name}' mapped to {pgs_id}; scoring file stored at {local_scoring_file}"
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
    

    #populate the prsc_job_parameters table with pgs_ids from pgs_reports_shop where report_name = 'cardiovascular_panel'
    panels = ['blood_panel_pgs_ids', 'cancer_pgs_ids', 'cardiovascular_pgs_ids', 'immune_and_autoimmune_pgs_ids', 'metabolic_and_endocrine_pgs_ids', 'neurological_and_psychiatric_pgs_ids', 'ophthalmology_pgs_ids', 'other_pgs_ids', 'respiratory_pgs_ids']
    

    
    random_panel_query = """
        INSERT INTO snpster_users.prsc_job_parameters (prsc_id, pgs_id)
        SELECT
            pj.prsc_id,
            prs.pgs_id
        FROM snpster_users.prsc_jobs AS pj
        CROSS JOIN LATERAL (
            SELECT panel_name
            FROM unnest(%s::text[]) AS panel_name
            ORDER BY random() + (pj.prsc_id * 0)
            LIMIT 1
        ) AS picked_panel
        JOIN snpster_users.pgs_reports_shop AS prs
            ON prs.report_name = picked_panel.panel_name
        WHERE pj.imputation_id = ANY(%s)
        ON CONFLICT DO NOTHING;
        """
    
    db_handler.execute_query(random_panel_query, (panels, imputation_ids))
    print("Populated prsc_job_parameters table with PGS IDs for reports.")
    
    
    
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
    
    db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
    db_handler.connect()
    #update_ftp_links_to_grch38()
    #transfer_files(RAW_DATA_DIR, TARGET_DIR)
    setup_pgs_reports()
    setup_prsc_jobs()
