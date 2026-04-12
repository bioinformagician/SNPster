import os
import subprocess
import sys
import time
from pathlib import Path

from config import BASE_OUTPUT_DIR, IMPUTATION_DEPENDENCIES, HARMONIZER_DEPENDENCIES, NEXTFLOW_BIN, SAMPLESHEET_PATH



sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))

from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT


if __name__ == "__main__":

    if not BASE_OUTPUT_DIR:
        raise ValueError("BASE_OUTPUT_DIR is not set in config.py")
    if not IMPUTATION_DEPENDENCIES:
        raise ValueError("IMPUTATION_DEPENDENCIES is not set in config.py")
    if not HARMONIZER_DEPENDENCIES:
        raise ValueError("HARMONIZER_DEPENDENCIES is not set in config.py")

    while True:

        db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
        db_handler.connect()

        db_utils = DbUtils(db_handler)

        query = f"""SELECT ij.user_id, ij.imputation_id, ui.genefile_location 
                    FROM snpster_users.user_information ui
                    JOIN snpster_users.imputation_jobs ij ON ui.user_id = ij.user_id
                    WHERE ij.imputation_status = 'queued'
                    LIMIT 5;"""
        
        results = db_utils.get_pd_dataframe_from_query(query)

        print(f"Found {len(results)} queued jobs to process.")

        results["output_dir"] = (
            BASE_OUTPUT_DIR
            + "/"
            + results["user_id"].astype(str)
            + "/"
            + results["imputation_id"].astype(str)
        )

        #write samplesheet which has the columns: identifier, output_dir, file_path (identifier = imputation_id)

        samplesheet_df = results[
            ["imputation_id", "output_dir", "genefile_location"]
        ].rename(
            columns={
                "imputation_id": "identifier",
                "genefile_location": "file_path",
            }
        )

        samplesheet_df.to_csv(SAMPLESHEET_PATH, index=False)
        user_id = results.iloc[0]["user_id"]
        #output_dir = f"{BASE_OUTPUT_DIR}/{user_id}"

        subprocess.run([NEXTFLOW_BIN, "run", "/home/frederik/github_projects/SNPster/data pipeline/nextflow/pipeline_orchestrator.nf",
                        "--samplesheet", SAMPLESHEET_PATH,
                        "--imputation_dependencies", IMPUTATION_DEPENDENCIES,
                        "--harmonizer_dependencies", HARMONIZER_DEPENDENCIES], check=True)

        
        update_query = f"""UPDATE snpster_users.imputation_jobs
                            SET imputation_status = 'processed' 
                            WHERE user_id = '{user_id}' AND imputation_id in ({', '.join([f"'{imputation_id}'" for imputation_id in results['imputation_id']])});"""
        db_handler.execute_query(update_query)