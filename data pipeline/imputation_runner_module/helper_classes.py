import os
import subprocess
import sys
import time
import shutil
import pandas as pd
from config import BASE_OUTPUT_DIR, IMPUTATION_DEPENDENCIES, HARMONIZER_DEPENDENCIES, NEXTFLOW_BIN, SAMPLESHEET_PATH, PIPELINE_PATH, NEXTFLOW_CONFIG, NEXTFLOW_WORK_DIR
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT




class EnvironmentHandler:
    def __init__(self, 
                 base_output_dir=BASE_OUTPUT_DIR,
                    imputation_dependencies=IMPUTATION_DEPENDENCIES,
                    harmonizer_dependencies=HARMONIZER_DEPENDENCIES,
                    nextflow_bin=NEXTFLOW_BIN,
                    nextflow_workdir=NEXTFLOW_WORK_DIR,
                    samplesheet_path=SAMPLESHEET_PATH,
                    pipeline_path=PIPELINE_PATH,
                    nextflow_config=NEXTFLOW_CONFIG,
                    samplesheet_df=None):
        
        self.base_output_dir = base_output_dir
        self.imputation_dependencies = imputation_dependencies
        self.harmonizer_dependencies = harmonizer_dependencies
        self.nextflow_bin = nextflow_bin
        self.samplesheet_path = samplesheet_path
        self.pipeline_path = pipeline_path
        self.nextflow_config = nextflow_config
        self.samplesheet_df = samplesheet_df
        self.nextflow_workdir = nextflow_workdir
        
    def validate_environment(self):
        if not self.base_output_dir:
            raise ValueError("BASE_OUTPUT_DIR is not set in config.py")
        if not self.imputation_dependencies:
            raise ValueError("IMPUTATION_DEPENDENCIES is not set in config.py")
        if not self.harmonizer_dependencies:
            raise ValueError("HARMONIZER_DEPENDENCIES is not set in config.py")

        docker_bin = shutil.which("docker")
        if docker_bin is None:
            raise RuntimeError(
                "Docker CLI not found in job_runner container. Rebuild the image so docker client is installed."
            )

        if not os.path.exists("/var/run/docker.sock"):
            raise RuntimeError(
                "Docker socket /var/run/docker.sock is not mounted. Check docker-compose volume mounts."
            )
    
        
    
class DatabaseQueryHandler:
    def __init__(self, db_utils: DbUtils = None):

        
        if db_utils is None:
            db_utils = DbUtils(db_handler = DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
        self.db_utils = db_utils
        
        
        
    def get_queued_jobs(self) -> pd.DataFrame:
        query = f"""SELECT ij.user_id, ij.imputation_id, ui.genefile_location 
                    FROM snpster_users.user_information ui
                    JOIN snpster_users.imputation_jobs ij ON ui.user_id = ij.user_id
                    WHERE ij.imputation_status = 'queued'
                    LIMIT 5;"""
        
        results = self.db_utils.get_pd_dataframe_from_query(query)
        
        return results
    
    

    def mark_jobs_running(self, imputation_ids: list) -> None:
        
        imputation_ids_sql = ", ".join(str(imputation_id) for imputation_id in imputation_ids)
        mark_running_query = f"""UPDATE snpster_users.imputation_jobs
                    SET imputation_status = 'running',
                        started_at = CURRENT_TIMESTAMP
                    WHERE imputation_id IN ({imputation_ids_sql})
                      AND imputation_status = 'queued';"""
        self.db_utils.db_handler.execute_query(mark_running_query)
        
        
    
    def mark_jobs_completed(self, imputation_ids: list) -> None:
        
        imputation_ids_sql = ", ".join(str(imputation_id) for imputation_id in imputation_ids)
        mark_completed_query = f"""UPDATE snpster_users.imputation_jobs
                                       SET imputation_status = 'completed',
                                           completed_at = CURRENT_TIMESTAMP
                                       WHERE imputation_id IN ({imputation_ids_sql});"""
                                       
        self.db_utils.db_handler.execute_query(mark_completed_query)
    
    
    
    def mark_jobs_failed(self, imputation_ids: list) -> None:
        
        imputation_ids_sql = ", ".join(str(imputation_id) for imputation_id in imputation_ids)
        mark_failed_query = f"""UPDATE snpster_users.imputation_jobs
                                    SET imputation_status = 'failed',
                                        completed_at = CURRENT_TIMESTAMP
                                    WHERE imputation_id IN ({imputation_ids_sql});"""
                                    
        self.db_utils.db_handler.execute_query(mark_failed_query)
        
        



class ImputationRunner:
    def __init__(self, env_handler: EnvironmentHandler,
                 query_handler: DatabaseQueryHandler,
                 job_df: pd.DataFrame = None):
        self.env_handler = env_handler
        self.query_handler = query_handler
        self.job_df = job_df
    
    
    def write_samplesheet(self) -> None:
        self.job_df["output_dir"] = (
            self.env_handler.base_output_dir
            + "/"
            + self.job_df["user_id"].astype(str)
            + "/"
            + self.job_df["imputation_id"].astype(str)
        )

        samplesheet_df = self.job_df[
            ["imputation_id", "output_dir", "genefile_location"]
        ].rename(
            columns={
                "imputation_id": "identifier",
                "genefile_location": "file_path",
            }
        )

        samplesheet_df.to_csv(self.env_handler.samplesheet_path, index=False)
        
        
    
    
    def run_imputation(self) -> None:
        
        print("Starting Nextflow run for imputation jobs with imputation IDs: ", self.job_df["imputation_id"].tolist())
        subprocess.run([
                        self.env_handler.nextflow_bin,
                        "run",
                        self.env_handler.pipeline_path,
                        "-c",
                        self.env_handler.nextflow_config,
                        "-work-dir",
                        self.env_handler.nextflow_workdir,
                        "--samplesheet",
                        self.env_handler.samplesheet_path,
                        "--imputation_dependencies",
                        self.env_handler.imputation_dependencies,
                        "--harmonizer_dependencies",
                        self.env_handler.harmonizer_dependencies,
                    ], check=True)

    
        
        
