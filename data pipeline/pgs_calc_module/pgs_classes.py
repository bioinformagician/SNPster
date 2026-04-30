import subprocess
import pandas as pd
import os
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT
from config import SCORING_FILE_SOURCE_DIR, SCORING_FILE_TARGET_DIR


class EnvironmentHandler:
    def __init__(self, 
                 output_dir: str, 
                 reference_data_path: str,
                 low_memory : bool, #split analysis into multiple runs of pgs calcs, calling NF pipeline multiple times (cache will be used)
                 scoring_file_str: str = None,
                 imputation_id: int = None,
                 prsc_id: int = None,
                 db_utils: DbUtils = None,
                 samplesheet_path: str = None,
                 scoring_file_source_dir: str = SCORING_FILE_SOURCE_DIR,
                 scoring_file_target_dir: str = SCORING_FILE_TARGET_DIR
                 ):
        
        self.samplesheet_path = samplesheet_path
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.low_memory = low_memory
        self.scoring_file_str = scoring_file_str
        self.imputation_id = imputation_id
        self.prsc_id = prsc_id
        if db_utils is None:
            db_utils = DbUtils(DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
        self.db_utils = db_utils
        self.scoring_file_source_dir = scoring_file_source_dir
        self.scoring_file_target_dir = scoring_file_target_dir
        
    
    def copy_scoring_files(self, scoring_file_list: list) -> None:
        """Copies scoring files from the mounted source directory to the scoring file directory."""
        for file in scoring_file_list:
            subprocess.run(["cp", file, self.scoring_file_target_dir], check=True)

    def connect_to_db(self) -> None:
        if not self.db_utils.db_handler.connect():
            raise ConnectionError(
                f"Failed to connect to database at {HOST}:{PORT}. "
                "Set DB_HOST and DB_PORT for the runtime environment."
            )
    
    def close_db_connection(self) -> None:
        self.db_utils.db_handler.close()
    
    def set_db_job_status(self, status: str) -> None:
        """Updates the job status in the database."""
        if self.prsc_id is not None:
            update_query = f"""UPDATE snpster_users.prsc_jobs
                            SET prsc_status = '{status}'
                            WHERE prsc_id = {self.prsc_id};"""
            self.db_utils.db_handler.execute_query(update_query)
        else:
            raise ValueError("prsc_id is not set. Cannot update job status.")


class PGSCalculator_Config:
    def __init__(self, 
                 environment_handler: EnvironmentHandler,
                 target_build: str = "GRCh38"):
        
        self.environment_handler = environment_handler
        self.target_build = target_build
        
    
        
        
        
        


class PGSCalculator:
    def __init__(self, 
                 environment_handler: EnvironmentHandler,
                 pgscalculator_config: PGSCalculator_Config):
        self.environment_handler = environment_handler
        self.pgscalculator_config = pgscalculator_config

    
    def set_job_parameters(self) -> None:


        """Fetches the list of imputation jobs from the database."""
        query = f"""WITH picked AS (
                    SELECT MIN(pj.prsc_id) AS prsc_id
                    FROM snpster_users.prsc_jobs pj
                    JOIN snpster_users.imputation_jobs ij
                    ON pj.imputation_id = ij.imputation_id
                    WHERE pj.prsc_status = 'queued'
                    AND ij.imputation_status = 'completed'
                )
                SELECT
                    pj.imputation_id,
                    pj.prsc_id,
                    ij.user_id,
                    pj.prsc_status,
                    pjp.pgs_id,
                    ij.imputation_status,
                    prs.scoring_file_path
                FROM picked p
                JOIN snpster_users.prsc_jobs pj
                ON pj.prsc_id = p.prsc_id
                JOIN snpster_users.prsc_job_parameters pjp
                ON pj.prsc_id = pjp.prsc_id
                JOIN snpster_users.imputation_jobs ij
                ON pj.imputation_id = ij.imputation_id
                JOIN snpster_users.pgs_reports_shop prs
                ON pjp.pgs_id = prs.pgs_id
                ORDER BY pj.prsc_id, pjp.pgs_id;"""
        #note to query: also need user_id to identify correct folder in srv/imputed folder. Also need to ofcourse not hardcode the imputation ID im looking for
        
        #also needs to fetch scoring file path from pgs_reports_shop

        results = self.environment_handler.db_utils.get_pd_dataframe_from_query(query)
        
        if results.empty:
            raise ValueError("No queued jobs found in the database.")
        
        prsc_id = results["prsc_id"].iloc[0]
        imputation_id = results["imputation_id"].iloc[0]
        #concatenate strings for scoring file path
        user_id = results["user_id"].iloc[0]
        self.environment_handler.prsc_id = prsc_id
        self.environment_handler.imputation_id = imputation_id
        self.environment_handler.samplesheet_path = f"/srv/imputed/{user_id}/{imputation_id}/output/samplesheet.csv"
        self.environment_handler.copy_scoring_files(results["scoring_file_path"].tolist())
        self.environment_handler.scoring_file_str = f"{self.environment_handler.scoring_file_target_dir}/*.txt.gz"
        #self.environment_handler.set_db_job_status("running") implement this when heartbeat function is implemented
        

    def upload_results(self) -> None:
        
        # path example: /output_dir/{imputation_id}/score/{imputation_id}_pgs.txt.gz
        """Uploads the PGS calculation results to the database and updates job status."""

        results_path = f"{self.environment_handler.output_dir}/{self.environment_handler.imputation_id}/score/{self.environment_handler.imputation_id}_pgs.txt.gz"

        results = pd.read_csv(results_path, sep="\t")
        results["PGS"] = results["PGS"].str.replace("_hmPOS_GRCh38", "", regex=False)
        results = results[results["sampleset"] != "reference"]
        
        for index, row in results.iterrows():
            
            #imputation_id = int(row["sampleset"])
            pgs_id = row["PGS"]
            percentile = float(row["percentile_MostSimilarPop"])
            z_most_similar_pop = float(row["Z_MostSimilarPop"])

            insert_query = f"""INSERT INTO snpster_users.prsc_job_results (prsc_id, pgs_id, percentile, z_most_similar_pop)
                            VALUES ({self.environment_handler.prsc_id}, '{pgs_id}', {percentile}, {z_most_similar_pop});"""
            
            self.environment_handler.db_utils.db_handler.execute_query(insert_query)
        
        
        


    
    def run_pgs_calculation(self) -> None:


        try:
            self.set_job_parameters()
        except ValueError as e:
            print(e)
            return


        if self.environment_handler.low_memory == "true":
            
            #get all scoring files in target dir and run NF pipeline for each of them separately, with same input and output dir (cache will be used so it should be faster than running it all together)
            scoring_file_list = [file_path for file_path in os.listdir(self.environment_handler.scoring_file_target_dir) if file_path.endswith(".txt.gz")]
            for score_file in scoring_file_list:

                command = [
                    "nextflow", "run", "/opt/pgsc_calc/main.nf",
                    "-profile", "conda",
                    "--input", self.environment_handler.samplesheet_path,
                    "--target_build", self.pgscalculator_config.target_build,
                    "--run_ancestry", self.environment_handler.reference_data_path,
                    "--outdir", self.environment_handler.output_dir,
                    "--min_overlap", "0.5",
                    "--scorefile", score_file
                ]

                print(f"Running command {command}")

                try:
                    result = subprocess.run(command, check=True)
                    print("PGS calculation completed successfully.")
                    self.upload_results() #might not work right now with current upload_results

                except subprocess.CalledProcessError as e:
                    print("Error during PGS calculation:")
                    print("Return code:", e.returncode)
        


        else:
            
            command = [
            "nextflow", "run", "/opt/pgsc_calc/main.nf",
            "-profile", "conda",
            "--input", self.environment_handler.samplesheet_path,
            "--target_build", self.pgscalculator_config.target_build,
            "--run_ancestry", self.environment_handler.reference_data_path,
            "--outdir", self.environment_handler.output_dir,
            "--min_overlap", "0.5",
            "--scorefile", self.environment_handler.scoring_file_str
        ]
            
            print(f"Running command {command}")
        
            try:
                result = subprocess.run(command, check=True)
                print("PGS calculation completed successfully.")
                self.upload_results()
                self.environment_handler.set_db_job_status("completed")  #maybe this should be a trigger function in db upon upload of data to prsc_job_results instead
                self.environment_handler.close_db_connection()
                

            except subprocess.CalledProcessError as e:
                print("Error during PGS calculation:")
                print("Return code:", e.returncode)
    