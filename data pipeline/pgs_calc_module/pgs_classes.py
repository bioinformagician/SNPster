import subprocess
import pandas as pd
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT



class EnvironmentHandler:
    def __init__(self, samplesheet_path: str, 
                 output_dir: str, 
                 reference_data_path: str,
                 low_memory : bool, #split analysis into multiple runs of pgs calcs, calling NF pipeline multiple times (cache will be used)
                 scoring_file_str: str = None,
                 imputation_id: int = None,
                 db_utils = DbUtils(DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
                 ):
        
        self.samplesheet_path = samplesheet_path
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.low_memory = low_memory
        self.scoring_file_str = scoring_file_str
        self.imputation_id = imputation_id
        self.db_utils = db_utils


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

    
    def get_jobs(self) -> pd.DataFrame:


        """Fetches the list of imputation jobs from the database."""
        query = f"""SELECT pj.imputation_id, pj.prsc_ud, pj.prsc_status, pjp.prsc_id, pjp.pgs_id, ij.imputation_status, prs.scoring_file_path
                    FROM snpster_users.prsc_jobs pj
                    JOIN snpster_users.prsc_job_parameters pjp ON pj.prsc_id = pjp.prsc_id
                    JOIN snpster_users.imputation_jobs ij ON pj.imputation_id = ij.imputation_id
                    JOIN snpster_users.pgs_reports_shop prs ON pjp.pgs_id = prs.pgs_id
                    WHERE pj.prsc_status = 'queued'
                    AND ij.imputation_status = 'completed'
                    LIMIT 1;"""
        
        #also needs to fetch scoring file path from pgs_reports_shop

        results = self.environment_handler.db_utils.get_pd_dataframe_from_query(query)
        
        if results.empty:
            raise ValueError("No queued jobs found in the database.")
        
        self.environment_handler.imputation_id = results["imputation_id"].iloc[0]
        #concatenate strings for scoring file path

        #self.environment_handler.scoring_file_str = results["scoring_file_path"].iloc[0] 

        self.environment_handler.scoring_file_str = "/srv/scoring_files/*.txt.gz" #i dont know if it will run all of them or pick correct ones
        return results

    def upload_results(self) -> None:
        
        # path example: /output_dir/{imputation_id}/score/{imputation_id}_pgs.txt.gz
        """Uploads the PGS calculation results to the database and updates job status."""

        results_path = f"{self.environment_handler.output_dir}/{self.environment_handler.imputation_id}/score/{self.environment_handler.imputation_id}_pgs.txt.gz"

        results = pd.read_csv(results_path, sep="\t")
        results["PGS"] = results["PGS"].str.replace("_hmPOS_GRCh38", "", regex=False)
        results = results[results["sampleset"] != "reference"]

        for index, row in results.iterrows():
            
            prsc_id = int(row["sampleset"])
            pgs_id = row["PGS"]
            percentile = float(row["percentile"])
            z_most_similar_pop = float(row["z_most_similar_pop"])

            insert_query = f"""INSERT INTO snpster_users.prsc_job_results (prsc_id, pgs_id, percentile, z_most_similar_pop)
                            VALUES ({prsc_id}, '{pgs_id}', {percentile}, {z_most_similar_pop});"""
            
            self.environment_handler.db_handler.execute_query(insert_query)

        


    
    def run_pgs_calculation(self) -> None:


        try:
            job_df = self.get_jobs()
        except ValueError as e:
            print(e)
            return
            
        pgs_id_str = ",".join(job_df['pgs_id'].tolist())
        
        """nextflow run pgscatalog/pgsc_calc \
            -profile <docker/singularity/conda> \
            --input samplesheet.csv --target_build GRCh37 \
            --pgs_id PGS001229 \
            --run_ancestry pgsc_HGDP+1kGP_v1.tar.zst"""
        
        command = [
            "nextflow", "run", "pgscatalog/pgsc_calc",
            "-profile", "singularity",
            "--input", self.environment_handler.samplesheet_path,
            "--target_build", self.pgscalculator_config.target_build,
            "--pgs_id", pgs_id_str,
            "--run_ancestry", self.environment_handler.reference_data_path,
            "--outdir", self.environment_handler.output_dir,
            "--min_overlap", "0.5",
            "--scorefile", self.environment_handler.scoring_file_str
        ]


        if self.environment_handler.low_memory == "true":

            pgs_id_list = pgs_id_str.split(",")

            for pgs_id in pgs_id_list:

                command = [
                    "nextflow", "run", "pgscatalog/pgsc_calc",
                    "-profile", "singularity",
                    "--input", self.environment_handler.samplesheet_path,
                    "--target_build", self.pgscalculator_config.target_build,
                    "--pgs_id", pgs_id,
                    "--run_ancestry", self.environment_handler.reference_data_path,
                    "--outdir", self.environment_handler.output_dir,
                    "--min_overlap", "0.5",
                    "--scorefile", self.environment_handler.scoring_file_str
                ]

                print(f"Running command {command}")

                try:
                    result = subprocess.run(command, check=True, capture_output=True, text=True)
                    print("PGS calculation completed successfully.")
                    print("Output:", result.stdout)

                except subprocess.CalledProcessError as e:
                    print("Error during PGS calculation:")
                    print("Return code:", e.returncode)
                    print("Output:", e.output)
                    print("Error message:", e.stderr)
            
                self.upload_results()
        


        else:
            
            command = [
            "nextflow", "run", "pgscatalog/pgsc_calc",
            "-profile", "singularity",
            "--input", self.environment_handler.samplesheet_path,
            "--target_build", self.pgscalculator_config.target_build,
            "--pgs_id", pgs_id_str,
            "--run_ancestry", self.environment_handler.reference_data_path,
            "--outdir", self.environment_handler.output_dir,
            "--min_overlap", "0.5",
            "--scorefile", self.environment_handler.scoring_file_str
        ]
            
            print(f"Running command {command}")
        
            try:
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                print("PGS calculation completed successfully.")
                print("Output:", result.stdout)

            except subprocess.CalledProcessError as e:
                print("Error during PGS calculation:")
                print("Return code:", e.returncode)
                print("Output:", e.output)
                print("Error message:", e.stderr)

            self.upload_results()
    