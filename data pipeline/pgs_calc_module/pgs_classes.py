import subprocess
import pandas as pd
import os
import shutil
import glob
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT
from config import SCORING_FILE_SOURCE_DIR, SCORING_FILE_TARGET_DIR, NF_WORK_DIR, OUTPUT_DIR, REFERENCE_DATA_PATH


class EnvironmentHandler:
    def __init__(self, 
                 scoring_file_str: str = None,
                 imputation_ids: list = None,
                 prsc_ids: list = None,
                 db_utils: DbUtils = None,
                 samplesheet_paths: list = None,
                 combined_samplesheet_path: str = None,
                 scoring_file_source_dir: str = SCORING_FILE_SOURCE_DIR,
                 scoring_file_target_dir: str = SCORING_FILE_TARGET_DIR,
                 nf_work_dir: str = NF_WORK_DIR,
                 output_dir: str = OUTPUT_DIR,
                 reference_data_path: str = REFERENCE_DATA_PATH
                 ):
        
        self.samplesheet_paths = samplesheet_paths
        self.combined_samplesheet_path = combined_samplesheet_path
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.scoring_file_str = scoring_file_str
        self.imputation_ids = imputation_ids
        self.prsc_ids = prsc_ids
        if db_utils is None:
            db_utils = DbUtils(DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
        self.db_utils = db_utils
        self.scoring_file_source_dir = scoring_file_source_dir
        self.scoring_file_target_dir = scoring_file_target_dir
        self.nf_work_dir = nf_work_dir
    
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
    
    def set_db_job_status(self, status: str, prsc_id: str) -> None:
        """Updates the job status in the database."""

        update_query = f"""UPDATE snpster_users.prsc_jobs
                        SET prsc_status = '{status}'
                        WHERE prsc_id = '{prsc_id}';"""
        self.db_utils.db_handler.execute_query(update_query)

    
    def clear_output_directory(self) -> None:
        """Clears the output directory."""
        if not self.output_dir or not os.path.isdir(self.output_dir):
            return

        for entry in os.listdir(self.output_dir):
            path = os.path.join(self.output_dir, entry)
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)
    
    def clear_directories(self) -> None:
        """Clears the output and scoring file directories."""
        directories_to_clear = [
            self.output_dir,
            self.scoring_file_target_dir,
            self.nf_work_dir,
        ]

        for directory in directories_to_clear:
            if not directory or not os.path.isdir(directory):
                continue

            for entry in os.listdir(directory):
                path = os.path.join(directory, entry)
                if os.path.isfile(path) or os.path.islink(path):
                    os.unlink(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)




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

    
    def set_job_parameters(self) -> bool:


        """The query gets the oldest prsc job with status 'queued' and imputation job with status 'completed', and then get all other queued prsc jobs requiring the same pgs_id(s) to run together in the same NF execution, to optimize computational time from reference_population*n_prsc_ids to reference_population+n_prsc_ids"""
        query = """WITH eligible AS (
                    SELECT
                        pj.imputation_id,
                        pj.prsc_id,
                        ij.user_id,
                        pj.prsc_status,
                        ij.imputation_status
                    FROM snpster_users.prsc_jobs pj
                    JOIN snpster_users.imputation_jobs ij
                        ON pj.imputation_id = ij.imputation_id
                    WHERE pj.prsc_status = 'queued'
                      AND ij.imputation_status = 'completed'
                ),
                pgs_sets AS (
                    SELECT
                        e.prsc_id,
                        ARRAY_AGG(DISTINCT pjp.pgs_id ORDER BY pjp.pgs_id) AS pgs_id_set
                    FROM eligible e
                    JOIN snpster_users.prsc_job_parameters pjp
                        ON e.prsc_id = pjp.prsc_id
                    GROUP BY e.prsc_id
                ),
                picked AS (
                    SELECT MIN(prsc_id) AS prsc_id
                    FROM pgs_sets
                ),
                picked_set AS (
                    SELECT ps.pgs_id_set
                    FROM pgs_sets ps
                    JOIN picked p
                        ON ps.prsc_id = p.prsc_id
                ),
                matching_jobs AS (
                    SELECT ps.prsc_id
                    FROM pgs_sets ps
                    JOIN picked_set pset
                        ON ps.pgs_id_set = pset.pgs_id_set
                )
                SELECT
                    e.imputation_id,
                    e.prsc_id,
                    e.user_id,
                    e.prsc_status,
                    pjp.pgs_id,
                    e.imputation_status,
                    prs.scoring_file_path
                FROM matching_jobs mj
                JOIN eligible e
                    ON e.prsc_id = mj.prsc_id
                JOIN snpster_users.prsc_job_parameters pjp
                    ON pjp.prsc_id = e.prsc_id
                JOIN snpster_users.pgs_reports_shop prs
                    ON prs.pgs_id = pjp.pgs_id
                ORDER BY e.prsc_id, pjp.pgs_id;"""

        results = self.environment_handler.db_utils.get_pd_dataframe_from_query(query)
        
        if results.empty:
            return False
        
        df_subset = results[["prsc_id", "imputation_id", "user_id"]].drop_duplicates()
        
        
        prsc_ids = df_subset["prsc_id"].unique().tolist()
        imputation_ids = df_subset["imputation_id"].unique().tolist()
        #concatenate strings for scoring file path
        user_ids = df_subset["user_id"].unique().tolist()
        self.environment_handler.prsc_ids = prsc_ids
        self.environment_handler.imputation_ids = imputation_ids
        self.environment_handler.user_ids = user_ids
        self.environment_handler.samplesheet_paths = [
        f"/srv/imputed/{row.user_id}/{row.imputation_id}/output/samplesheet.csv"
        for row in df_subset.itertuples(index=False)
                                                        ]
        self.environment_handler.copy_scoring_files(sorted(set(results["scoring_file_path"].tolist())))
        self.environment_handler.scoring_file_str = f"{self.environment_handler.scoring_file_target_dir}/*.txt.gz"

        return True
    
    def copy_imputed_files(self) -> None:
        """Copies imputed files from the mounted imputed directory flat into nf_work_dir."""
        for idx, imputation_id in enumerate(self.environment_handler.imputation_ids):
            source_dir = f"/srv/imputed/{self.environment_handler.user_ids[idx]}/{imputation_id}/output/"
            for entry in os.scandir(source_dir):
                if entry.is_file() and not entry.name.endswith(".csv"):
                    print(f"Copying {entry.path} to {self.environment_handler.nf_work_dir}")
                    shutil.copy2(entry.path, self.environment_handler.nf_work_dir)
    
    
    
    def create_combined_samplesheet(self) -> None:
        """Combines all samplesheets from the required imputation jobs into a single samplesheet for the NF pipeline."""

        combined_samplesheet_path = os.path.join(self.environment_handler.nf_work_dir, "samplesheet.csv")
        
        samplesheet_dfs = []
        for samplesheet_path in self.environment_handler.samplesheet_paths:
            print(f"Reading samplesheet from {samplesheet_path}")
            df = pd.read_csv(samplesheet_path)
            samplesheet_dfs.append(df)

        combined_df = pd.concat(samplesheet_dfs, ignore_index=True)
        combined_df.to_csv(combined_samplesheet_path, index=False)
        self.environment_handler.combined_samplesheet_path = combined_samplesheet_path

    

    def upload_results(self, imputation_id, prsc_id) -> None:
        
        # path example: /output_dir/{imputation_id}/score/{imputation_id}_pgs.txt.gz
        """Uploads the PGS calculation results to the database and updates job status."""

        results_path = f"{self.environment_handler.output_dir}/{imputation_id}/score/{imputation_id}_pgs.txt.gz"

        results = pd.read_csv(results_path, sep="\t")
        
        results["PGS"] = results["PGS"].str.replace("_hmPOS_GRCh38", "", regex=False)
        results = results[results["sampleset"] != "reference"]
        
        for index, row in results.iterrows():
            
            #imputation_id = int(row["sampleset"])
            pgs_id = row["PGS"]
            percentile = float(row["percentile_MostSimilarPop"])
            z_most_similar_pop = float(row["Z_MostSimilarPop"])

            insert_query = f"""INSERT INTO snpster_users.prsc_job_results (prsc_id, pgs_id, percentile, z_most_similar_pop)
                            VALUES ({prsc_id}, '{pgs_id}', {percentile}, {z_most_similar_pop});"""
            
            self.environment_handler.db_utils.db_handler.execute_query(insert_query)
        
    def validate_results(self, imputation_id) -> bool:
        """Validates the PGS calculation results before uploading."""
        
        results_path = f"{self.environment_handler.output_dir}/{imputation_id}/score/{imputation_id}_pgs.txt.gz"
        
        try:
            results = pd.read_csv(results_path, sep="\t")
        except Exception as e:
            print(f"Error reading results file: {e}")
            return False
            
        if results.empty:
            print("No results found in the output file.")
            return False
        
        return True



    
    def run_pgs_calculation(self, sample_sheet: str) -> None:
        
            
        command = [
        "nextflow", "run", "/opt/pgsc_calc/main.nf",
        "-profile", "conda",
        "--input", sample_sheet,
        "--target_build", self.pgscalculator_config.target_build,
        "--run_ancestry", self.environment_handler.reference_data_path,
        "--outdir", self.environment_handler.output_dir,
        "--min_overlap", "0.01",
        "--scorefile", self.environment_handler.scoring_file_str
    ]
        
        print(f"Running command {command}")
    
        try:
            subprocess.run(command, check=True)
            print("PGS calculation completed successfully.")
            
            

        except subprocess.CalledProcessError as e:
            print("Error during PGS calculation:")
            print("Return code:", e.returncode)
    