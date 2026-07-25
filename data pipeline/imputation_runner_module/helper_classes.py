import os
import subprocess
import shutil
import re
from pathlib import Path
import pandas as pd
from config import BASE_OUTPUT_DIR, IMPUTATION_DEPENDENCIES, HARMONIZER_DEPENDENCIES, NEXTFLOW_BIN, SAMPLESHEET_PATH, PIPELINE_PATH, NEXTFLOW_CONFIG, NEXTFLOW_WORK_DIR
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT
from vcf_classes import VCFUtilities




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
    
    def clear_environment(self) -> None:
        """Clears the output and scoring file directories."""
        
        
        self.samplesheet_df = None
        
        
        # Keep full Nextflow task history for debugging/resume.
        # Only clear runner-owned staging dirs under the work dir.
        directories_to_clear = [
            os.path.join(self.nextflow_workdir, "std"),
            os.path.join(self.nextflow_workdir, "harm"),
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
                    
        
        if os.path.exists(self.samplesheet_path):
            os.remove(self.samplesheet_path)
        
    
class DatabaseQueryHandler:
    def __init__(self, db_utils: DbUtils = None):

        
        if db_utils is None:
            db_utils = DbUtils(db_handler = DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
        self.db_utils = db_utils
        
        
        
    def get_queued_jobs(self) -> pd.DataFrame:
        
        
        query = f"""SELECT DISTINCT ON (ij.imputation_id)
                        uf.user_id,
                        ij.imputation_id,
                        ijp.file_id,
                        uf.genefile_location
                    FROM snpster_users.imputation_jobs ij
                    JOIN snpster_users.imputation_job_parameters ijp
                        ON ij.imputation_id = ijp.imputation_id
                    JOIN snpster_users.user_files uf
                        ON uf.file_id = ijp.file_id
                    WHERE ij.imputation_status = 'queued'
                    ORDER BY ij.imputation_id ASC, ijp.file_id ASC
                    LIMIT 33;"""
        
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
    
    
    def upload_imputation_info(self, imputation_id: int, file_path: str, number_of_variants: int, chromosome: int) -> None:
        """Uploads imputation info to the database."""
        insert_query = f"""
            INSERT INTO snpster_users.imputed_data (imputation_id, file_path, number_of_variants, chromosome)
            VALUES ({imputation_id}, '{file_path}', {number_of_variants}, {chromosome})
            ON CONFLICT (imputation_id, chromosome) DO UPDATE
            SET file_path = EXCLUDED.file_path,
                number_of_variants = EXCLUDED.number_of_variants;
        """
        
        self.db_utils.db_handler.execute_query(insert_query)

        
        
        
        
        
        


class ImputationRunner:
    
    def __init__(self, env_handler: EnvironmentHandler,
                 query_handler: DatabaseQueryHandler,
                 vcf_utilities: VCFUtilities,
                 job_df: pd.DataFrame = None):
        
        self.env_handler = env_handler
        self.query_handler = query_handler
        self.job_df = job_df
        self.vcf_utilities = vcf_utilities
        
    
    def set_job_df(self) -> None:
        self.job_df = self.query_handler.get_queued_jobs()
    
    
    def write_samplesheet(self) -> None:
        
        """Write samplesheet like:
            identifier,output_dir,file_path
            76,/srv/imputed/75/76,/srv/raw/genome_Aaron_Hill_v3_Full_20191101162607.zip
            77,/srv/imputed/76/77,/srv/raw/23andMe_results.zip
            78,/srv/imputed/77/78,/srv/raw/genome_Orlando_Montalvo_Full_20140522160413.zip
            79,/srv/imputed/78/79,/srv/raw/dna-data-2013-11-09.zip
            80,/srv/imputed/79/80,/srv/raw/genome_Full_20141024183341.zip
        """
        
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
        pipeline_path = self.env_handler.pipeline_path

        standardizer_output_dir = os.path.join(self.env_handler.nextflow_workdir, "std")
        harmonizer_output_dir = os.path.join(self.env_handler.nextflow_workdir, "harm")
        os.makedirs(standardizer_output_dir, exist_ok=True)
        os.makedirs(harmonizer_output_dir, exist_ok=True)

        subprocess.run([
                        self.env_handler.nextflow_bin,
                        "run",
                        pipeline_path,
                        "-c",
                        self.env_handler.nextflow_config,
                        "-w",
                        self.env_handler.nextflow_workdir,
                        "--samplesheet",
                        self.env_handler.samplesheet_path,
                        "--imputation_dependencies",
                        self.env_handler.imputation_dependencies,
                        "--harmonizer_dependencies",
                        self.env_handler.harmonizer_dependencies,
                        "--standardizer_output_dir",
                        standardizer_output_dir,
                        "--harmonizer_output_dir",
                        harmonizer_output_dir,
                    ], check=True)


    
    def evaluate_results(self) -> tuple[list[int], list[int]]:
        """Evaluate the results of the imputation by checking output dir for files, if 22 files present -> success otherwise fail"""

        successful_ids = []
        failed_ids = []
        for _, row in self.job_df.iterrows():
            imputation_id = row["imputation_id"]
            output_dir = row["output_dir"]
            output_dir = os.path.join(output_dir, "qc_output")

            if not os.path.isdir(output_dir):
                failed_ids.append(imputation_id)
                continue

            file_count = sum(len(files) for _, _, files in os.walk(output_dir))
            if file_count == 22:
                successful_ids.append(imputation_id)
            else:
                failed_ids.append(imputation_id)
        
        #remove failed ids from self.job_df
        self.job_df = self.job_df[~self.job_df["imputation_id"].isin(failed_ids)]
        
        
        return successful_ids, failed_ids
    
    
    

    def get_output_imputation_ids(self) -> set[int]:
        """Extract imputation IDs from final output file names (e.g., ...ImpID5...)."""
        found_ids: set[int] = set()
        pattern = re.compile(r"ImpID(\d+)")

        for path in self.get_final_output_files():
            match = pattern.search(os.path.basename(path))
            if match:
                found_ids.add(int(match.group(1)))

        return found_ids
    
    
    
    def upload_number_of_variants(self) -> pd.DataFrame:
        """Count variants in imputed VCFs and upsert the results into snpster_users.imputed_data."""

        candidate_files = [
            vcf_path
            for result_dir in self.job_df["output_dir"].astype(str).tolist()
            for vcf_path in sorted((Path(result_dir) / "qc_output").glob("*.vcf.gz"))
        ]

        rows = []

        for vcf_path in candidate_files:
            number_of_variants = self.vcf_utilities.get_number_of_variants_in_vcf(str(vcf_path))

            rows.append({
                "imputation_id": int(self.vcf_utilities._get_imputation_id_from_vcf(str(vcf_path))),
                "file_path": str(vcf_path),
                "number_of_variants": int(number_of_variants),
                "chromosome": int(self.vcf_utilities._get_chromosome_from_vcf(str(vcf_path))),
            })

        upload_df = pd.DataFrame(rows, columns=["imputation_id", "file_path", "number_of_variants", "chromosome"])

        self.query_handler.db_utils.upsert_dataframe_to_db(
            dataframe=upload_df,
            table_name="imputed_data",
            schema="snpster_users",
            conflict_columns=["imputation_id", "chromosome"],
        )

        print(f"Upserted {len(upload_df)} rows into snpster_users.imputed_data")
        return upload_df
        
        
        
        
    
    


    
        
        
