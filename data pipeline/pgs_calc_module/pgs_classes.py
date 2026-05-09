import subprocess
import pandas as pd
import os
import shutil
import glob
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT
from config import SCORING_FILE_SOURCE_DIR, SCORING_FILE_TARGET_DIR, NF_WORK_DIR, OUTPUT_DIR, REFERENCE_DATA_PATH, BCFTOOLS_THREADS, SAMPLESET_NAME


class EnvironmentHandler:
    def __init__(self, 
                 merged_sample_sheet: str = None,
                 scoring_file_str: str = None,
                 user_ids: list = None,
                 imputation_ids: list = None,
                 prsc_ids: list = None,
                 id_map: pd.DataFrame = None,
                 db_utils: DbUtils = None,
                 samplesheet_paths: list = None,
                 combined_samplesheet_path: str = None,
                 sampleset_name: str = SAMPLESET_NAME,
                 scoring_file_source_dir: str = SCORING_FILE_SOURCE_DIR,
                 scoring_file_target_dir: str = SCORING_FILE_TARGET_DIR,
                 nf_work_dir: str = NF_WORK_DIR,
                 output_dir: str = OUTPUT_DIR,
                 reference_data_path: str = REFERENCE_DATA_PATH
                 ):
        
        self.merged_sample_sheet = merged_sample_sheet
        self.samplesheet_paths = samplesheet_paths
        self.combined_samplesheet_path = combined_samplesheet_path
        self.user_ids = user_ids
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.scoring_file_str = scoring_file_str
        self.imputation_ids = imputation_ids
        self.prsc_ids = prsc_ids
        self.id_map = id_map
        if db_utils is None:
            db_utils = DbUtils(DbHandler(user=USERNAME, password=PASSWORD, host=HOST, port=PORT))
        self.db_utils = db_utils
        self.scoring_file_source_dir = scoring_file_source_dir
        self.scoring_file_target_dir = scoring_file_target_dir
        self.nf_work_dir = nf_work_dir
        self.sampleset_name = sampleset_name
        
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

        update_query = f"""UPDATE snpster_users.prsc_jobs
                        SET prsc_status = '{status}'
                        WHERE prsc_id in ({', '.join(map(str, self.prsc_ids))});"""
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
            self.nf_work_dir
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
        self.environment_handler.id_map = df_subset
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

    

    def upload_results(self, path) -> None:
        
        """Uploads the PGS calculation results to the database and updates job status."""

        results = pd.read_csv(path, sep="\t")
        
        results["PGS"] = results["PGS"].str.replace("_hmPOS_GRCh38", "", regex=False)
        results = results[results["sampleset"] != "reference"]
        
        #group dataframe by sampleset
        
        results = results.groupby("sampleset")
        
        for sampleset, group in results:
            
            for index, row in group.iterrows():
                imputation_id = int(row['FID'])
                
                matched_rows = self.environment_handler.id_map[
                    self.environment_handler.id_map["imputation_id"] == imputation_id
                ]
                
                prsc_id = matched_rows["prsc_id"].values[0]
                pgs_id = row["PGS"]
                percentile_most_similar_pop = float(row["percentile_MostSimilarPop"])
                z_most_similar_pop = float(row["Z_MostSimilarPop"])
                z_norm1 = float(row["Z_norm1"])
                z_norm2 = float(row["Z_norm2"])
                score_sum = float(row["SUM"])
                percent_variants_matched = 'NULL'


                insert_query = f"""INSERT INTO snpster_users.prsc_job_results (prsc_id, pgs_id, percentile_most_similar_pop, z_norm1, z_norm2, z_most_similar_pop, score_sum, percent_variants_matched)
                                VALUES ({prsc_id}, '{pgs_id}', {percentile_most_similar_pop}, {z_norm1}, {z_norm2}, {z_most_similar_pop}, {score_sum}, {percent_variants_matched});"""
                
                self.environment_handler.db_utils.db_handler.execute_query(insert_query)
            
    def validate_results(self, path) -> bool:
        """Validates the PGS calculation results before uploading."""
        
        
        try:
            results = pd.read_csv(path, sep="\t")
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




class VCFHandler:
    def __init__(
        self,
        environment_handler: EnvironmentHandler,
        user_imputation_id_dict: dict[str, str] | None = None,
        full_vcf_file_path_dict: dict[str, list[str]] | None = None,
    ):
        self.environment_handler = environment_handler
        self.user_imputation_id_dict = user_imputation_id_dict
        self.full_vcf_file_path_dict = full_vcf_file_path_dict

    def set_user_imputation_id_dict(self) -> None:
        # Set mapping: imputation_id -> user_id from shared environment state.
        user_ids = self.environment_handler.user_ids
        imputation_ids = self.environment_handler.imputation_ids

        self.user_imputation_id_dict = {
            str(imputation_id): str(user_id)
            for imputation_id, user_id in zip(imputation_ids, user_ids)
        }

    def set_vcf_file_dict(self) -> None:
        # Build: chromosome -> list of full VCF paths across all selected jobs.

        sample_sheet_pd = pd.read_csv(self.environment_handler.combined_samplesheet_path)
        self.full_vcf_file_path_dict = {}

        for row in sample_sheet_pd.itertuples(index=False):
            chrom = str(row.chrom)
            path_prefix = str(row.path_prefix)
            fmt = str(row.format)
            imputation_id = str(row.sampleset)

            user_id = self.user_imputation_id_dict.get(imputation_id)
            if user_id is None:
                raise KeyError(f"No user_id mapping found for imputation_id={imputation_id}")

            file_name = f"{path_prefix}.{fmt}"
            full_path = f"/srv/imputed/{user_id}/{imputation_id}/output/{file_name}"

            self.full_vcf_file_path_dict.setdefault(chrom, []).append(full_path)

    def _prepare_merge_input(self, input_path: str) -> str:
        prepared_dir = os.path.join(self.environment_handler.nf_work_dir, "merge_inputs")
        os.makedirs(prepared_dir, exist_ok=True)

        if input_path.endswith(".vcf.gz"):
            if not os.path.exists(f"{input_path}.tbi"):
                subprocess.run(["bcftools", "index", "-t", input_path], check=True)
            return input_path

        file_name = os.path.basename(input_path)
        prepared_path = os.path.join(prepared_dir, f"{file_name}.gz")

        if not os.path.exists(prepared_path):
            subprocess.run(["bcftools", "view", input_path, "-Oz", "-o", prepared_path], check=True)

        if not os.path.exists(f"{prepared_path}.tbi"):
            subprocess.run(["bcftools", "index", "-t", prepared_path], check=True)

        return prepared_path

    def merge_vcf_files(self) -> None:

        """sampleset,path_prefix,chrom,format"""
        merged_file_sample_sheet = pd.DataFrame(columns=["sampleset", "path_prefix", "chrom", "format"])

        for chrom, files in self.full_vcf_file_path_dict.items():

            print(f"VCF files for chromosome {chrom}: {files}")
            output_path = f"{self.environment_handler.nf_work_dir}/cohort_chr{chrom}.vcf.gz"

            missing = [path for path in files if not os.path.exists(path)]
            if missing:
                raise FileNotFoundError(f"Missing VCF files for chromosome {chrom}: {missing}")

            prepared_files = [self._prepare_merge_input(path) for path in files]

            
            print(f"Merging {len(prepared_files)} files for chromosome {chrom} into {output_path}")


            subprocess.run(["bcftools", "merge", *prepared_files, "-Oz", "-o", output_path, '--threads', str(BCFTOOLS_THREADS)], check=True)
            
            base_filename = f"cohort_chr{chrom}"
            merged_file_sample_sheet.loc[len(merged_file_sample_sheet)] = {
                "sampleset": self.environment_handler.sampleset_name,
                "path_prefix": base_filename,
                "chrom": chrom,
                "format": "vcf"
            }
            
        merged_file_sample_sheet.to_csv(f"{self.environment_handler.nf_work_dir}/merged_samplesheet.csv", index=False)
        self.environment_handler.merged_sample_sheet = f"{self.environment_handler.nf_work_dir}/merged_samplesheet.csv"

        