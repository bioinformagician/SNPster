import subprocess
import pandas as pd
import os
import shutil
import glob
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, HOST, PORT
from config import PGS_CHUNK_SIZE, SCORING_FILE_SOURCE_DIR, SCORING_FILE_TARGET_DIR, NF_WORK_DIR, OUTPUT_DIR, REFERENCE_DATA_PATH, BCFTOOLS_THREADS, SAMPLESET_NAME, VCF_MERGE_SHEET_DIR, PGS_RESULT_DIR
from vcf_classes import VCFUtilities


vcf_utilities = VCFUtilities()

class EnvironmentHandler:
    def __init__(self, 
                 merged_sample_sheet: str = None,
                 scoring_file_strs: list[list] = None,
                 user_ids: list = None,
                 imputation_ids: list = None,
                 prsc_ids: list = None,
                 id_map: pd.DataFrame = None,
                 db_utils: DbUtils = None,
                 samplesheet_paths: list = None,
                 vcf_merge_sheet_dir: str = VCF_MERGE_SHEET_DIR,
                 vcf_merge_sheet: str = None,
                 sampleset_name: str = SAMPLESET_NAME,
                 scoring_file_source_dir: str = SCORING_FILE_SOURCE_DIR,
                 scoring_file_target_dir: str = SCORING_FILE_TARGET_DIR,
                 nf_work_dir: str = NF_WORK_DIR,
                 output_dir: str = OUTPUT_DIR,
                 reference_data_path: str = REFERENCE_DATA_PATH,
                 pgs_chunk_size: int = PGS_CHUNK_SIZE,
                 pgs_result_dir: str = PGS_RESULT_DIR
                 ):
        
        self.merged_sample_sheet = merged_sample_sheet
        self.samplesheet_paths = samplesheet_paths if samplesheet_paths is not None else []
        self.vcf_merge_sheet = vcf_merge_sheet
        self.vcf_merge_sheet_dir = vcf_merge_sheet_dir
        self.user_ids = user_ids
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.scoring_file_strs = scoring_file_strs
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
        self.pgs_chunk_size = pgs_chunk_size
        self.pgs_result_dir = pgs_result_dir
        
    def copy_scoring_files(self, target_dir, scoring_file_list: list) -> None:
        """Copies scoring files from the mounted source directory to the scoring file directory."""
    
        for file in scoring_file_list:
            shutil.copy2(file, target_dir)

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

    
    
    def clear_directories(self) -> None:
        """Clears the output and scoring file directories."""
        directories_to_clear = [
            self.output_dir,
            self.scoring_file_target_dir,
            self.nf_work_dir,
            self.vcf_merge_sheet_dir
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
        
        for pattern in ["/app/.nextflow*", "/app/.nextflow.log*"]:
            for path in glob.glob(pattern):
                if os.path.isfile(path) or os.path.islink(path):
                    os.unlink(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
    
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
    
    def create_samplesheet_for_directory(self, directory: str) -> None:
        """Creates a sample sheet for the VCF merging by listing all VCF files in the specified directory.

            frederik@frederik-MS-02-Ultra:/srv/imputed/0/1/qc_output$ ls
            IMPID1.chr10.split_QC.vcf.gz  IMPID1.chr16.split_QC.vcf.gz  IMPID1.chr21.split_QC.vcf.gz  IMPID1.chr6.split_QC.vcf.gz
            IMPID1.chr11.split_QC.vcf.gz  IMPID1.chr17.split_QC.vcf.gz  IMPID1.chr22.split_QC.vcf.gz  IMPID1.chr7.split_QC.vcf.gz
            IMPID1.chr12.split_QC.vcf.gz  IMPID1.chr18.split_QC.vcf.gz  IMPID1.chr2.split_QC.vcf.gz   IMPID1.chr8.split_QC.vcf.gz
            IMPID1.chr13.split_QC.vcf.gz  IMPID1.chr19.split_QC.vcf.gz  IMPID1.chr3.split_QC.vcf.gz   IMPID1.chr9.split_QC.vcf.gz
            IMPID1.chr14.split_QC.vcf.gz  IMPID1.chr1.split_QC.vcf.gz   IMPID1.chr4.split_QC.vcf.gz
            IMPID1.chr15.split_QC.vcf.gz  IMPID1.chr20.split_QC.vcf.gz  IMPID1.chr5.split_QC.vcf.gz
        
        """
        
        
        
        rows = []
        
        for entry in os.scandir(directory):
            
            base_path = os.path.dirname(entry.path)
            
            if entry.is_file() and entry.name.endswith(".vcf.gz"):
                full_path = os.path.join(base_path, entry.name)
                
                # Extract imputation ID from filename pattern: IMPID{id}.chr{chrom}.split_QC.vcf.gz
                if entry.name.startswith("IMPID") and ".chr" in entry.name:
                    imputation_id = entry.name.split(".")[0].replace("IMPID", "")
                    chrom = entry.name.split(".chr")[1].split(".")[0]
                else:
                    # Fallback to reading from VCF if filename doesn't match expected pattern
                    imputation_id = vcf_utilities._get_imputation_id_from_vcf(full_path)
                    chrom = vcf_utilities._get_chromosome_from_vcf(full_path)
                
                path_prefix = os.path.join(base_path, entry.name.split(".chr")[0])
                format = "vcf" #doesnt accept the compression extension
                
                rows.append({
                    "sampleset": imputation_id, #imputation_id
                    "path_prefix": path_prefix, #/path/to/target_genomes/cineca_synthetic_subset
                    "chrom": chrom, #22,21,20....
                    "format": format #vcf
                })
        
        samplesheet = pd.DataFrame(rows, columns=["sampleset", "path_prefix", "chrom", "format"])
        # Create unique filename using directory hash to avoid collisions
        samplesheet_path = os.path.join(self.nf_work_dir, f"samplesheet_IMPID{imputation_id}.csv")
        samplesheet.to_csv(samplesheet_path, index=False)
        self.samplesheet_paths.append(samplesheet_path)
    
    
    def create_full_path_samplesheets(self) -> None:
        """
        Creates a sample sheet for the VCF merging by 
        adding the full path to the prefix column of the samplesheets
        
        This is needed for the vcf_combiner.py script executed by nextflow to find the files to merge
        """
        
        combined_df = pd.DataFrame(columns=["sampleset", "full_vcf_path", "chrom", "format"])
        
        for samplesheet_path in self.samplesheet_paths:
            df = pd.read_csv(samplesheet_path)
            
            # path_prefix already contains the full path + base name, just need to add file extension
            df["full_vcf_path"] = df.apply(lambda row: f"{row['path_prefix']}.chr{row['chrom']}.split_QC.vcf.gz", axis=1)
            combined_df = pd.concat([combined_df, df], ignore_index=True)
        
        #write one file per chromosome
        
        for chrom in combined_df["chrom"].unique():
            chrom_df = combined_df[combined_df["chrom"] == chrom]
            output_path = f"{self.vcf_merge_sheet_dir}/vcf_merge_sheet_chr{chrom}.csv"
            chrom_df.to_csv(output_path, index=False)
            print(f"Created VCF merge samplesheet: {output_path}")
    
    
    def create_merged_vcf_samplesheet(self) -> None:
        """ create the new samplesheet with merged vcf files and set the environment variable"""
        
        rows = []
        
        for filename in os.listdir(self.vcf_merge_sheet_dir):
            if not filename.endswith(".vcf.gz"):
                continue

            chrom = None
            if filename.startswith("merged_vcf_chr"):
                # Legacy filename format: merged_vcf_chr<chrom>.vcf.gz
                chrom = filename[len("merged_vcf_chr"):].replace(".vcf.gz", "")
            elif filename.endswith(".merged.vcf.gz"):
                # Current filename format from vcf_classes: <chrom>.merged.vcf.gz
                chrom = filename.replace(".merged.vcf.gz", "")

            if chrom is not None:
                chrom = chrom.replace("chr", "")
                sampleset = self.sampleset_name
                path_prefix = os.path.join(self.vcf_merge_sheet_dir, filename.replace(".vcf.gz", ""))
                format = "vcf" #doesnt accept the compression extension
                
                rows.append({
                    "sampleset": sampleset,
                    "path_prefix": path_prefix,
                    "chrom": chrom,
                    "format": format
                })
        
        samplesheet = pd.DataFrame(rows, columns=["sampleset", "path_prefix", "chrom", "format"])
        samplesheet.to_csv(f"{self.vcf_merge_sheet_dir}/merged_sample_sheet.csv", index=False)
        self.merged_sample_sheet = f"{self.vcf_merge_sheet_dir}/merged_sample_sheet.csv"

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
                        SELECT DISTINCT
                            pj.imputation_id,
                            pj.prsc_id,
                            uf.user_id,
                            pj.prsc_status,
                            ij.imputation_status
                        FROM snpster_users.prsc_jobs pj
                        JOIN snpster_users.imputation_jobs ij
                            ON pj.imputation_id = ij.imputation_id
                        JOIN snpster_users.imputation_job_parameters ijp
                            ON ij.imputation_id = ijp.imputation_id
                        JOIN snpster_users.user_files uf
                            ON uf.file_id = ijp.file_id
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
                    ),
                    limited_matching_jobs AS (
                        SELECT prsc_id
                        FROM matching_jobs
                        ORDER BY prsc_id
                        LIMIT 30
                    )
                    SELECT
                        e.imputation_id,
                        e.prsc_id,
                        e.user_id,
                        e.prsc_status,
                        pjp.pgs_id,
                        e.imputation_status,
                        prs.scoring_file_path
                    FROM limited_matching_jobs mj
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
        
        self.environment_handler.base_folder_paths = [
            f"/srv/imputed/{row.user_id}/{row.imputation_id}/qc_output/"
            for row in df_subset.itertuples(index=False)
        ]
        self.environment_handler.samplesheet_paths = []  # Initialize empty, will be populated by create_samplesheets()
        
        scoring_files = sorted(set(results["scoring_file_path"].tolist()))
        
        #divide into chunks of 50 to avoid memory issues
        
        scoring_files_chunks = [scoring_files[i:i + self.environment_handler.pgs_chunk_size] for i in range(0, len(scoring_files), self.environment_handler.pgs_chunk_size)]
        
        #create subdir for each chunk
        
        subdirs = [os.path.join(self.environment_handler.scoring_file_target_dir, f"chunk_{i}") for i in range(len(scoring_files_chunks))]
        
        for subdir in subdirs:
            os.makedirs(subdir, exist_ok=True)
        
        for subdir, chunk in zip(subdirs, scoring_files_chunks):
            self.environment_handler.copy_scoring_files(subdir, chunk)

        self.environment_handler.scoring_file_strs = [subdir + "/*.txt.gz" for subdir in subdirs]

        return True
    
    def copy_imputed_files(self) -> None:
        """Copies imputed files from the mounted imputed directory flat into nf_work_dir."""
        for idx, imputation_id in enumerate(self.environment_handler.imputation_ids):
            source_dir = f"/srv/imputed/{self.environment_handler.user_ids[idx]}/{imputation_id}/qc_output/"
            for entry in os.scandir(source_dir):
                if entry.is_file() and not entry.name.endswith(".csv"):
                    print(f"Copying {entry.path} to {self.environment_handler.nf_work_dir}")
                    shutil.copy2(entry.path, self.environment_handler.nf_work_dir)
    
    
    
    def create_chr_specific_samplesheet(self) -> None:
        """Combines all samplesheets from the required imputation jobs 
        into a single samplesheet for the NF pipeline. first gather all
        and then split them"""

        samplesheet_dfs = []
        for samplesheet_path in self.environment_handler.samplesheet_paths:
            print(f"Reading samplesheet from {samplesheet_path}")
            df = pd.read_csv(samplesheet_path)
            samplesheet_dfs.append(df)

        combined_df = pd.concat(samplesheet_dfs, ignore_index=True)
        
        for chrom, group_df in combined_df.groupby("chrom"):
            #split df by chromosome and write one file per chromosome
            group_df.to_csv(os.path.join(self.environment_handler.vcf_merge_sheet_dir, f"samplesheet_chr{chrom}.csv"), index=False)
        
    
    
    def move_pgs_results(self, scoring_file:str, summary_file:str) -> None:
        """Moves the PGS calculation results to the designated result directory.
           This is needed for batched runs so the full data can be uploaded to the db
           after all batches are compeleted, so jobs are not half finished so the server
           can just be shut down without issues"""
           
        if not os.path.exists(self.environment_handler.pgs_result_dir):
            os.makedirs(self.environment_handler.pgs_result_dir, exist_ok=True)
        
        #if dir is empty, move the files, else append to existing files
        
        if os.listdir(self.environment_handler.pgs_result_dir):
            # Append to existing files
            existing_score_file = os.path.join(self.environment_handler.pgs_result_dir, os.path.basename(scoring_file))
            existing_summary_file = os.path.join(self.environment_handler.pgs_result_dir, os.path.basename(summary_file))
            
            existing_score_file_pd = pd.read_csv(existing_score_file, sep="\t")
            existing_summary_file_pd = pd.read_csv(existing_summary_file)
            
            new_score_file_pd = pd.read_csv(scoring_file, sep="\t")
            new_summary_file_pd = pd.read_csv(summary_file)
            
            # Append new data to existing data
            combined_score_file_pd = pd.concat([existing_score_file_pd, new_score_file_pd], ignore_index=True)
            combined_summary_file_pd = pd.concat([existing_summary_file_pd, new_summary_file_pd], ignore_index=True)
            
            # Write combined data back to the existing files
            combined_score_file_pd.to_csv(existing_score_file, sep="\t", index=False)
            combined_summary_file_pd.to_csv(existing_summary_file, index=False)
        
        else:
            # Move new files to the result directory
            shutil.move(scoring_file, os.path.join(self.environment_handler.pgs_result_dir, os.path.basename(scoring_file)))
            shutil.move(summary_file, os.path.join(self.environment_handler.pgs_result_dir, os.path.basename(summary_file)))
            
        
        
        
    def upload_results(self, scoring_file:str, summary_file:str) -> None:
        
        """Uploads the PGS calculation results to the database and updates job status.
        This has been enshittified by copilot, make clean version later"""

        results = pd.read_csv(scoring_file, sep="\t")
        summary_statistics = pd.read_csv(summary_file)
        
        # Debug: Print available columns
        print(f"Summary file columns: {list(summary_statistics.columns)}")
        print(f"Results file columns: {list(results.columns)}")
        print(f"Summary file shape: {summary_statistics.shape}")
        print(f"First few rows of summary:\n{summary_statistics.head()}")
        
        # Filter summary statistics to only include matched scores with proper IDs
        if "match_status" in summary_statistics.columns and "match_IDs" in summary_statistics.columns:
            # Convert match_IDs to string for comparison (could be bool or string)
            summary_statistics["match_IDs"] = summary_statistics["match_IDs"].astype(str)
            summary_statistics = summary_statistics[
                (summary_statistics["match_status"] == "matched") & 
                (summary_statistics["match_IDs"].str.lower() == "true")
            ]
            print(f"Filtered summary to {len(summary_statistics)} matched score rows (before aggregation)")
        
        # Clean PGS names in summary for matching
        if "accession" in summary_statistics.columns:
            summary_statistics["accession_clean"] = summary_statistics["accession"].str.replace("_hmPOS_GRCh38", "", regex=False)
        
        # Aggregate by dataset and accession to combine match_flipped=true and match_flipped=false rows
        # Sum counts and average percentages weighted by counts
        if "count" in summary_statistics.columns and "percent" in summary_statistics.columns:
            summary_agg = summary_statistics.groupby(["dataset", "accession_clean"]).agg({
                "count": "sum",
                "percent": lambda x: (x * summary_statistics.loc[x.index, "count"]).sum() / summary_statistics.loc[x.index, "count"].sum() if summary_statistics.loc[x.index, "count"].sum() > 0 else 0
            }).reset_index()
            print(f"Aggregated to {len(summary_agg)} unique PGS scores per dataset")
            summary_statistics = summary_agg
        
        # Remove reference samples
        results = results[results["sampleset"] != "reference"]
        
        # Clean PGS names in results BEFORE merging
        results["PGS_clean"] = results["PGS"].str.replace("_hmPOS_GRCh38", "", regex=False)
        
        # Debug: Show what we're trying to match
        print(f"Sample PGS values in results: {results['PGS'].head(3).tolist()}")
        print(f"Sample PGS_clean values: {results['PGS_clean'].head(3).tolist()}")
        print(f"Sample accession_clean values: {summary_statistics['accession_clean'].head(3).tolist()}")
        print(f"Unique samplesets in results: {results['sampleset'].unique()}")
        print(f"Unique datasets in summary: {summary_statistics['dataset'].unique()}")
        
        # Merge results with summary statistics
        if "dataset" in summary_statistics.columns and "accession_clean" in summary_statistics.columns:
            results = results.merge(
                summary_statistics, 
                left_on=["sampleset", "PGS_clean"], 
                right_on=["dataset", "accession_clean"], 
                how="inner"  # Changed to inner join to only keep matched scores
            )
            print(f"After merge: {len(results)} result rows")
        else:
            raise ValueError("Summary file missing expected columns 'dataset' or 'accession'")
        
        # Drop rows with NaN in critical columns
        results = results.dropna(subset=["count", "percent"])

        # Build a fast lookup for imputation_id -> prsc_id (keep first mapping, same behavior as before)
        id_map_df = self.environment_handler.id_map[["imputation_id", "prsc_id"]].drop_duplicates(subset=["imputation_id"])
        imputation_to_prsc = dict(zip(id_map_df["imputation_id"], id_map_df["prsc_id"]))

        insert_rows = []
        missing_imputation_ids = set()

        for _, row in results.iterrows():
            # Extract imputation ID from FID
            imputation_id = int(str(row["FID"]).split("_")[-1])
            prsc_id = imputation_to_prsc.get(imputation_id)

            if prsc_id is None:
                missing_imputation_ids.add(imputation_id)
                continue

            insert_rows.append({
                "prsc_id": int(prsc_id),
                "pgs_id": row["PGS_clean"],
                "percentile_most_similar_pop": float(row["percentile_MostSimilarPop"]),
                "z_norm1": float(row["Z_norm1"]),
                "z_norm2": float(row["Z_norm2"]),
                "z_most_similar_pop": float(row["Z_MostSimilarPop"]),
                "score_sum": float(row["SUM"]),
                "percent_variants_matched": float(row["percent"]),
                "n_variants_matched": int(row["count"]),
            })

        if missing_imputation_ids:
            available_ids = sorted(self.environment_handler.id_map["imputation_id"].drop_duplicates().tolist())
            raise ValueError(
                "No matching prsc_id found for imputation_id(s): "
                f"{sorted(missing_imputation_ids)}. Available imputation_ids in id_map: {available_ids}"
            )

        if not insert_rows:
            print("No rows to insert into snpster_users.prsc_job_results after filtering/merging.")
            return

        insert_df = pd.DataFrame(insert_rows)
        self.environment_handler.db_utils.insert_dataframe_to_db(
            dataframe=insert_df,
            table_name="prsc_job_results",
            schema="snpster_users",
        )
        print(f"Inserted {len(insert_df)} rows into snpster_users.prsc_job_results.")
            
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
    
    
    def create_samplesheets(self) -> None:
        
        for idx, imputation_id in enumerate(self.environment_handler.imputation_ids):
            imputed_dir = f"/srv/imputed/{self.environment_handler.user_ids[idx]}/{imputation_id}/qc_output"
            self.environment_handler.create_samplesheet_for_directory(imputed_dir)
    
    
    def run_vcf_merging(self) -> None:
        """nextflow run vcf_combiner_pipeline.nf --vcf_samplesheet_dir ./nf_data --output_dir ./nf_data"""
        
        print("Running VCF merging with Nextflow...")
        print(f"Looking for files in: {self.environment_handler.vcf_merge_sheet_dir}")
        
        # Debug: list files in the directory
        files_in_dir = os.listdir(self.environment_handler.vcf_merge_sheet_dir)
        merge_sheets = [f for f in files_in_dir if f.startswith("vcf_merge_sheet_chr") and f.endswith(".csv")]
        print(f"Found {len(merge_sheets)} merge sheet files: {merge_sheets}")
        
        command = [
            "nextflow", "run", "/app/vcf_combiner_pipeline.nf",
            "-work-dir", self.environment_handler.nf_work_dir,
            "--vcf_samplesheet_dir", self.environment_handler.vcf_merge_sheet_dir,
            "--output_dir", self.environment_handler.vcf_merge_sheet_dir,
        ]
        
        print(f"running command {command}")
        
        subprocess.run(command, check=True)

    
    def run_pgs_calculation(self, sample_sheet: str, scoring_file_str) -> None:
        if not sample_sheet or not os.path.exists(sample_sheet):
            raise FileNotFoundError(f"PGS input samplesheet not found: {sample_sheet}")

        input_df = pd.read_csv(sample_sheet, dtype=str)
        if input_df.empty:
            raise ValueError(f"PGS input samplesheet has no data rows: {sample_sheet}")

        if not self.environment_handler.reference_data_path or not os.path.exists(self.environment_handler.reference_data_path):
            raise FileNotFoundError(
                "PGS reference data archive not found at "
                f"{self.environment_handler.reference_data_path}. "
                "Check the pgs_calc volume mount and REFERENCE_DATA_PATH."
            )
    
        command = [
                "nextflow", "run", "/opt/pgsc_calc/main.nf",
                "-work-dir", self.environment_handler.nf_work_dir,
                "-profile", "conda",
                "--input", sample_sheet,
                "--target_build", self.pgscalculator_config.target_build,
                "--run_ancestry", self.environment_handler.reference_data_path,
                "--outdir", self.environment_handler.output_dir,
                "--min_overlap", "0.01",
                "--scorefile", scoring_file_str,
            ]
        
        print(f"Running PGS calculation with command: {' '.join(command)}")

        subprocess.run(command, check=True)
        print("PGS calculation completed successfully.")
            



