import os
import subprocess
import pandas as pd
from pathlib import Path
import tempfile
import shutil
from config import (
    REFERENCE_VCF_DIR, REFERENCE_VCF_PATTERN, POPULATION_PANEL_FILE, 
    DEFAULT_CHROMOSOME, K_POPULATIONS, ANCESTRY_METHOD, REFERENCE_PANEL,
    ANCESTRY_LABELS, LEGACY_POPULATION_MAP, DB_ANCESTRY_COLUMNS
)
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH
from db_handler import DbHandler, DbUtils

class AncestryEnvironmentHandler:
    def __init__(
        self,
        vcf_file: str | None = None,
        bed_file: str | None = None,
        output_dir: str = ".",
        reference_vcf_dir: str = REFERENCE_VCF_DIR,  # Your existing 1000G HGDP VCFs
        reference_vcf_pattern: str = REFERENCE_VCF_PATTERN,
        population_panel_file: str = POPULATION_PANEL_FILE,  # Population assignments
        use_chromosomes: list = None,  # List of chromosomes to use (default: [1, 2, 21, 22])
        k_populations: int = K_POPULATIONS,  # Number of ancestral populations
        ancestry_results = None
    ):
        self.vcf_file = vcf_file
        self.bed_file = bed_file
        self.output_dir = output_dir
        self.reference_vcf_dir = reference_vcf_dir
        self.reference_vcf_pattern = reference_vcf_pattern
        self.population_panel_file = population_panel_file
        self.use_chromosomes = use_chromosomes if use_chromosomes else ['1', '2', '21', '22']
        self.k_populations = k_populations
        self.ancestry_results = ancestry_results
        os.makedirs(output_dir, exist_ok=True)


class AncestryInference:

    def __init__(
        self,
        environment_handler: AncestryEnvironmentHandler,
        db_utils: DbUtils | None = None
    ):
        self.env = environment_handler

        if db_utils is None:
            db_handler = DbHandler(
                port=PORT,
                user=USERNAME,
                password=PASSWORD,
                host=HOST
            )
            self.db_utils = DbUtils(db_handler=db_handler)
        else:
            self.db_utils = db_utils
        
        
    def _get_reference_vcf_for_chr(self, chrom: str) -> str:
        """Find the reference VCF file for a given chromosome."""
        pattern = self.env.reference_vcf_pattern.format(chrom=chrom)
        ref_path = os.path.join(self.env.reference_vcf_dir, pattern)
        
        if not os.path.exists(ref_path):
            raise FileNotFoundError(f"Reference VCF not found: {ref_path}")
        
        return ref_path
    
    def _extract_sample_chr_vcf(self, input_vcf: str, chrom: str, output_vcf: str) -> str:
        """Prepare sample VCF - recompress with bgzip if needed."""
        print(f"Preparing sample VCF for chromosome {chrom}...")
        
        temp_plain = None
        temp_idx = input_vcf + ".tbi"
        
        try:
            # Test if it's bgzip by trying to index it
            # bgzip files can be indexed, regular gzip cannot
            test_result = subprocess.run(
                ["bcftools", "index", "-t", "-f", input_vcf],
                capture_output=True,
                check=False
            )
            
            if test_result.returncode != 0:
                # File is not bgzip compressed, need to recompress
                print(f"  Recompressing with bgzip (file was gzip compressed)...")
                temp_plain = input_vcf.replace(".vcf.gz", ".temp.vcf")
                
                # Decompress to plain text
                with open(temp_plain, 'w') as f_out:
                    subprocess.run(["gunzip", "-c", input_vcf], 
                                 stdout=f_out, 
                                 check=True)
                
                # Recompress with bgzip to output
                subprocess.run([
                    "bcftools", "view",
                    "-Oz", "-o", output_vcf,
                    temp_plain
                ], check=True)
            else:
                # Already bgzip, just copy to output location
                import shutil
                shutil.copy(input_vcf, output_vcf)
            
            # Index the output file
            subprocess.run(["bcftools", "index", "-t", output_vcf], check=True)
            
        finally:
            # Clean up temp files if created
            if temp_plain and os.path.exists(temp_plain):
                os.remove(temp_plain)
            if temp_idx and os.path.exists(temp_idx):
                os.remove(temp_idx)
        
        return output_vcf
    
    def _merge_sample_with_reference(self, sample_vcf: str, ref_vcf: str, output_prefix: str) -> str:
        """Merge sample and reference VCFs using only exact matching variants."""
        print(f"Merging sample with 1000 Genomes reference...")
        
        merged_vcf = f"{output_prefix}.merged.vcf.gz"
        
        with tempfile.TemporaryDirectory(
            dir=self.env.output_dir,
            prefix="ancestry_intersection_"
        ) as intersection_dir:
            subprocess.run([
                "bcftools", "isec",
                "--nfiles", "=2",
                "--collapse", "none",
                "-Oz", "-p", intersection_dir,
                sample_vcf,
                ref_vcf
            ], check=True)

            subprocess.run([
                "bcftools", "merge",
                "--force-samples",  # Allow sample name conflicts
                os.path.join(intersection_dir, "0000.vcf.gz"),
                os.path.join(intersection_dir, "0001.vcf.gz"),
                "-Oz", "-o", merged_vcf
            ], check=True)
        
        subprocess.run(["bcftools", "index", "-t", merged_vcf], check=True)
        
        return merged_vcf
    
    def _convert_vcf_to_plink(self, vcf_file: str, output_prefix: str) -> str:
        """Convert VCF to PLINK format for ADMIXTURE."""
        print(f"Converting merged VCF to PLINK BED format...")
        
        subprocess.run([
            "plink2",
            "--vcf", vcf_file,
            "--make-bed",
            "--out", output_prefix,
            "--allow-extra-chr",
            "--max-alleles", "2"  # Keep only biallelic SNPs
        ], check=True)
        
        return f"{output_prefix}.bed"
    
    def _create_population_file(self, merged_bed_prefix: str, panel_file: str) -> str:
        """Create .pop file for ADMIXTURE supervised mode with dynamic population support."""
        print("Creating population assignment file...")
        print(f"Reading panel file from: {panel_file}")
        
        # Read sample IDs from .fam file
        fam_file = f"{merged_bed_prefix}.fam"
        fam_df = pd.read_csv(fam_file, sep=r'\s+', header=None, usecols=[0, 1])
        fam_df.columns = ['family_id', 'sample_id']
        
        # Read 1000G HGDP population panel
        panel_df = pd.read_csv(panel_file, sep='\t')
        # Expected columns: sample, pop, super_pop, gender
        print(f"Panel file loaded: {len(panel_df)} samples")
        print(f"Panel columns: {list(panel_df.columns)}")
        print(f"Sample super_pop values (first 20): {panel_df['super_pop'].head(20).tolist()}")
        
        # Get unique populations from panel and create dynamic mapping
        unique_pops = panel_df['super_pop'].unique()
        super_pop_map = {pop: pop for pop in unique_pops}
        
        # Apply legacy mappings (e.g., AFR → SSA)
        for old_code, new_code in LEGACY_POPULATION_MAP.items():
            if old_code in super_pop_map:
                super_pop_map[old_code] = new_code
        
        print(f"Found {len(unique_pops)} populations in reference panel: {sorted(unique_pops)}")
        
        # Create population assignments
        pop_assignments = []
        for _, row in fam_df.iterrows():
            sample_id = row['sample_id']
            
            # Check if this is a reference sample (in panel) or user sample
            if sample_id in panel_df['sample'].values:
                super_pop = panel_df[panel_df['sample'] == sample_id]['super_pop'].values[0]
                mapped_pop = super_pop_map.get(super_pop, super_pop)
                pop_assignments.append(mapped_pop)
            else:
                # User sample - marked as unknown for supervised learning
                pop_assignments.append('-')
        
        # Write .pop file
        pop_file = f"{merged_bed_prefix}.pop"
        with open(pop_file, 'w') as f:
            for pop in pop_assignments:
                f.write(f"{pop}\n")
        
        print(f"Population file created: {pop_file}")
        return pop_file
    
    def _run_admixture(self, bed_file: str, k: int = 5) -> str:
        """Run ADMIXTURE in supervised mode."""
        print(f"Running ADMIXTURE (K={k}, supervised mode)...")
        
        bed_path = Path(bed_file)
        
        # ADMIXTURE must be run from the directory containing the files
        original_dir = os.getcwd()
        os.chdir(bed_path.parent)
        
        try:
            subprocess.run([
                "admixture",
                "--supervised",
                "-j4",
                bed_path.name,
                str(k)
            ], check=True)
        finally:
            os.chdir(original_dir)
        
        # Output: {bed_prefix}.{K}.Q
        q_file = str(bed_path.with_suffix(f'.{k}.Q'))
        return q_file
    
    def _parse_admixture_results(self, q_file: str, fam_file: str, panel_file: str) -> pd.DataFrame:
        """Parse ADMIXTURE results with dynamic K population support."""
        
        # Read sample IDs
        fam_df = pd.read_csv(fam_file, sep=r'\s+', header=None, usecols=[0, 1])
        fam_df.columns = ['family_id', 'sample_id']
        
        # Read ancestry fractions
        q_df = pd.read_csv(q_file, sep=r'\s+', header=None)
        
        # Determine column order from .pop file
        # ADMIXTURE assigns columns based on order populations first appear in .pop file
        pop_file = q_file.rsplit('.', 2)[0] + '.pop'  # e.g., merged_plink.5.Q -> merged_plink.pop
        with open(pop_file, 'r') as f:
            pop_labels = [line.strip() for line in f if line.strip() != '-']
        
        # Get unique populations in order of first appearance
        seen = set()
        ancestry_cols = []
        for pop in pop_labels:
            if pop not in seen and pop != '-':
                seen.add(pop)
                ancestry_cols.append(pop)
        
        # Ensure we have exactly K columns
        ancestry_cols = ancestry_cols[:q_df.shape[1]]
        
        # Validate population codes
        for pop_code in ancestry_cols:
            if pop_code not in ANCESTRY_LABELS:
                print(f"Warning: Unknown population code '{pop_code}' - adding to results as-is")
        
        q_df.columns = ancestry_cols
        
        print(f"Detected {len(ancestry_cols)} populations in ADMIXTURE results: {ancestry_cols}")
        
        # Combine
        result_df = pd.concat([fam_df, q_df], axis=1)
        
        # Determine primary ancestry using the ancestry columns
        result_df['primary_ancestry'] = result_df[ancestry_cols].idxmax(axis=1)
        result_df['primary_ancestry_percentage'] = result_df[ancestry_cols].max(axis=1)
        
        return result_df




    def run_ancestry_inference(self) -> pd.DataFrame:
        """Main workflow: infer ancestry using multiple chromosomes from 1000G HGDP VCF references."""
        
        chromosomes = self.env.use_chromosomes
        merged_vcfs = []
        
        # Step 1-3: Process each chromosome separately
        for chrom in chromosomes:
            print(f"\n=== Processing chromosome {chrom} ===")
            
            # Get reference VCF for this chromosome
            ref_vcf = self._get_reference_vcf_for_chr(chrom)
            print(f"Using reference: {ref_vcf}")
            
            # Extract same chromosome from sample VCF
            sample_chr_vcf = os.path.join(self.env.output_dir, f"sample_chr{chrom}.vcf.gz")
            
            if self.env.vcf_file:
                self._extract_sample_chr_vcf(self.env.vcf_file, chrom, sample_chr_vcf)
            else:
                raise ValueError("vcf_file is required")
            
            # Merge sample with reference at common positions
            merged_prefix = os.path.join(self.env.output_dir, f"merged_chr{chrom}")
            merged_vcf = self._merge_sample_with_reference(sample_chr_vcf, ref_vcf, merged_prefix)
            merged_vcfs.append(merged_vcf)
        
        # Step 4: Concatenate all chromosome VCFs
        print(f"\n=== Concatenating {len(merged_vcfs)} chromosome VCFs ===")
        combined_vcf = os.path.join(self.env.output_dir, "merged_all_chrs.vcf.gz")
        subprocess.run([
            "bcftools", "concat",
            "-Oz", "-o", combined_vcf,
            *merged_vcfs
        ], check=True)
        subprocess.run(["bcftools", "index", "-t", combined_vcf], check=True)
        
        # Step 5: Convert to PLINK format
        plink_prefix = os.path.join(self.env.output_dir, "merged_plink")
        self._convert_vcf_to_plink(combined_vcf, plink_prefix)
        
        # Step 6: Create population file for supervised ADMIXTURE
        self._create_population_file(plink_prefix, self.env.population_panel_file)
        
        # Step 7: Run ADMIXTURE
        bed_file = f"{plink_prefix}.bed"
        q_file = self._run_admixture(bed_file, k=self.env.k_populations)
        
        # Step 8: Parse results (user samples only)
        fam_file = f"{plink_prefix}.fam"
        results_df = self._parse_admixture_results(q_file, fam_file, self.env.population_panel_file)
        
        self.env.ancestry_results = results_df
    
    
    def upload_ancestry_results(self):
        """Upload ancestry results to database or save to file.
            target table:
            CREATE TABLE snpster_users.user_ancestry (
                user_id VARCHAR(100) PRIMARY KEY
                    REFERENCES snpster_users.user_information(user_id)
                    ON DELETE CASCADE,

                eur NUMERIC(8,6) CHECK (eur >= 0 AND eur <= 1),
                afr NUMERIC(8,6) CHECK (afr >= 0 AND afr <= 1),
                eas NUMERIC(8,6) CHECK (eas >= 0 AND eas <= 1),
                sas NUMERIC(8,6) CHECK (sas >= 0 AND sas <= 1),
                amr NUMERIC(8,6) CHECK (amr >= 0 AND amr <= 1),

                primary_ancestry VARCHAR(10)
                    CHECK (primary_ancestry IN ('EUR', 'AFR', 'EAS', 'SAS', 'AMR')),

                primary_ancestry_percentage NUMERIC(8,6)
                    CHECK (primary_ancestry_percentage >= 0 AND primary_ancestry_percentage <= 1),

                ancestry_method VARCHAR(100),
                reference_panel VARCHAR(100),

                created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
            );
            
            output file example: 
            
            family_id,sample_id,EUR,AFR,EAS,SAS,AMR,primary_ancestry,primary_ancestry_percentage
            0,IMPID_5,0.827631,1e-05,0.009602,0.081771,0.080986,EUR,0.827631
            0,IMPID_1,0.746731,1e-05,1e-05,0.147048,0.106201,EUR,0.746731
            0,IMPID_3,0.767679,0.004482,1e-05,0.1187,0.109129,EUR,0.767679
            0,IMPID_2,0.809307,0.011837,1e-05,0.069393,0.109452,EUR,0.809307
            0,IMPID_4,0.783832,0.023304,0.01798,0.083357,0.091527,EUR,0.783832
            0,HG00405,1.4e-05,0.999956,1e-05,1e-05,1e-05,AFR,0.999956
            0,HG00408,1e-05,0.99996,1e-05,1e-05,1e-05,AFR,0.99996
            0,HG00418,1e-05,0.99996,1e-05,1e-05,1e-05,AFR,0.99996

        """
        
        #modify the df to fit the requirements for the db table, e.g. rename columns, filter out non IMPID samples, etc.
        
        #prune non IMPID samples
        results_df = self.env.ancestry_results[self.env.ancestry_results['sample_id'].str.startswith('IMPID')].copy()
        
        #remove IMPID_ from string
        results_df.loc[:, 'sample_id'] = results_df['sample_id'].str.replace('IMPID_', '', regex=False)
        
        #convert to integer to match database type
        results_df.loc[:, 'sample_id'] = results_df['sample_id'].astype(int)
        
        #query db to get user_id
        
        imputation_ids = results_df['sample_id'].tolist()
        
        # Try to connect with fast fail (2 retries, 5 second wait = max 10 seconds)
        if not self.db_utils.db_handler.connect(retries=2, wait_time=5):
            raise RuntimeError("Failed to connect to database for ancestry upload. Check DB_HOST, DB_PORT, and network connectivity.")
        
        try:
            ids = self.db_utils.get_user_id_from_imputation_id(imputation_ids)
            
            if ids is None or ids.empty:
                raise ValueError(f"No matching imputation_ids found in database: {imputation_ids}")
            
            #join ids on results_df by the ids column imputation_id and sample_id
            results_df = results_df.merge(ids, left_on='sample_id', right_on='imputation_id', how='left')
            
            #drop sample_id, imputation_id, and family_id columns
            results_df = results_df.drop(columns=['sample_id', 'imputation_id', 'family_id'])
            
            #make all colnames lowercase
            results_df.columns = [col.lower() for col in results_df.columns]
            
            # Ensure all expected database columns exist (set to NULL if missing)
            for col in DB_ANCESTRY_COLUMNS:
                if col not in results_df.columns:
                    results_df[col] = None
                    print(f"Note: Population '{col.upper()}' not in inference results - setting to NULL")
            
            # Add metadata columns
            results_df['ancestry_method'] = ANCESTRY_METHOD
            results_df['reference_panel'] = REFERENCE_PANEL
            
            # Select only required columns in correct order
            output_cols = ['user_id'] + DB_ANCESTRY_COLUMNS + ['primary_ancestry', 'primary_ancestry_percentage', 'ancestry_method', 'reference_panel']
            results_df = results_df[output_cols]
            
            #upload
            print("\nAncestry results to upload:")
            print(results_df)
            
            self.db_utils.upsert_dataframe_to_db(
                results_df, 
                table_name="user_ancestry", 
                schema="snpster_users",
                conflict_columns=['user_id']
            )
            print(f"✓ Successfully uploaded {len(results_df)} ancestry results to database")
            print(f"  Populations included: {[col.upper() for col in DB_ANCESTRY_COLUMNS if results_df[col].notna().any()]}")
            
        finally:
            self.db_utils.db_handler.close()
