import os
import subprocess
import pandas as pd
from pathlib import Path
import tempfile
import shutil
from config import (
    REFERENCE_VCF_DIR, REFERENCE_VCF_PATTERN, POPULATION_PANEL_FILE, 
    DEFAULT_CHROMOSOME, K_POPULATIONS, ANCESTRY_METHOD, REFERENCE_PANEL,
    ANCESTRY_LABELS, LEGACY_POPULATION_MAP, DB_ANCESTRY_COLUMNS, MIN_ANCESTRY_MARKERS
)
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH
from db_handler import DbHandler, DbUtils

class AncestryEnvironmentHandler:
    def __init__(
        self,
        vcf_files: dict[str, str] | None = None,  # chromosome -> per-chromosome VCF path
        bed_file: str | None = None,
        output_dir: str = ".",
        reference_vcf_dir: str = REFERENCE_VCF_DIR,  # Your existing 1000G HGDP VCFs
        reference_vcf_pattern: str = REFERENCE_VCF_PATTERN,
        population_panel_file: str = POPULATION_PANEL_FILE,  # Population assignments
        k_populations: int = K_POPULATIONS,  # Number of ancestral populations
        ancestry_results = None
    ):
        if not vcf_files:
            raise ValueError("vcf_files must contain at least one chromosome -> VCF path mapping")
        self.vcf_files = vcf_files
        self.bed_file = bed_file
        self.output_dir = output_dir
        self.reference_vcf_dir = reference_vcf_dir
        self.reference_vcf_pattern = reference_vcf_pattern
        self.population_panel_file = population_panel_file
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

    def _index_vcf(self, vcf_file: str) -> str:
        """Ensure a tabix index exists for vcf_file."""
        subprocess.run(["bcftools", "index", "-t", "-f", vcf_file], check=True)
        return vcf_file
    
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

    def _normalize_sample_order(self, vcf_files: list[str]) -> list[str]:
        """Ensure chromosome VCFs use the same sample order before concatenation."""
        def get_sample_ids(vcf_file: str) -> list[str]:
            result = subprocess.run(
                ["bcftools", "query", "-l", vcf_file],
                check=True,
                capture_output=True,
                text=True
            )
            return result.stdout.splitlines()

        canonical_ids = get_sample_ids(vcf_files[0])
        canonical_set = set(canonical_ids)
        normalized_vcfs = [vcf_files[0]]

        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=self.env.output_dir,
            prefix="ancestry_sample_order_",
            suffix=".txt"
        ) as sample_order_file:
            sample_order_file.write("\n".join(canonical_ids) + "\n")
            sample_order_file.flush()

            for vcf_file in vcf_files[1:]:
                sample_ids = get_sample_ids(vcf_file)
                if sample_ids == canonical_ids:
                    normalized_vcfs.append(vcf_file)
                    continue

                sample_set = set(sample_ids)
                if len(sample_ids) != len(canonical_ids) or sample_set != canonical_set:
                    missing = sorted(canonical_set - sample_set)
                    extra = sorted(sample_set - canonical_set)
                    raise ValueError(
                        f"Chromosome VCF sample set differs for {vcf_file}: "
                        f"missing={missing[:10]}, extra={extra[:10]}"
                    )

                normalized_vcf = os.path.join(
                    self.env.output_dir,
                    f"{Path(vcf_file).name}.sample_ordered.vcf.gz"
                )
                subprocess.run([
                    "bcftools", "view",
                    "-S", sample_order_file.name,
                    "-Oz", "-o", normalized_vcf,
                    vcf_file
                ], check=True)
                subprocess.run(["bcftools", "index", "-t", normalized_vcf], check=True)
                normalized_vcfs.append(normalized_vcf)

        return normalized_vcfs
    
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

    def _validate_user_genotypes(self, plink_prefix: str) -> None:
        """Reject user samples whose observed calls are effectively all heterozygous."""
        qc_prefix = f"{plink_prefix}.genotype_qc"
        subprocess.run([
            "plink2",
            "--bfile", plink_prefix,
            "--het",
            "--out", qc_prefix
        ], check=True)

        het_df = pd.read_csv(f"{qc_prefix}.het", sep=r'\s+')
        user_df = het_df[het_df['IID'].str.startswith('IMPID_')].copy()
        if user_df.empty:
            raise ValueError("No IMPID user samples were found in the ancestry dataset")

        user_df['homozygous_fraction'] = user_df['O(HOM)'] / user_df['OBS_CT']
        invalid = user_df[
            (user_df['OBS_CT'] >= 1000) &
            (user_df['homozygous_fraction'] < 0.05)
        ]
        if not invalid.empty:
            sample_ids = ', '.join(invalid['IID'].head(10))
            raise ValueError(
                "Invalid ancestry input: user genotypes are almost entirely heterozygous "
                f"for {len(invalid)} sample(s), including {sample_ids}. "
                "Check VCF standardization and merging before running ADMIXTURE."
            )
    
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
        
        merged_vcfs = []

        # Step 1-3: Process each chromosome's already-merged VCF separately
        for chrom, sample_chr_vcf in self.env.vcf_files.items():
            print(f"\n=== Processing chromosome {chrom} ===")
            
            # Get reference VCF for this chromosome
            ref_vcf = self._get_reference_vcf_for_chr(chrom)
            print(f"Using reference: {ref_vcf}")
            
            self._index_vcf(sample_chr_vcf)
            
            # Merge sample with reference at common positions
            merged_prefix = os.path.join(self.env.output_dir, f"merged_chr{chrom}")
            merged_vcf = self._merge_sample_with_reference(sample_chr_vcf, ref_vcf, merged_prefix)
            merged_vcfs.append(merged_vcf)
        
        # Step 4: Concatenate all chromosome VCFs
        print(f"\n=== Concatenating {len(merged_vcfs)} chromosome VCFs ===")
        merged_vcfs = self._normalize_sample_order(merged_vcfs)
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
        self._validate_user_genotypes(plink_prefix)

        # Guard against running (supervised) ADMIXTURE on too few markers, which
        # produces unstable/meaningless ancestry proportions that would otherwise
        # be silently uploaded to the database.
        with open(f"{plink_prefix}.bim") as bim_file:
            marker_count = sum(1 for _ in bim_file)
        print(f"Total markers available for ADMIXTURE: {marker_count}")
        if marker_count < MIN_ANCESTRY_MARKERS:
            raise ValueError(
                f"Only {marker_count} markers survived intersection with the reference "
                f"panel (minimum required: {MIN_ANCESTRY_MARKERS}). Ancestry results would "
                "be unreliable - check chromosome coverage and REF/ALT harmonization."
            )
        
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
