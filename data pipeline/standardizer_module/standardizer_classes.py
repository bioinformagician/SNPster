from data_models import DataContainer, FileHandler, pd, os
from pyliftover import LiftOver
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess
import gzip

class EnvironmentHandler:
    
    def __init__(self,
                 chain_file_dict: dict[str, str],
                 grch_to_hg_identifier_dict: dict[str, str],
                 user_file:str,
                 output_dir: str,
                 pvar_ref_file: str
                 ) -> None:
        
        self.output_dir = output_dir
        self.chain_file_dict = chain_file_dict
        self.user_file = user_file
        self.grch_to_hg_identifier_dict = grch_to_hg_identifier_dict
        self.pvar_ref_file = pvar_ref_file
    
    def validate_environment(self) -> None:
        
        os.makedirs(self.output_dir, exist_ok=True)
            
        paths = list(self.chain_file_dict.values()) + [self.user_file, self.pvar_ref_file]
        #remove None values (for GRCh38 chain file)
        paths = [path for path in paths if path is not None]
        for path in paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")
    

class WorkflowOrchestrator:
    def __init__(self,
                    data_container: DataContainer,
                    file_handler: FileHandler,
                    environment_handler: EnvironmentHandler
                    ):
        
        self.data_container = data_container
        self.file_handler = file_handler
        self.environment_handler = environment_handler

    def run_command(self, command: list[str]) -> None:
        """Run a command in the subprocess and handle errors."""
        try:
            if command and os.path.isfile(command[0]) and not os.access(command[0], os.X_OK):
                try:
                    os.chmod(command[0], os.stat(command[0]).st_mode | 0o111)
                except PermissionError:
                    pass
            print(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Error running command: {' '.join(command)}\n"
                f"Return code: {e.returncode}\n"
                f"STDOUT:\n{e.stdout}\n"
                f"STDERR:\n{e.stderr}\n"
            )

    def check_dictionary_coherence(self) -> None:
        """Check that the keys and values in the build extraction tool matches the liftover files dicts"""
        for key in self.file_handler.genome_build_dict.values():
            if key not in self.environment_handler.chain_file_dict.keys():
                raise ValueError(f"Genome build '{key}' in genome_build_dict does not have a corresponding chain file in environment_handler.chain_file_dict")

        for key in self.environment_handler.grch_to_hg_identifier_dict.keys():
            if key not in self.file_handler.genome_build_dict.values():
                raise ValueError(f"Genome build '{key}' in grch_to_hg_identifier_dict does not have a corresponding entry in file_handler.genome_build_dict")
        
        for vendor in self.file_handler.forward_strand_vendors:
            if vendor not in self.file_handler.accepted_vendors_dict.keys():
                raise ValueError(f"Vendor '{vendor}' in forward_strand_vendors is not in accepted_vendors_dict")
        
    def set_vendor(self) -> None:
        self.data_container.vendor = self.file_handler.identify_vendor()
        print(f"Identified vendor: {self.data_container.vendor}")
    
    def set_genome_build(self) -> None:
        self.data_container.genome_build = self.file_handler.identify_genome_build()
        print(f"Identified genome build: {self.data_container.genome_build}")
        
        if self.data_container.genome_build != 'GRCh38':
            self.data_container.lift_over = True
    
    def set_strand_direction(self) -> None:
        if self.data_container.vendor is None:
            raise ValueError("Vendor must be set before determining strand direction.")
        self.data_container.is_forward_strand = self.data_container.vendor in self.file_handler.forward_strand_vendors
            
    
    def set_microarray_data(self) -> pd.DataFrame:
        self.data_container.microarray_data = self.file_handler.normalize_file()
    
    def evaluate_liftover(self) -> None:
        
        if self.data_container.lift_over:
            print(f"Performing liftover from {self.data_container.genome_build} to GRCh38...")
            lo = LiftOver(self.environment_handler.grch_to_hg_identifier_dict[self.data_container.genome_build], 
                          self.environment_handler.grch_to_hg_identifier_dict['GRCh38'])
            self.data_container.lift_over_data(lo)
            # Update genome_build metadata to reflect the liftover
            self.data_container.genome_build = 'GRCh38'
            print(f"Liftover complete. Genome build updated to GRCh38")
    
    def evaluate_zipping(self) -> None:
        if self.file_handler.is_zipped_file():
            unzipped_file = self.file_handler.unzip_file()
            self.file_handler.user_file = unzipped_file
            self.environment_handler.user_file = unzipped_file
        
    
    def write_parquet_output(self) -> None:
        """ Naming of file should follow {IMPIDx}.FILEIDx.chr{ALL}.{STAGE}.filextension"""
        df = self.data_container.harmonized_data
        if df is None:
            raise ValueError("No harmonized data available to write.")

        if self.data_container.file_id is None:
            raise ValueError("file_id must be set before writing parquet output.")

        output_path = os.path.join(
            self.environment_handler.output_dir,
            f"IMPID{self.data_container.imputation_id}.FILEID{self.data_container.file_id}.chrALL.standardizedMicroarray.parquet",
        )

        # Convert DataFrame to PyArrow Table
        table = pa.Table.from_pandas(df)

        # Add custom metadata
        metadata = {
            b'vendor': self.data_container.vendor.encode('utf-8'),
            b'genome_build': self.data_container.genome_build.encode('utf-8'),
            b'is_forward_strand': str(self.data_container.is_forward_strand).encode('utf-8'),
            b'imputation_id': str(self.data_container.imputation_id).encode('utf-8'),
            b'file_id': str(self.data_container.file_id).encode('utf-8')
        }

        # Merge with existing schema metadata
        existing_metadata = table.schema.metadata or {}
        combined_metadata = {**existing_metadata, **metadata}

        # Create new schema with metadata
        new_schema = table.schema.with_metadata(combined_metadata)
        table = table.cast(new_schema)

        # Write parquet with metadata
        pq.write_table(table, output_path)
        print(f"Standardized data written to {output_path}")
        print(f"Metadata: vendor={self.data_container.vendor}, genome_build={self.data_container.genome_build}, is_forward_strand={self.data_container.is_forward_strand}, imputation_id={self.data_container.imputation_id}, file_id={self.data_container.file_id}")
            
            
            
            
            
    
    #----------------- Harmonization methods -----------------#
    
    
    
    
    
    
    def create_user_snp_list(self) -> None:
        user_snps_path = os.path.join(self.environment_handler.output_dir, "user_snps.txt")
        self.data_container.microarray_data["# rsid"].to_csv(user_snps_path, sep="\t", index=False, header=False)
        self.environment_handler.user_snp_list_path = user_snps_path
    
    
    def extract_reference_data(self) -> None:
        r"""example: C:\Users\frezz>"C:\Users\frezz\Desktop\SNPster\gwas_catalog_data_cleaning\plink2_win_avx2\plink2.exe" 
        --pvar "C:\Users\frezz\Desktop\SNPster\gwas_catalog_data_cleaning\get_forward_alleles\all_phase3.pvar.zst" 
        --extract "C:\Users\frezz\Desktop\SNPster\gwas_catalog_data_cleaning\get_forward_alleles\risk_alleles.txt" 
        --make-just-pvar --out "C:\Users\frezz\Downloads\subset_hg37"""
        
        
        
        command = [
            "plink2",
            "--pvar", self.environment_handler.pvar_ref_file,
            "--extract", self.environment_handler.user_snp_list_path,
            "--make-just-pvar",
            "--threads", "1",
            "--memory", "8000", #self.environment_handler.plink_1_9_memory_mb,
            "--out", f"{self.environment_handler.output_dir}/subset_hg38"
        ]
        
        self.run_command(command)
        self.environment_handler.reference_data_path = f"{self.environment_handler.output_dir}/subset_hg38.pvar"
        
            
        
    def read_vcf_like_to_df(self) -> pd.DataFrame:
        """
        Read a VCF/pvar into a DataFrame using pandas.
        - Skips header lines starting with '#'
        - Uses the '#CHROM ...' line for column names
        - Optionally expands selected INFO keys into separate columns
        
        Move this to vcf_utilities module
        """
        # 1) find the column header line
        opener = gzip.open if str(self.environment_handler.reference_data_path).endswith(".gz") else open
        cols = None
        with opener(self.environment_handler.reference_data_path, "rt", encoding="utf-8", newline="") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    # keep the '#CHROM' name or drop the leading '#', either is fine
                    cols = line.rstrip("\n").lstrip("#").split("\t")
                    break
        if cols is None:
            raise ValueError("No #CHROM header line found in file.")

        # 2) read the body; skip all lines beginning with '#'
        df = pd.read_csv(
            self.environment_handler.reference_data_path,
            sep="\t",
            comment="#",
            header=None,
            names=cols,
            dtype={cols[0]: "string"}  # keep chromosome as string (handles 'X','MT')
        )

        return df
    
    def run_harmonization_workflow(self) -> None:
        
        """
        Harmonizes strand orientation by flipping alleles to match the GRCh38 reference.
        Note: Genome build liftover should be handled by the standardizer module.
        This step harmonizes alleles to match the GRCh38 reference, regardless of strand.
        """
        
        # Always harmonize to ensure consistent GRCh38 positions and alleles
        print("Harmonizing data to GRCh38 reference...")
        self.create_user_snp_list()
        print("User SNP list created.")
        print("Extracting reference data...")
        self.extract_reference_data() #extracts reference data based on snp list
        self.data_container.reference_data = self.read_vcf_like_to_df() #reads extracted reference data
        print("Reference data extracted and set in data container")
        print("Harmonizing data...")
        self.data_container.harmonize_data()
        print("Data harmonization complete. Harmonized data:")
        print(self.data_container.harmonization_stats)
        print(self.data_container.harmonized_data.head())

