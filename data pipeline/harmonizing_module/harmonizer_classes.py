import os
import subprocess
import gzip
from data_models import *
import pyarrow.parquet as pq

class EnvironmentHandler:
    def __init__(self,
                 output_dir: str,
                 user_upload_file: str,
                 plink_1_9_path: str,
                 plink_2_0_path: str,
                 pvar_ref_file: str,
                 PLINK_PREFIX: str,
                 plink_reference_fasta: str,
                 plink_1_9_memory_mb: str,
                 plink_1_9_threads: str,
                 chromosome_split_files: dict = None,
                 user_snp_list_path: str = None,
                 reference_data_path: str = None,
                 split_harmonized_file_paths: dict[str, str] = None,
                 bed_file_paths: dict[str, str] = None,
                 vcf_file_paths: dict[str, str] = None
                 ):
        
        self.output_dir = output_dir
        self.user_upload_file = user_upload_file
        self.plink_1_9_path = plink_1_9_path
        self.plink_2_0_path = plink_2_0_path
        self.PLINK_PREFIX = PLINK_PREFIX
        self.plink_reference_fasta = plink_reference_fasta
        self.chromosome_split_files = chromosome_split_files
        self.user_snp_list_path = user_snp_list_path
        self.pvar_ref_file = pvar_ref_file
        self.reference_data_path = reference_data_path
        self.split_harmonized_file_paths = split_harmonized_file_paths
        self.bed_file_paths = bed_file_paths
        self.vcf_file_paths = vcf_file_paths
        self.plink_1_9_memory_mb = plink_1_9_memory_mb
        self.plink_1_9_threads = plink_1_9_threads
        self.create_output_directory()
        self.validate_paths()
        
        
    def create_output_directory(self) -> None:
        os.makedirs(self.output_dir, exist_ok=True)
    
    
    def validate_paths(self) -> None:
        
        for path in [
            self.output_dir,
            self.user_upload_file,
            self.plink_1_9_path,
            self.plink_2_0_path,
            self.pvar_ref_file,
            self.plink_reference_fasta
        ]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")
    
    
        



class WorkflowOrchestrator:
    def __init__(self,
                    environment_handler: EnvironmentHandler,
                    data_container: DataContainer,
                    ):
        
        self.environment_handler = environment_handler
        self.data_container = data_container

    def read_parquet(self) -> None:
        
        
        # Read parquet file with metadata
        parquet_file = pq.read_table(self.environment_handler.user_upload_file)
        self.data_container.microarray_data = parquet_file.to_pandas()
        
        # Extract custom metadata from schema
        metadata = parquet_file.schema.metadata

        self.data_container.vendor = metadata[b'vendor'].decode('utf-8')
        self.data_container.genome_build = metadata[b'genome_build'].decode('utf-8')
        self.data_container.identifier = metadata[b'identifier'].decode('utf-8')
        
        self.data_container.is_forward_strand = metadata[b'is_forward_strand'].decode('utf-8')
        if self.data_container.is_forward_strand == 'True':
            self.data_container.is_forward_strand = True
        else:
            self.data_container.is_forward_strand = False

        
        print(f"Loaded parquet with metadata: vendor={self.data_container.vendor}, genome_build={self.data_container.genome_build}, is_forward_strand={self.data_container.is_forward_strand}")
        
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
            self.environment_handler.plink_2_0_path,
            "--pvar", self.environment_handler.pvar_ref_file,
            "--extract", self.environment_handler.user_snp_list_path,
            "--make-just-pvar",
            "--threads", "1",
            "--memory", "8000", #self.environment_handler.plink_1_9_memory_mb,
            "--out", f"{self.environment_handler.output_dir}/subset_hg37"
        ]
        
        self.run_command(command)
        self.environment_handler.reference_data_path = f"{self.environment_handler.output_dir}/subset_hg37.pvar"
    
    
    
    def run_command(self, command: list[str]) -> None:
        """Run a command in the subprocess and handle errors."""
        try:
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
    
    def read_vcf_like_to_df(self) -> pd.DataFrame:
        """
        Read a VCF/pvar into a DataFrame using pandas.
        - Skips header lines starting with '#'
        - Uses the '#CHROM ...' line for column names
        - Optionally expands selected INFO keys into separate columns
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
    
    
    
    def create_harmonized_chromosome_files(self) -> None:
        
        #make new dir for harmonized chromosome files
        harmonized_dir = os.path.join(self.environment_handler.output_dir, "harmonized_chromosomes")
        os.makedirs(harmonized_dir, exist_ok=True)
        
        output_filepaths = {}

        dataframe_split = {chrom: df for chrom, df in self.data_container.harmonized_data.groupby('chromosome')}

        for chrom, df in dataframe_split.items():
            
        #write each dataframe to a separate file in the temp_dir with the same header as the original file
            
            output_file = os.path.join(harmonized_dir, f"chr{chrom}.txt")
            df.to_csv(output_file, sep="\t", index=False)
            print(f"Wrote chromosome {chrom} to {output_file}")
            output_filepaths[chrom] = output_file

        self.environment_handler.split_harmonized_file_paths = output_filepaths
    
    def convert_23andme_to_bed(self) -> None:
        """Convert 23andMe text file to PLINK binary format."""


        output_dir = os.path.join(self.environment_handler.output_dir, "bed_files")
        os.makedirs(output_dir, exist_ok=True)
        output_paths = {}
        
        for chr_number, file in self.environment_handler.split_harmonized_file_paths.items():
            print(f"Converting file: {file} to BED format")

            filename = os.path.basename(file)
            output_path = os.path.join(output_dir, f"{self.environment_handler.PLINK_PREFIX}_{filename.replace('.txt', '')}")
            
            command = [
                self.environment_handler.plink_1_9_path,
                "--23file", file,
                "FAM001", "ID001",
                "--memory", self.environment_handler.plink_1_9_memory_mb,
                "--threads", self.environment_handler.plink_1_9_threads,
                "--make-bed",
                "--allow-no-sex",
                "--out", output_path
            ]
            
            self.run_command(command)
            output_paths[chr_number] = output_path+".bed"
            
        self.environment_handler.bed_file_paths = output_paths
    
    
    
    def convert_bed_to_vcf(self) -> None:
        output_vcf_files = {}
        for chr_number, file in self.environment_handler.bed_file_paths.items():
                print(f"Processing BED file: {file} to VCF format")
                
                filename = file.replace(".bed", "") #the extension is removed because the plink program will need the .bed, .bim, and .fam files, therefore we provide the generic filename. It will find all three files automatically from the generic filename 
                """Convert PLINK binary files to VCF format using a reference genome."""
                command = [
                    self.environment_handler.plink_2_0_path,
                    "--bfile", filename,
                    "--split-par", "b38", #hg38 for grch38 genome build 
                    "--fa", self.environment_handler.plink_reference_fasta,
                    "--ref-from-fa",
                    "--export", "vcf", "bgz",
                    "--out", rf"{filename}"
                ]
                
                self.run_command(command)
                output_vcf_files[chr_number] = rf"{filename}_{self.data_container.identifier}.vcf.gz"

                
        self.environment_handler.vcf_file_paths = output_vcf_files
    
    def confirm_paths_exist(self, paths: dict[str, str]) -> None:
        for path in paths.values():
            if not os.path.exists(path):
                raise FileNotFoundError(f"Harmonized chromosome file not found: {path}")
    
    
    def run_harmonization_workflow(self) -> None:
        
        """In the future the genome build should also be evaluated to ensure proper harmonization"""
        
        if self.data_container.is_forward_strand is False:
            print("Data is not in forward strand orientation. Proceeding with harmonization workflow...")
            self.create_user_snp_list()
            print("User SNP list created.")
            print("Extracting reference data...")
            self.extract_reference_data() #extracts refernce data based on snp list
            self.data_container.reference_data = self.read_vcf_like_to_df() #reads extracted refernce data
            print("Reference data extracted and set in data container")
            print("Harmonizing data...")
            self.data_container.harmonize_data()
            print("Data harmonization complete. Harmonized data:")
            print(self.data_container.harmonization_stats)
            print(self.data_container.harmonized_data.head())
        else:
            print("Data is already in forward strand orientation. Skipping harmonization workflow...")
            self.data_container.harmonized_data = self.data_container.microarray_data

    


    
    
    
    
    
    