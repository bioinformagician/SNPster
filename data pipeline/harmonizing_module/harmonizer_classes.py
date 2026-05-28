import os
import subprocess
import gzip
from data_models import *
import pyarrow.parquet as pq
from vcf_classes import VCFUtilities 


vcf_utilities = VCFUtilities()

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
        
        for binary_path in [self.plink_1_9_path, self.plink_2_0_path]:
            if os.path.isfile(binary_path) and not os.access(binary_path, os.X_OK):
                try:
                    os.chmod(binary_path, os.stat(binary_path).st_mode | 0o111)
                except PermissionError:
                    pass
    
    
        



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
            "--out", f"{self.environment_handler.output_dir}/subset_hg38"
        ]
        
        self.run_command(command)
        self.environment_handler.reference_data_path = f"{self.environment_handler.output_dir}/subset_hg38.pvar"
    
    def prepare_reference_fasta(self) -> str:
        """Prepare reference FASTA for bcftools norm (decompress and index if needed).
        
        Returns: Path to the prepared (decompressed and indexed) FASTA file.
        """
        fasta_gz = self.environment_handler.plink_reference_fasta
        
        # Check if it's gzipped
        if not fasta_gz.endswith('.gz'):
            # Already decompressed, just need to index
            fasta_path = fasta_gz
        else:
            # Need to decompress
            fasta_path = fasta_gz.replace('.gz', '')
            
        # Check if decompressed file exists
        if not os.path.exists(fasta_path):
            print(f"Decompressing reference FASTA: {fasta_gz} -> {fasta_path}")
            with gzip.open(fasta_gz, 'rb') as f_in:
                with open(fasta_path, 'wb') as f_out:
                    # Decompress in chunks to handle large files
                    chunk_size = 1024 * 1024  # 1MB chunks
                    while True:
                        chunk = f_in.read(chunk_size)
                        if not chunk:
                            break
                        f_out.write(chunk)
            print(f"Decompression complete: {fasta_path}")
        
        # Check if index exists
        fai_path = fasta_path + '.fai'
        if not os.path.exists(fai_path):
            print(f"Creating FASTA index: {fai_path}")
            index_command = ["samtools", "faidx", fasta_path]
            self.run_command(index_command)
            print(f"Index created: {fai_path}")
        
        return fasta_path
    
    
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

        # Define valid chromosomes (standard autosomes, sex chromosomes, and mitochondrial)
        valid_chromosomes = {str(i) for i in range(1, 23)} | {'X', 'Y', 'MT', 'M'}
        # Also accept chr-prefixed versions
        valid_chromosomes |= {f'chr{c}' for c in list(valid_chromosomes)}
        
        dataframe_split = {chrom: df for chrom, df in self.data_container.harmonized_data.groupby('chromosome')}

        for chrom, df in dataframe_split.items():
            # Skip alternate contigs, patches, and unplaced sequences
            chrom_str = str(chrom)
            # Extract base chromosome (remove chr prefix if present for comparison)
            base_chrom = chrom_str.replace('chr', '') if chrom_str.startswith('chr') else chrom_str
            
            # Check if this is a standard chromosome (not alt, patch, or unplaced)
            if base_chrom not in valid_chromosomes and not any(base_chrom.startswith(str(i)) and '_' not in base_chrom for i in range(1, 23)):
                print(f"Skipping non-standard chromosome: {chrom_str}")
                continue
            
            # Ensure proper chromosome naming (add chr prefix if not present)
            if not chrom_str.startswith('chr'):
                chrom_str = f'chr{chrom_str}'
            
            output_file = os.path.join(harmonized_dir, f"{chrom_str}.txt")
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
            output_path = os.path.join(output_dir, f"{self.environment_handler.PLINK_PREFIX}_{filename.replace('.txt', '')}_{self.data_container.identifier}")
            
            command = [
                self.environment_handler.plink_1_9_path,
                "--23file", file,
                "FAM001", self.data_container.identifier,
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
        # Prepare reference FASTA (decompress and index if needed) once before processing
        prepared_fasta = self.prepare_reference_fasta()
        
        output_vcf_files = {}
        for chr_number, file in self.environment_handler.bed_file_paths.items():
                print(f"Processing BED file: {file} to VCF format")
                
                filename = file.replace(".bed", "") #the extension is removed because the plink program will need the .bed, .bim, and .fam files, therefore we provide the generic filename. It will find all three files automatically from the generic filename
                 
                """Convert PLINK binary files to VCF format using a reference genome.
                Using --ref-from-fa to ensure all samples have consistent REF alleles matching the reference FASTA.
                Then normalize with bcftools to ensure consistent REF orientation across all samples.
                """
                command = [
                    self.environment_handler.plink_2_0_path,
                    "--bfile", filename,
                    "--split-par", "b38", #hg38 for grch38 genome build 
                    "--fa", self.environment_handler.plink_reference_fasta,
                    "--ref-from-fa", "force",
                    "--export", "vcf", "bgz",
                    "--out", rf"{filename}"
                ]
                
                self.run_command(command)
                
                # Normalize VCF to ensure consistent REF alleles across all samples
                # This fixes cases where PLINK sets different REF for homozygous vs heterozygous sites
                vcf_path = rf"{filename}.vcf.gz"
                normalized_vcf = rf"{filename}.normalized.vcf.gz"
                
                print(f"Normalizing VCF file: {vcf_path}")
                norm_command = [
                    "bcftools", "norm",
                    "--check-ref", "s",  # Swap REF/ALT if REF doesn't match FASTA
                    "--fasta-ref", prepared_fasta,  # Use decompressed and indexed FASTA
                    "-Oz", "-o", normalized_vcf,
                    vcf_path
                ]
                
                self.run_command(norm_command)
                
                # Replace original VCF with normalized version
                os.replace(normalized_vcf, vcf_path)
                
                # Reindex the normalized VCF
                index_command = ["bcftools", "index", "-f", vcf_path]
                self.run_command(index_command)
                
                output_vcf_files[chr_number] = vcf_path

                
        self.environment_handler.vcf_file_paths = output_vcf_files
    
    def confirm_paths_exist(self, paths: dict[str, str]) -> None:
        for path in paths.values():
            if not os.path.exists(path):
                raise FileNotFoundError(f"Harmonized chromosome file not found: {path}")
    
    def add_imputation_id_to_vcfs(self) -> None:
        
        for chrom, vcf_path in self.environment_handler.vcf_file_paths.items():
            vcf_utilities.add_imputation_id_to_vcf(vcf_file = vcf_path, imputation_id=self.data_container.identifier)
    
    
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

    


    
    
    
    
    
    