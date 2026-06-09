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
                 PLINK_PREFIX: str,
                 plink_reference_fasta: str,
                 plink_1_9_memory_mb: str,
                 plink_1_9_threads: str,
                 bed_file_path: str = None,
                 vcf_file_path: str = None
                 ):
        
        self.output_dir = output_dir
        self.user_upload_file = user_upload_file
        self.plink_1_9_path = plink_1_9_path
        self.plink_2_0_path = plink_2_0_path
        self.PLINK_PREFIX = PLINK_PREFIX
        self.plink_reference_fasta = plink_reference_fasta
        self.bed_file_path = bed_file_path
        self.vcf_file_path = vcf_file_path
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
        self.data_container.harmonized_data = parquet_file.to_pandas()
        
        # Extract custom metadata from schema
        metadata = parquet_file.schema.metadata

        self.data_container.vendor = metadata[b'vendor'].decode('utf-8')
        self.data_container.genome_build = metadata[b'genome_build'].decode('utf-8')
        self.data_container.imputation_id = metadata[b'imputation_id'].decode('utf-8')
        
        self.data_container.is_forward_strand = metadata[b'is_forward_strand'].decode('utf-8')
        if self.data_container.is_forward_strand == 'True':
            self.data_container.is_forward_strand = True
        else:
            self.data_container.is_forward_strand = False

        
        print(f"Loaded parquet with metadata: vendor={self.data_container.vendor}, genome_build={self.data_container.genome_build}, is_forward_strand={self.data_container.is_forward_strand}")
        

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
    
    
    
    
    def convert_23andme_to_bed(self) -> None:
        """
        Convert 23andMe text file to PLINK binary format.
        """

        output_dir = self.environment_handler.output_dir
        os.makedirs(output_dir, exist_ok=True)

        harmonized_df = self.data_container.harmonized_data.copy()
        
        harmonized_df["chromosome"] = (
            harmonized_df["chromosome"]
            .astype(str)
            .str.replace(r"\.0$", "", regex=True)
        )
        harmonized_df["position"] = (
            pd.to_numeric(harmonized_df["position"], errors="raise")
            .round()
            .astype(int)
        )
        
        chr_number = harmonized_df['chromosome'].iloc[0]
        
        harmonized_df = harmonized_df.sort_values(
            by=["position", "# rsid"],
            kind="mergesort"
        ).reset_index(drop=True)
        
        csv_path = os.path.join(
            output_dir,
            f"IMPID{self.data_container.imputation_id}.chr{chr_number}.HarmonizedMicroarray.csv"
        )

        harmonized_df.to_csv(
            csv_path,
            sep="\t",
            index=False,
            header=False
        )
        
        print(f"Converting file: {csv_path} to BED format")

        filename = os.path.basename(csv_path)
        output_prefix = filename.replace(".HarmonizedMicroarray.csv", "")
        output_path = os.path.join(output_dir, f"{output_prefix}.{self.environment_handler.PLINK_PREFIX}")
        
        command = [
            self.environment_handler.plink_1_9_path,
            "--23file", csv_path,
            "IMPDID", self.data_container.imputation_id,
            "--memory", self.environment_handler.plink_1_9_memory_mb,
            "--threads", self.environment_handler.plink_1_9_threads,
            "--make-bed",
            "--allow-no-sex",
            "--out", output_path
        ]
        
        self.run_command(command)
            
        self.environment_handler.bed_file_path = output_path+".bed"
    
    
    
    def convert_bed_to_vcf(self) -> None:
        # Prepare reference FASTA (decompress and index if needed) once before processing
        prepared_fasta = self.prepare_reference_fasta()
        
        
        file = self.environment_handler.bed_file_path
        
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

                
        self.environment_handler.vcf_file_path = vcf_path
    
    def confirm_path_exist(self, path) -> None:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Harmonized chromosome file not found: {path}")

    def add_imputation_id_to_vcfs(self) -> None:
        
        vcf_utilities.add_imputation_id_to_vcf(vcf_file=self.environment_handler.vcf_file_path, imputation_id=self.data_container.imputation_id)
            
    


    


    
    
    
    
    
    