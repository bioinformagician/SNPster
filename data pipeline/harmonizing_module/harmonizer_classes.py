import os
import re
import subprocess
import gzip
from data_models import *
        
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
    
    
    
    
    
    def normalize_file(self) -> pd.DataFrame:
            
            #'#rsidchromosomepositiongenotype'
            normalized_header = ['# rsid', 'chromosome', 'position', 'genotype']
            allowed_headers = [
                '#rsidchromosomepositiongenotype', #23andme
                'rsidchromosomepositiongenotype', #23andme without #
                'rsidchromosomepositionallele1allele2', #ancestry
                'RSID,CHROMOSOME,POSITION,RESULT', #ftdna
                'Name,Variation,Chromosome,Position,Strand,YourCode' #myheritage
            ]
            
            header = self.extract_header()
            
            if header not in allowed_headers:
                raise ValueError(f"Header format not recognized. Found: {header}")
            
            data = self.extract_data()
            
            
            if header == '#rsidchromosomepositiongenotype': #is 23andme
                data.columns = normalized_header
            
            if header == 'rsidchromosomepositiongenotype': #is 23andme without #
                #rename header to standardized header
                data.rename(columns={'rsid': '# rsid'}, inplace=True)

            elif header == 'rsidchromosomepositionallele1allele2': #is ancestry
                data.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2']
                #combine allele1 and allele2 into genotype
                data['genotype'] = data['allele1'] + data['allele2']
                data = data.drop(columns=['allele1', 'allele2'])
            
            elif header == 'RSID,CHROMOSOME,POSITION,RESULT': #is ftdna
                data= data.rename(columns={'RESULT': 'genotype'})
                data.columns = [col.lower() for col in data.columns]
            
            elif header == 'Name,Variation,Chromosome,Position,Strand,YourCode': #is myheritage
                data.drop(columns=['Strand', 'Variation'], inplace=True)
                data.rename(columns={'Position': 'position', 'Chromosome': 'chromosome', 'YourCode': 'genotype', 'Name': '# rsid'}, inplace=True)
                
            
            
            #keep only clean genotypes
            data = data[data['genotype'].str.match("^[ATCG]{1,2}$")]

            # clean rsids
            
            data = data.rename(columns={'rsid': '# rsid'})

            data = data[data['# rsid'].str.match("^rs[0-9]{1,}$")]


            #last step to ensure correct column order
            data = data[['# rsid','chromosome','position','genotype']]


            return data
        
    
    
    def extract_header(self) -> str:
        """
        Find the first line that contains an rsxxxx pattern (rs followed by digits)
        and return the line immediately above it with all whitespace removed.
        """
        pattern = re.compile(r'\brs\d+\b', flags=re.IGNORECASE)

        with open(self.user_upload_file, 'r') as file:
            lines = file.readlines()

        # Find the first line that matches the rs pattern
        rs_line_index = None
        for i, line in enumerate(lines):
            if pattern.search(line):
                rs_line_index = i
                break

        if rs_line_index is None:
            raise ValueError("No line containing an rsxxxx pattern was found.")

        if rs_line_index == 0:
            raise ValueError("The rsxxxx pattern appears on the first line; no line above to extract.")

        header_line = lines[rs_line_index - 1]

        # Flatten: remove all whitespace characters (spaces, tabs, newlines)
        header_line = re.sub(r'\s+', '', header_line)

        return header_line
    
    def extract_data(self) -> pd.DataFrame:
        try: 
            #check if a # is present in the first line
            with open(self.user_upload_file, 'r') as file:
                first_line = file.readline()
                has_hash = first_line.startswith('#')
            
            
            if has_hash:
                df = pd.read_csv(self.user_upload_file, sep=r'\s+', comment="#", dtype=str, low_memory=False)
            else:
                #no comment, has header
                df = pd.read_csv(self.user_upload_file, sep=r'\s+', dtype=str, low_memory=False)
            
            #check number of columns
            if df.shape[1] < 4:

                    df = pd.read_csv(
                        self.user_upload_file,
                        sep=',',
                        engine='c',
                        dtype=str,
                        low_memory=True,
                        header=0,       # the CSV header line (e.g., RSID,CHROMOSOME,POSITION,RESULT)
                    )
                    # strip quotes if present
                    for c in df.columns:
                        df[c] = df[c].astype(str).str.replace('"', '', regex=False)

            if df.shape[1] < 4:
                raise ValueError("File does not have enough columns to extract rsid, chromosome, position, and genotype.")
            
        except pd.errors.ParserError:
                df = pd.read_csv(
                    self.user_upload_file,
                    comment="#",          # skip the header blurb
                    sep=",",              # CSV
                    quotechar='"',        # quoted fields like "rs4477212"
                    dtype=str,            # keep everything as strings
                    keep_default_na=False,# don't convert strings like "NA" to NaN
                    engine="c",           # fast
                )

        return df
        



class WorkflowOrchestrator:
    def __init__(self,
                    environment_handler: EnvironmentHandler,
                    data_container: DataContainer
                    ):
        
        self.environment_handler = environment_handler
        self.data_container = data_container

    
    def set_microarray_data(self) -> pd.DataFrame:
        self.data_container.microarray_data = self.environment_handler.normalize_file()
        
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
        #beware of the 
        
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
                    "--split-par", "b37", #hg38 for grch38 genome build 
                    "--fa", self.environment_handler.plink_reference_fasta,
                    "--ref-from-fa",
                    "--export", "vcf", "bgz",
                    "--out", rf"{filename}"
                ]
                
                self.run_command(command)
                output_vcf_files[chr_number] = rf"{filename}.vcf.gz"

                
        self.environment_handler.vcf_file_paths = output_vcf_files
    
    def confirm_paths_exist(self, paths: dict[str, str]) -> None:
        for path in paths.values():
            if not os.path.exists(path):
                raise FileNotFoundError(f"Harmonized chromosome file not found: {path}")
    


    
    
    
    
    
    