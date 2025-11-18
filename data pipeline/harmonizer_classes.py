import pandas as pd
import random
import os
import re
from dataclasses import dataclass
from typing import List, Optional
import subprocess
import sys 
import gzip

@dataclass
class DataContainer:
    microarray_data: Optional[pd.DataFrame] = None
    reference_data: Optional[pd.DataFrame] = None
    harmonized_data: Optional[pd.DataFrame] = None
    harmonization_stats: Optional[pd.DataFrame] = None
    

    
    def _fix_one_allele(self, allele, ref, alts):
        """
        Return (new_allele, flipped) where:
        - new_allele is the original allele or its complement if that matches
        - flipped is True iff we changed to the complement
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        if pd.isna(allele):
            return allele, False

        allele = str(allele).upper()

        if allele == ref or allele in alts:
            return allele, False

        comp = complement.get(allele)
        if comp == ref or comp in alts:
            return comp, True

        return allele, False


    def flip_alleles_to_match_ref_alt(self) -> None:
        """
        For each row, ensure allele1/allele2 match REF or ALT (multi-allelic supported).
        If not, try complement; if complement matches, replace allele and mark *_flipped=True.
        """
        out = self.harmonized_data

        # Normalize strings (only for columns that exist)
        for col in ("REF", "ALT", "allele1", "allele2"):
            if col in out.columns:
                out[col] = out[col].astype("string").str.upper()
                out[col] = out[col].replace(["NAN", "NONE", "NULL", ""], pd.NA)

        def upper_or_empty(value):
            return "" if pd.isna(value) else str(value).upper()

        def allele_ok(allele, ref, alts):
            return (not pd.isna(allele)) and (allele == ref or allele in alts)

        def process_row(r):
            ref = upper_or_empty(r.get("REF"))
            alt_str = upper_or_empty(r.get("ALT"))
            alts = {x for x in alt_str.split(",") if x and x != "."}

            a1_new, a1_flip = self._fix_one_allele(r.get("allele1"), ref, alts)
            a2_new, a2_flip = self._fix_one_allele(r.get("allele2"), ref, alts)

            a1_ok = allele_ok(a1_new, ref, alts)
            a2_ok = allele_ok(a2_new, ref, alts)

            return pd.Series({
                "allele1": a1_new,
                "allele2": a2_new,
                "allele1_flipped": a1_flip,
                "allele2_flipped": a2_flip,
                "allele1_ok": a1_ok,
                "allele2_ok": a2_ok,
                "both_ok": a1_ok and a2_ok,
            })

        cols = [
            "allele1", "allele2",
            "allele1_flipped", "allele2_flipped",
            "allele1_ok", "allele2_ok", "both_ok"
        ]
        out[cols] = out.apply(process_row, axis=1)

        self.harmonized_data = out
    
    def set_harmonization_stats(self) -> None:
        total = len(self.harmonized_data)

        # Booleans sum to counts; cast to int just in case columns are 0/1 ints
        flipped_allele1 = int(self.harmonized_data["allele1_flipped"].sum())
        flipped_allele2 = int(self.harmonized_data["allele2_flipped"].sum())
        ok_allele1      = int(self.harmonized_data["allele1_ok"].sum())
        ok_allele2      = int(self.harmonized_data["allele2_ok"].sum())
        ok_both         = int(self.harmonized_data["both_ok"].sum())

        pct = (lambda x: (x / total * 100.0) if total else 0.0)

        data = pd.DataFrame()


        data = pd.DataFrame([{
            "total_snps": total,
            "flipped_allele1": flipped_allele1,
            "flipped_allele1_pct": round(pct(flipped_allele1), 2),
            "flipped_allele2": flipped_allele2,
            "flipped_allele2_pct": round(pct(flipped_allele2), 2),
            "ok_allele1": ok_allele1,
            "ok_allele1_pct": round(pct(ok_allele1), 2),
            "ok_allele2": ok_allele2,
            "ok_allele2_pct": round(pct(ok_allele2), 2),
            "ok_both": ok_both,
            "ok_both_pct": round(pct(ok_both), 2),
        }])

        self.harmonization_stats = data
    
    def harmonize_data(self) -> None:
        

        self.microarray_data[['allele1', 'allele2']] = self.microarray_data['genotype'].str.extract(r'([ACGT])([ACGT])')

    
        #left join reference_data onto user_data on rsid
        self.harmonized_data = pd.merge(self.microarray_data, self.reference_data, left_on="# rsid", right_on="ID", how="left")


        self.flip_alleles_to_match_ref_alt()
        self.set_harmonization_stats()


        self.harmonized_data['genotype'] = self.harmonized_data['allele1'] + self.harmonized_data['allele2']
        self.harmonized_data = self.harmonized_data[self.harmonized_data["both_ok"] == True].drop(columns=["ID", "allele1_flipped", "allele2_flipped", "allele1_ok", "allele2_ok", "both_ok", "REF", "ALT", "allele1", "allele2"])


        
class EnvironmentHandler:
    def __init__(self,
                 working_dir: str,
                 user_upload_file: str,
                 plink_1_9_path: str,
                 plink_2_0_path: str,
                 plink_map_file: str,
                 pvar_ref_file: str,
                 PLINK_PREFIX: str,
                 plink_reference_fasta: str,
                 beagle_references: str,
                 chromosome_split_files: dict = None,
                 vcf_plink_reference_mapping: pd.DataFrame = None,
                 user_snp_list_path: str = None,
                 reference_data_path: str = None,
                 split_harmonized_file_paths: dict[str, str] = None,
                 bed_file_paths: dict[str, str] = None,
                 vcf_file_paths: dict[str, str] = None
                 ):
        
        self.working_dir = working_dir
        self.user_upload_file = user_upload_file
        self.plink_1_9_path = plink_1_9_path
        self.plink_2_0_path = plink_2_0_path
        self.plink_map_file = plink_map_file
        self.PLINK_PREFIX = PLINK_PREFIX
        self.plink_reference_fasta = plink_reference_fasta
        self.chromosome_split_files = chromosome_split_files
        self.vcf_plink_reference_mapping = vcf_plink_reference_mapping
        self.user_snp_list_path = user_snp_list_path
        self.pvar_ref_file = pvar_ref_file
        self.reference_data_path = reference_data_path
        self.split_harmonized_file_paths = split_harmonized_file_paths
        self.bed_file_paths = bed_file_paths
        self.vcf_file_paths = vcf_file_paths
        self.beagle_references = beagle_references
        self.validate_paths()
    
    
    def validate_paths(self) -> None:
        
        for path in [
            self.working_dir,
            self.user_upload_file,
            self.plink_1_9_path,
            self.plink_2_0_path,
            self.plink_map_file,
            self.pvar_ref_file,
            self.plink_reference_fasta
        ]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")
    
    def set_beagle_files(self) -> dict:
        beagle_files = {}
        for file in os.listdir(self.beagle_references):
            if file.endswith(".bref3"):
                chrom = re.search(r"chr(\d+)\.", file)
                if chrom:
                    beagle_files[chrom.group(1)] = os.path.join(self.beagle_references, file)
            
        self.beagle_references = beagle_files
    
    
    
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
                    data_container: DataContainer,
                    working_dir: str
                    ):
        
        self.environment_handler = environment_handler
        self.data_container = data_container
        self.working_dir = working_dir


    def initiate_working_directory(self) -> None:
        
        random_dir = f"temp_working_dir_{random.randint(1000,9999)}"
        if not os.path.exists(f"{self.environment_handler.working_dir}/{random_dir}"):
            os.makedirs(f"{self.environment_handler.working_dir}/{random_dir}")
            self.working_dir = f"{self.environment_handler.working_dir}/{random_dir}"
        else:
            self.initiate_working_directory()
    
    def set_microarray_data(self) -> pd.DataFrame:
        self.data_container.microarray_data = self.environment_handler.normalize_file()
        
    def create_user_snp_list(self) -> None:
        user_snps_path = os.path.join(self.working_dir, "user_snps.txt")
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
            "--out", f"{self.working_dir}/subset_hg37"
        ]
        
        self.run_command(command)
        self.environment_handler.reference_data_path = f"{self.working_dir}/subset_hg37.pvar"
    
    
    
    def run_command(self, command: list[str]) -> None:
        """Run a command in the subprocess and handle errors."""
        try:
            print(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Error running command: {' '.join(command)}\n{e.stderr}"
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
        harmonized_dir = os.path.join(self.working_dir, "harmonized_chromosomes")
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


        output_dir = os.path.join(self.working_dir, "bed_files")
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
    
    def create_vcf_reference_mapping(self) -> pd.DataFrame:
        
        mapping_df = pd.DataFrame()
        
        
        print(list(self.environment_handler.vcf_file_paths.items()))
        print("-----------------")
        print(list(self.environment_handler.beagle_references.items()))
        vcf_files_df=pd.DataFrame(
            list(self.environment_handler.vcf_file_paths.items()),
            columns=["chromosome_number", "vcf_file"]
        )
        
        beagle_reference_df=pd.DataFrame(
            list(self.environment_handler.beagle_references.items()),
            columns=["chromosome_number", "reference_file"]
        )
        
        plink_map_file = [self.environment_handler.plink_map_file] * len(vcf_files_df)
        
        #join vcf and mapping_df
        mapping_df = pd.merge(
            vcf_files_df,
            beagle_reference_df,
            on="chromosome_number",
            how="inner"
        )
        
        mapping_df["plink_map_file"] = plink_map_file
        
        return mapping_df

    
    

if __name__ == "__main__":
    from config import TEMP_DIR, PLINK_PREFIX, PLINK_MAP_PATH, BEAGLE_REFERENCE_DIR, test_file, GWAS_CATALOG_SNPS, PLINK_1_9_PATH, PLINK_2_0_PATH, PVAR_REF_FILE, PLINK_REFERENCE_FASTA
    
    environment_handler = EnvironmentHandler(
        working_dir=TEMP_DIR,
        user_upload_file=test_file,
        plink_1_9_path=PLINK_1_9_PATH,
        plink_2_0_path=PLINK_2_0_PATH,
        plink_map_file=PLINK_MAP_PATH,
        PLINK_PREFIX=PLINK_PREFIX,
        plink_reference_fasta=PLINK_REFERENCE_FASTA,
        pvar_ref_file=PVAR_REF_FILE,
        beagle_references=BEAGLE_REFERENCE_DIR
    )
    
    environment_handler.set_beagle_files()

    data_container = DataContainer()
    
    workflow_orchestrator = WorkflowOrchestrator(
        
        environment_handler=environment_handler,
        data_container=data_container,
        working_dir=TEMP_DIR
        
    )
    

    workflow_orchestrator.initiate_working_directory()
    print(f"Working directory initialized at: {workflow_orchestrator.working_dir}")
    
    workflow_orchestrator.set_microarray_data()
    print("Microarray data set in data container:")
    print(workflow_orchestrator.data_container.microarray_data.head())
    
    workflow_orchestrator.create_user_snp_list()
    print("User SNP list created.")
    
    print("Extracting reference data...")
    workflow_orchestrator.extract_reference_data()
    workflow_orchestrator.data_container.reference_data = workflow_orchestrator.read_vcf_like_to_df()
    print("Reference data extracted and set in data container")
    
    print("Harmonizing data...")
    workflow_orchestrator.data_container.harmonize_data()
    print("Data harmonization complete. Harmonized data:")
    print(workflow_orchestrator.data_container.harmonization_stats)
    print(workflow_orchestrator.data_container.harmonized_data.head())
    
    print("Splitting harmonized data into chromosome-specific files...")
    workflow_orchestrator.create_harmonized_chromosome_files()
    print("Splitting complete. Chromosome-specific file paths:")
    
    print("Converting harmonized files to BED format...")
    workflow_orchestrator.convert_23andme_to_bed()
    workflow_orchestrator.confirm_paths_exist(workflow_orchestrator.environment_handler.bed_file_paths)
    print("Conversion to BED format complete. BED file paths:")
    
    print("Converting BED files to VCF format...")
    workflow_orchestrator.convert_bed_to_vcf()
    workflow_orchestrator.confirm_paths_exist(workflow_orchestrator.environment_handler.vcf_file_paths)
    print("Conversion to VCF format complete. VCF file paths:")
    print(workflow_orchestrator.environment_handler.vcf_file_paths)
    print(workflow_orchestrator.environment_handler.split_harmonized_file_paths)
    

    print("Creating VCF to reference mapping...")
    vcf_reference_mapping_df = workflow_orchestrator.create_vcf_reference_mapping()
    print("VCF to reference mapping created:")
    print(vcf_reference_mapping_df)
    
    
    
    
    
    
    
    
    
    
    