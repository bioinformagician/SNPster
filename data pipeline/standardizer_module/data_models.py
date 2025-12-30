from typing import Optional
from dataclasses import dataclass
import pandas as pd
import re
import os 
import zipfile

@dataclass
class DataContainer:
    vendor: Optional[str] = None
    genome_build: Optional[str] = None
    lift_over: bool = False
    microarray_data: Optional[pd.DataFrame] = None

    

    def lift_over_data(self, lo) -> None:
        df = self.microarray_data

        chrom = df["chromosome"].astype(str).to_numpy()
        pos   = df["position"].astype(int).to_numpy()

        new_chrom = chrom.copy()
        new_pos   = pos.copy()

        for i in range(len(df)):
            lifted = lo.convert_coordinate("chr" + chrom[i], int(pos[i]))
            if lifted:
                #remove chr prefix
                new_chrom[i] = lifted[0][0].replace("chr", "")
                # lifted[0] is typically: (new_chr, new_pos, strand, ... )
                new_pos[i]   = int(lifted[0][1])

        df["chromosome"] = new_chrom
        df["position"]   = new_pos
        
        print(df)
        
        


class FileHandler:

    def __init__(self, 
                 user_file: Optional[str],
                 accepted_vendors_dict: dict[str, str],
                 genome_build_dict: dict[str, str]
                 ):
        
        self.user_file = user_file
        self.accepted_vendors_dict = accepted_vendors_dict
        self.genome_build_dict = genome_build_dict
    
    def is_zipped_file(self) -> bool:
        return zipfile.is_zipfile(self.user_file)
    
    def unzip_file(self) -> str:

        with zipfile.ZipFile(self.user_file, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(self.user_file))
            extracted_files = zip_ref.namelist()
        
        # Filter out directories, keep only files
        extracted_files = [f for f in extracted_files if not f.endswith('/')]
        
        if len(extracted_files) != 1:
            raise ValueError(f"Zip file should contain exactly one data file. Found: {extracted_files}")
        
        return os.path.join(os.path.dirname(self.user_file), extracted_files[0])
        
        
    def identify_vendor(self) -> str:
    
        with open(self.user_file, "r", errors='ignore') as f:
            for line in f:
                if not line.startswith("#"):
                    break
                for vendor_name, substring in self.accepted_vendors_dict.items():
                    if substring.lower() in line.lower():
                        
                        return vendor_name
        
        raise ValueError("Vendor could not be identified from the file header.")


    def identify_genome_build(self) -> str:
        
        with open(self.user_file, "r", errors='ignore') as f:
            for line in f:
                if not line.startswith("#"):
                    break
                for build_name, build_code in self.genome_build_dict.items():
                    if build_name.lower() in line.lower():
                        
                        return build_code
    
        raise ValueError("Genome build could not be identified from the file header.")
    
    
    
    def normalize_file(self) -> pd.DataFrame:
            
            
            #'#rsidchromosomepositiongenotype'
            normalized_header = ['# rsid', 'chromosome', 'position', 'genotype']
            allowed_headers = [
                '#rsidchromosomepositiongenotype', #23andme
                'rsidchromosomepositiongenotype', #23andme without # and livingdna
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
                
                data.rename(columns={'rsid': '# rsid'}, inplace=True)

            elif header == 'rsidchromosomepositionallele1allele2': #is ancestry
                data.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2']
                #combine allele1 and allele2 into genotype
                data['genotype'] = data['allele1'] + data['allele2']
                data = data.drop(columns=['allele1', 'allele2'])
            
            elif header == 'RSID,CHROMOSOME,POSITION,RESULT': #is ftdna
                data= data.rename(columns={'RESULT': 'genotype'})
                data.columns = [col.lower() for col in data.columns]
            
            elif header == 'Name,Variation,Chromosome,Position,Strand,YourCode': #bs company

                m = data["Strand"].eq("-")

                trans = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"})

                data.loc[m, "YourCode"] = (
                    data.loc[m, "YourCode"]
                        .astype("string")          
                        .str.translate(trans) 
                )
   
                
                
                data.drop(columns=['Strand', 'Variation'], inplace=True)
                data.rename(columns={'Position': 'position', 'Chromosome': 'chromosome', 'YourCode': 'genotype', 'Name': '# rsid'}, inplace=True)
                
            
            
            #keep only clean genotypes
            data = data[data['genotype'].str.match("^[ATCG]{1,2}$")]

            # clean rsids
            
            data = data.rename(columns={'rsid': '# rsid'})

            data = data[data['# rsid'].str.match("^rs[0-9]{1,}$")]


            #last step to ensure correct column order
            data = data[['# rsid','chromosome','position','genotype']]
            
            #only keep numeric chromosmes
            
            data = data[data['chromosome'].str.match("^[0-9]{1,2}$")]
            
            #rename # rsid to rsid
            data = data.rename(columns={'# rsid': 'rsid'})


            return data
        
    
    
    def extract_header(self) -> str:
        """
        Find the first line that contains an rsxxxx pattern (rs followed by digits)
        and return the line immediately above it with all whitespace removed.
        """
        pattern = re.compile(r'\brs\d+\b', flags=re.IGNORECASE)

        with open(self.user_file, 'r') as file:
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
            with open(self.user_file, 'r') as file:
                first_line = file.readline()
                has_hash = first_line.startswith('#')
            
            
            if has_hash:
                df = pd.read_csv(self.user_file, sep=r'\s+', comment="#", dtype=str, low_memory=False)
            else:
                #no comment, has header
                df = pd.read_csv(self.user_file, sep=r'\s+', dtype=str, low_memory=False)
            
            #check number of columns
            if df.shape[1] < 4:

                    df = pd.read_csv(
                        self.user_file,
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
                    self.user_file,
                    comment="#",          # skip the header blurb
                    sep=",",              # CSV
                    quotechar='"',        # quoted fields like "rs4477212"
                    dtype=str,            # keep everything as strings
                    keep_default_na=False,# don't convert strings like "NA" to NaN
                    engine="c",           # fast
                )

        return df
    
