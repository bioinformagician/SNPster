from typing import Optional
from dataclasses import dataclass
import pandas as pd
import re
import os 
import zipfile

@dataclass
class DataContainer:
    imputation_id: Optional[str] = None
    vendor: Optional[str] = None
    genome_build: Optional[str] = None
    is_forward_strand: Optional[bool] = None
    lift_over: bool = False
    microarray_data: Optional[pd.DataFrame] = None
    split_microarray_data: Optional[pd.DataFrame] = None
    reference_data: Optional[pd.DataFrame] = None
    harmonized_data: Optional[pd.DataFrame] = None
    harmonization_stats: Optional[pd.DataFrame] = None

    

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
        
        # Filter to keep only standard chromosomes (1-22, X, Y, MT)
        # This removes alternate contigs like chr14_GL000009v2_random
        standard_chroms_pattern = "^([0-9]{1,2}|X|Y|MT)$"
        df = df[df["chromosome"].str.match(standard_chroms_pattern, na=False)]
        self.microarray_data = df
        
        print(df)
        
        
        #------------ Harmonization methods ------------
        
        
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
        
        # Filter out alternate contigs, patches, and unplaced sequences
        # Keep only standard chromosomes: 1-22, X, Y, MT
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']

        initial_count = len(self.harmonized_data)
        self.harmonized_data['chrom_clean'] = self.harmonized_data['chromosome'].astype(str).str.replace('chr', '')
        self.harmonized_data = self.harmonized_data[self.harmonized_data['chrom_clean'].isin(valid_chromosomes)]
        self.harmonized_data = self.harmonized_data.drop(columns=['chrom_clean'])
        filtered_count = initial_count - len(self.harmonized_data)
        
        if filtered_count > 0:
            print(f"Filtered out {filtered_count} variants on non-standard chromosomes (alt contigs, patches, etc.)")
        
        # Keep only the 4 columns needed for PLINK --23file format
        self.harmonized_data = self.harmonized_data[['# rsid', 'chromosome', 'position', 'genotype']]
        
        


class FileHandler:

    def __init__(self, 
                 user_file: Optional[str],
                 accepted_vendors_dict: dict[str, str],
                 genome_build_dict: dict[str, str],
                 forward_strand_vendors: list[str]
                 ):
        
        self.user_file = user_file
        self.accepted_vendors_dict = accepted_vendors_dict
        self.genome_build_dict = genome_build_dict
        self.forward_strand_vendors = forward_strand_vendors
    
    def is_zipped_file(self) -> bool:
        return zipfile.is_zipfile(self.user_file)
    
    def unzip_file(self) -> str:

        with zipfile.ZipFile(self.user_file, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(self.user_file))
            extracted_files = zip_ref.namelist()
        
        # Filter out directories and metadata artifacts (e.g. __MACOSX, .DS_Store, AppleDouble)
        cleaned_files = []
        for file_name in extracted_files:
            if file_name.endswith('/'):
                continue

            path_parts = file_name.split('/')
            base_name = os.path.basename(file_name)

            if '__MACOSX' in path_parts:
                continue
            if base_name.startswith('._') or base_name.startswith('.'):
                continue

            cleaned_files.append(file_name)

        extracted_files = cleaned_files
        
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
    
