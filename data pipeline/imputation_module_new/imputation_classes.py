import os
import pandas as pd
from dataclasses import dataclass
import gzip
import subprocess
import numpy as np
import re
import gc
import polars as pl




class EnvironmentHandler:
    def __init__(self,
                 working_dir: str,
                 java_exe: str,
                 beagle_jar: str,
                 heap_gb: int,
                 beagle_threads: int,
                 output_dir: str,
                 beagle_reference_dir: str,
                 plink_map_dir: str,
                 vcf_file: str = None,
                 beagle_references: dict[str, str] = None,
                 plink_map_files: dict[str, str] = None,
                 vcf_plink_reference_mapping: pd.DataFrame = None #columns: chromosome_number, vcf_file, reference_file, plink_map_file
                 ):
        
        self.working_dir = working_dir
        self.java_exe = java_exe
        self.beagle_jar = beagle_jar
        self.vcf_plink_reference_mapping = vcf_plink_reference_mapping
        self.heap_gb = heap_gb
        self.beagle_threads = beagle_threads
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.beagle_reference_dir = beagle_reference_dir
        self.plink_map_dir = plink_map_dir
        self.plink_map_files = plink_map_files
        self.beagle_references = beagle_references
        self.validate_paths()
        self.set_beagle_files()
        self.set_plink_map_files()


        
    def set_beagle_files(self) -> dict:
        beagle_files = {}
        for file in os.listdir(self.beagle_reference_dir):
            if file.endswith(".vcf.gz") or file.endswith(".bref3"):
                chrom = re.search(r"chr(\d+)\.", file)
                if chrom:
                    beagle_files[chrom.group(1)] = os.path.join(self.beagle_reference_dir, file)
            
        self.beagle_references = beagle_files
    
    def set_plink_map_files(self) -> dict:
        plink_map_files = {}
        for file in os.listdir(self.plink_map_dir):
            if file.endswith(".map"):
                chrom = re.search(r"chr(\d+)[^/]*\.map$", file)
                if chrom:
                    plink_map_files[chrom.group(1)] = os.path.join(self.plink_map_dir, file)
        
        self.plink_map_files = plink_map_files


    def validate_paths(self) -> None:
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        for path in [
            self.working_dir,
            self.java_exe,
            self.beagle_jar,
            self.beagle_reference_dir,
            self.plink_map_dir,
        ]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")














class ImputationRunner:
    def __init__(self,
                    environment_handler: EnvironmentHandler,
                    ):
        
        self.environment_handler = environment_handler
        
    
    def _get_imputation_id_from_vcf(self, vcf_path: str) -> str:
        filename = os.path.basename(vcf_path)
        match = re.search(r"_(\d+)\.vcf\.gz$", filename)
        if match:
            return match.group(1)
        else:
            raise ValueError(f"Could not extract imputation ID from VCF filename: {filename}")
    
    
    
    def create_vcf_reference_mapping(self) -> None:
        vcf_file = self.environment_handler.vcf_file
        if not vcf_file:
            raise ValueError("No input VCF provided. Expected environment_handler.vcf_file to be set.")

        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"Input VCF does not exist: {vcf_file}")

        chrom_value = self._get_chrom_number_from_vcf(vcf_file)
        chrom_match = re.search(r"(\d+)", str(chrom_value))
        if chrom_match is None:
            raise ValueError(f"Could not extract chromosome number from VCF content: {vcf_file}")

        vcf_files_df = pd.DataFrame([
            {"chromosome_number": chrom_match.group(1), "vcf_file": vcf_file}
        ])
        
        beagle_reference_df=pd.DataFrame(
            list(self.environment_handler.beagle_references.items()),
            columns=["chromosome_number", "reference_file"]
        )
        
        plink_map_files=pd.DataFrame(
            list(self.environment_handler.plink_map_files.items()),
            columns=["chromosome_number", "plink_map_file"]
        )
        
        #join vcf and mapping_df
        mapping_df = pd.merge(
            vcf_files_df,
            beagle_reference_df,
            on="chromosome_number",
            how="inner"
        )
        
        mapping_df = pd.merge(
            mapping_df,
            plink_map_files,
            on="chromosome_number",
            how="inner"
        )
        
        if mapping_df.empty:
            raise ValueError("File mismatch when merging VCF, Beagle reference, and PLINK map files.")
    
        #add imputed_file
        mapping_df["imputed_file"] = None
        self.environment_handler.vcf_plink_reference_mapping = mapping_df
    
    
                
        
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



    def impute_data(self, gt_vcf: str, map_file: str,
                        ref_panel: str, chr_number) -> None:
        

        out = os.path.join(f"{self.environment_handler.output_dir}/imputed_chr{chr_number}.risk")
        
        cmd = [
            self.environment_handler.java_exe, f"-Xmx{self.environment_handler.heap_gb}g", "-jar", str(self.environment_handler.beagle_jar),
            f"gt={gt_vcf}",
            f"ref={ref_panel}",
            f"map={map_file}",
            f"out={out}",
            f"nthreads={self.environment_handler.beagle_threads}",
            f"gp=true",
            # Beagle imputes by default when ref= is provided
        ]
        
        self.run_command(cmd)
        
        out = out + ".vcf.gz"
        
        return out
    
        
    
    def impute_vcf_file(self) -> None:

        
        map_df = self.environment_handler.vcf_plink_reference_mapping
        
        vcf_path = map_df["vcf_file"].iloc[0]
        map_file = map_df["plink_map_file"].iloc[0]
        ref_file = map_df["reference_file"].iloc[0]
        chr_number = map_df["chromosome_number"].iloc[0]
        
        print(f"Imputing chromosome {chr_number}...")
        
        self.impute_data(vcf_path, map_file, ref_file, chr_number)
        
            

        
        
    
    def create_samplesheet(self) -> None:
        
        """Write a sample sheet for each imputation id + chromosome number, with the columns: sampleset, path_prefix, chrom, format. 
        This format is needed for the nextflow pgs_calc module downstream """        

        
        for imputation_id, group in self.environment_handler.split_imputed_vcf_mapping.groupby("imputation_id"):
        
            vcf_paths= list(group["split_vcf_filepath"])
            
            vcf_filenames = [os.path.basename(path) for path in vcf_paths]
            
            vcf_base_names = [name.split(".", 1)[0] for name in vcf_filenames]
            
            
            chr_numbers = list(group["chromosome_number"])

            samplesheet_df = pd.DataFrame({
                "sampleset": "cohort",
                "path_prefix": vcf_base_names,
                "chrom": chr_numbers,
                "format": "vcf",
            })
            
            samplesheet_df.to_csv(
                os.path.join(self.environment_handler.output_dir, f"samplesheet_IMPID{imputation_id}.csv"),
                index=False
            )
            
            
        
    
    def _open_text(self, vcf_file):
        vcf_file = str(vcf_file)
        return gzip.open(vcf_file, "rt") if vcf_file.endswith(".gz") else open(vcf_file, "rt")


    def _get_sample_ids(self, vcf_file) -> list[str]:
        """takes a vcf file path as input at returns the sample ids in that file(will be more than one if the file is a merged vcf file)"""
        with self._open_text(vcf_file) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    columns = line.rstrip("\n").split("\t")
                    return columns[9:]

        raise ValueError(f"No #CHROM header found in VCF file: {vcf_file}")


    def _get_chrom_number_from_vcf(self, vcf_file) -> str:
        with self._open_text(vcf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                columns = line.rstrip("\n").split("\t")
                return columns[0]

        raise ValueError(f"No variant rows found in VCF file: {vcf_file}")


    def _get_imputation_id_from_sample_id(self, sample_id) -> str:
        return sample_id.split("_", 1)[1]


    def _add_minimal_contig_header(self, input_vcf: str, output_vcf: str, chrom:str) -> None:
        
        """BCF tools require the contig header to be present in the VCF so need to add it before splitting vcf files by sample"""
        
        with gzip.open(input_vcf, "rt") as f:
            lines = f.readlines()

        contig_line = f"##contig=<ID={chrom}>\n"
        has_contig = any(line.startswith(f"##contig=<ID={chrom}") for line in lines)

        with gzip.open(output_vcf, "wt") as f:
            for line in lines:
                if line.startswith("#CHROM") and not has_contig:
                    f.write(contig_line)
                f.write(line)




    def _get_imputation_id_from_split_imputed_filename(self, filename) -> str:
        match = re.search(r"split_imputed_ImpID(\d+)_chr\d+\.vcf\.gz$", filename)
        if match:
            return match.group(1)
        else:
            raise ValueError(f"Could not extract imputation ID from split imputed VCF filename: {filename}")
    

            
            
        

    

