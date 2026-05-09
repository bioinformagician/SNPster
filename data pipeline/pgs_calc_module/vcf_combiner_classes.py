import os
import subprocess
import pandas as pd

class VCFEnvironmentHandler:
    def __init__(self, 
                 vcf_samplesheet_path: str | None = None, #this csv file contains columns: sampleset, path_prefix, chrom, format
                 output_dir: str | None = None,
                 ):
        
        self.vcf_samplesheet_path = vcf_samplesheet_path
        self.output_dir = output_dir
    
    
    
    
class VCFHandler:
    def __init__(
        self,
        vcf_environment_handler: VCFEnvironmentHandler,
    ):
        self.vcf_environment_handler = vcf_environment_handler


    def _prepare_merge_input(self, input_path: str) -> str:
        prepared_dir = os.path.join(self.vcf_environment_handler.output_dir, "merge_inputs")
        os.makedirs(prepared_dir, exist_ok=True)

        if input_path.endswith(".vcf.gz"):
            if not os.path.exists(f"{input_path}.tbi"):
                subprocess.run(["bcftools", "index", "-t", input_path], check=True)
            return input_path

        file_name = os.path.basename(input_path)
        prepared_path = os.path.join(prepared_dir, f"{file_name}.gz")

        if not os.path.exists(prepared_path):
            subprocess.run(["bcftools", "view", input_path, "-Oz", "-o", prepared_path], check=True)

        if not os.path.exists(f"{prepared_path}.tbi"):
            subprocess.run(["bcftools", "index", "-t", prepared_path], check=True)

        return prepared_path

    def merge_vcf_files(self) -> None:

        """columns in sample_sheet: sampleset,full_vcf_path,chrom,format there will only be data from a single chrom when parsed to this module"""
        merged_file_sample_sheet = pd.read_csv(self.vcf_environment_handler.vcf_samplesheet_path, dtype=str)

        chrom = merged_file_sample_sheet["chrom"].iloc[0]
        combined_vcf_path_string = merged_file_sample_sheet["full_vcf_path"].str.cat(sep=",")
        sampleset = merged_file_sample_sheet["sampleset"].iloc[0]


        print(f"VCF files for chromosome {chrom}: {combined_vcf_path_string}")
        output_path = os.path.join(self.vcf_environment_handler.output_dir, f"merged_vcf_chr_{chrom}.vcf.gz")

        prepared_files = [self._prepare_merge_input(path) for path in merged_file_sample_sheet["full_vcf_path"].tolist()]


        subprocess.run(["bcftools", "merge", *prepared_files, "-Oz", "-o", output_path], check=True)
