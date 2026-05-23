import os
import subprocess
import tempfile
import shutil
import pandas as pd
import gzip

class VCFEnvironmentHandler:
    def __init__(self, 
                 vcf_samplesheet_path: str | None = None, #this csv file contains columns: full_vcf_path, chrom, imputation_id, sample_id
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


    def _prepare_merge_input(self, input_path: str, tmp_dir: str) -> str:
        if input_path.endswith(".vcf.gz"):
            if not os.path.exists(f"{input_path}.tbi"):
                subprocess.run(["bcftools", "index", "-t", input_path], check=True)
            return input_path

        file_name = os.path.basename(input_path)
        prepared_path = os.path.join(tmp_dir, f"{file_name}.gz")

        subprocess.run(["bcftools", "view", input_path, "-Oz", "-o", prepared_path], check=True)
        subprocess.run(["bcftools", "index", "-t", prepared_path], check=True)

        return prepared_path

    def merge_vcf_files(self) -> None:

        """columns in sample_sheet: full_vcf_path,chrom there will only be data from a single chrom when parsed to this module"""
        merged_file_sample_sheet = pd.read_csv(self.vcf_environment_handler.vcf_samplesheet_path, dtype=str)

        chrom = merged_file_sample_sheet["chrom"].iloc[0]
        combined_vcf_path_string = merged_file_sample_sheet["full_vcf_path"].str.cat(sep=",")


        print(f"VCF files for chromosome {chrom}: {combined_vcf_path_string}")
        output_path = os.path.join(self.vcf_environment_handler.output_dir, f"merged_vcf_chr_{chrom}.vcf.gz")

        tmp_dir = tempfile.mkdtemp()
        try:
            prepared_files = [self._prepare_merge_input(path, tmp_dir) for path in merged_file_sample_sheet["full_vcf_path"].tolist()]
            subprocess.run(["bcftools", "merge", *prepared_files, "-Oz", "-o", output_path], check=True)
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    
    
    
    def run_vcf_merging(self, merged_samplesheet_dir: str) -> None: #this is duplicate code from pgs_calc module, fix later
        """nextflow run vcf_combiner_pipeline.nf --samplsheet_dir ./nf_data --output_dir ./nf_data"""
        
        print("Running VCF merging with Nextflow...")
        
        command = [
            "nextflow", "run", "vcf_combiner_pipeline.nf",
            "--samplsheet_dir", merged_samplesheet_dir,
            "--output_dir", merged_samplesheet_dir
        ]
        
        print(f"running command {command}")
        
        subprocess.run(command, check=True)
    
    
    
    
    
    
    
    
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
    
    def _get_imputation_id_from_sample_id(self, sample_id) -> str:
        return sample_id.split("_", 1)[1]
    
    
    def split_vcf_files(self) -> pd.DataFrame:
        """Splits the merged VCF file into chromosome-specific files.
        
            columns in sample_sheet: full_vcf_path, chrom, sample_id
        """
        
        imputed_samplesheet = pd.read_csv(self.vcf_environment_handler.vcf_samplesheet_path, dtype=str)

        vcf_file = imputed_samplesheet["full_vcf_path"].iloc[0]
        chrom = imputed_samplesheet["chrom"].iloc[0]
        sample_ids = imputed_samplesheet["sample_id"].iloc[0].split(",")
        self._add_minimal_contig_header(vcf_file, vcf_file, chrom)  # in-place to add contig header, otherwise bcftools will not accept file format
        
        output_df = pd.DataFrame()
        
        for sample_id in sample_ids:
            
            output_path = os.path.join(self.vcf_environment_handler.output_dir, f"split_imputed_ImpID{self._get_imputation_id_from_sample_id(sample_id)}_chr{chrom}.vcf.gz")
            imputation_id = self._get_imputation_id_from_sample_id(sample_id)
            
            print(f"Extracting sample {sample_id} (Imputation ID: {imputation_id}) from {vcf_file} to {output_path} using bcftools...")
            
            subprocess.run([
                "bcftools", "view",
                "-s", sample_id,
                "-Oz",
                "-o", output_path,
                vcf_file
            ], check=True)
            
            df = pd.DataFrame({
                "chrom": chrom,
                "imputation_id": imputation_id,
                "split_vcf_filepath": output_path,
                "qc_imputed_file": None
            }, index=[0])
            
            output_df = pd.concat([output_df, df], ignore_index=True)
    
        return output_df
    
    
    
    def run_vcf_splitting(self, samplesheet_dir: str) -> None:
        
        
        print("Running VCF splitting with Nextflow...")
        
        command = [
            "nextflow", "run", "vcf_splitter_pipeline.nf",
            "--samplsheet_dir", samplesheet_dir,
            "--output_dir", samplesheet_dir
        ]
        
        print(f"running command {command}")
        
        subprocess.run(command, check=True)

    