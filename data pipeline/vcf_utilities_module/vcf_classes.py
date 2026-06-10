import os
from pathlib import Path
import subprocess
import tempfile
import shutil
import pandas as pd
import gzip

class VCFEnvironmentHandler:
    def __init__(self, 
                 vcf_samplesheet_path: str | None = None, #this csv file contains columns: full_vcf_path, chrom, imputation_id, sample_id
                 output_dir: str | None = None,
                 vcf_file: str | None = None
                 ):
        
        self.vcf_samplesheet_path = vcf_samplesheet_path
        self.output_dir = output_dir
        self.vcf_file = vcf_file
    
    
    
    
class VCFHandler: #split this into a VCFMerger and VCFSplitter class later, but for now it is easier to keep bcftools merging and splitting code together since they both use the same sample sheet format and are run sequentially in the workflow
    def __init__(
        self,
        vcf_environment_handler: VCFEnvironmentHandler,
    ):
        self.vcf_environment_handler = vcf_environment_handler
        self.vcf_utilities = VCFUtilities()


    def _get_imputation_id_from_sample_id(self, sample_id) -> str:
        return self.vcf_utilities._get_imputation_id_from_sample_id(sample_id)


    def _get_sample_ids_from_vcf(self, vcf_file_path: str) -> list[str]:
        return self.vcf_utilities._get_sample_ids_from_vcf(vcf_file_path)


    def _get_imputation_id_from_vcf(self, vcf_file_path: str) -> str:
        return self.vcf_utilities._get_imputation_id_from_vcf(vcf_file_path)


    def _get_chromosome_from_vcf(self, vcf_file_path: str) -> str:
        return self.vcf_utilities._get_chromosome_from_vcf(vcf_file_path)


    def _prepare_merge_input(self, input_path: str, tmp_dir: str) -> str:
        # Always convert to fresh BGZF in temp space and merge with --no-index.
        # This avoids failures from plain gzip files mislabeled as .vcf.gz.
        # Also drop records with ALT='.' because malformed files can contain
        # non-reference genotypes at such rows, which causes bcftools merge to fail.
        file_name = os.path.basename(input_path)
        if file_name.endswith(".vcf.gz"):
            file_name = file_name[:-7]
        elif file_name.endswith(".vcf"):
            file_name = file_name[:-4]

        prepared_path = os.path.join(tmp_dir, f"{file_name}.bgzf.vcf.gz")
        subprocess.run([
            "bcftools", "view",
            "-e", "ALT='.'",
            input_path,
            "-Oz", "-o", prepared_path,
        ], check=True)

        return prepared_path

    def merge_vcf_files(self) -> None:

        """columns in sample_sheet: full_vcf_path,chrom there will only be data from a single chrom when parsed to this module"""
        merged_file_sample_sheet = pd.read_csv(self.vcf_environment_handler.vcf_samplesheet_path, dtype=str)

        chrom = merged_file_sample_sheet["chrom"].iloc[0] #this has the format "chr1", "chr2"...
        combined_vcf_path_string = merged_file_sample_sheet["full_vcf_path"].str.cat(sep=",")


        print(f"VCF files for chromosome {chrom}: {combined_vcf_path_string}")
        os.makedirs(self.vcf_environment_handler.output_dir, exist_ok=True)
        output_path = os.path.join(self.vcf_environment_handler.output_dir, f"{chrom}.merged.vcf.gz")

        tmp_dir = tempfile.mkdtemp()
        try:
            prepared_files = [self._prepare_merge_input(path, tmp_dir) for path in merged_file_sample_sheet["full_vcf_path"].tolist()]
            subprocess.run(["bcftools", "merge", "--no-index", *prepared_files, "-Oz", "-o", output_path], check=True)
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
        is_gzipped = input_vcf.endswith(".gz")
        input_opener = gzip.open if is_gzipped else open
        contig_line = f"##contig=<ID={chrom}>\n"

        with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as plain_tmp:
            plain_tmp_path = plain_tmp.name

        try:
            with input_opener(input_vcf, "rt") as f:
                lines = f.readlines()

            has_contig = any(line.startswith(f"##contig=<ID={chrom}") for line in lines)

            with open(plain_tmp_path, "wt") as f:
                for line in lines:
                    if line.startswith("#CHROM") and not has_contig:
                        f.write(contig_line)
                    f.write(line)

            if is_gzipped:
                subprocess.run([
                    "bcftools", "view",
                    plain_tmp_path,
                    "-Oz", "-o", output_vcf
                ], check=True)
            else:
                os.replace(plain_tmp_path, output_vcf)
                plain_tmp_path = None
        finally:
            if plain_tmp_path is not None and os.path.exists(plain_tmp_path):
                os.remove(plain_tmp_path)
    

    
    
    def split_vcf_files(self) -> None:
        """Splits the merged VCF file into chromosome-specific files.
        
            columns in sample_sheet: full_vcf_path, chrom
        """
        
        vcf_file = self.vcf_environment_handler.vcf_file
        os.makedirs(self.vcf_environment_handler.output_dir, exist_ok=True)
        chrom = self._get_chromosome_from_vcf(vcf_file)
        self._add_minimal_contig_header(vcf_file, vcf_file, chrom)  # in-place to add contig header, otherwise bcftools will not accept file format
        
        print(f"Splitting vcf file using bcftools, outputting to {self.vcf_environment_handler.output_dir}...")
        
        subprocess.run([
            "bcftools", "+split",
            "-Oz",
            "-o", self.vcf_environment_handler.output_dir, vcf_file
        ], check=True)
        
        outputted_files = [f for f in os.listdir(self.vcf_environment_handler.output_dir) if f.endswith(".vcf.gz")]
        
        #rename files from sample id to {IMPIDx}.{chrx}.{split}.vcf.gz
        
        for file in outputted_files:
            sample_id = file.split(".vcf.gz")[0]
            imputation_id = self._get_imputation_id_from_sample_id(sample_id)
            new_file_name = f"IMPID{imputation_id}.chr{chrom}.split.vcf.gz"
            os.rename(os.path.join(self.vcf_environment_handler.output_dir, file), 
                      os.path.join(self.vcf_environment_handler.output_dir, new_file_name))
        
    
    
    
    def run_vcf_splitting(self, samplesheet_dir: str) -> None:
        
        
        print("Running VCF splitting with Nextflow...")
        
        command = [
            "nextflow", "run", "vcf_splitter_pipeline.nf",
            "--samplsheet_dir", samplesheet_dir,
            "--output_dir", samplesheet_dir
        ]
        
        print(f"running command {command}")
        
        subprocess.run(command, check=True)




class VCFUtilities:
    
    """Class to provide utilities such as adding and reading metadata to/from VCF files """
    
    IMPUTATION_ID_KEY = "IMPUTATION_ID"

    def _rewrite_vcf_preserving_bgzf(self, vcf_path: Path, line_rewriter) -> None:
        """Rewrite VCF text while keeping BGZF framing for .vcf.gz files."""
        is_gzipped = self._is_gzipped(vcf_path)
        input_opener = gzip.open if is_gzipped else open

        with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as plain_tmp:
            plain_tmp_path = Path(plain_tmp.name)

        try:
            with input_opener(vcf_path, "rt") as fin, open(plain_tmp_path, "wt") as fout:
                line_rewriter(fin, fout)

            if is_gzipped:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz") as gz_tmp:
                    gz_tmp_path = Path(gz_tmp.name)
                try:
                    subprocess.run([
                        "bcftools", "view",
                        str(plain_tmp_path),
                        "-Oz", "-o", str(gz_tmp_path)
                    ], check=True)
                    os.replace(gz_tmp_path, vcf_path)
                except Exception:
                    gz_tmp_path.unlink(missing_ok=True)
                    raise
            else:
                os.replace(plain_tmp_path, vcf_path)
                plain_tmp_path = None
        finally:
            if plain_tmp_path is not None:
                plain_tmp_path.unlink(missing_ok=True)
    
    def _get_imputation_id_from_sample_id(self, sample_id) -> str:
        return sample_id.split("_", 1)[1]


    def _get_sample_ids_from_vcf(self, vcf_file_path: str) -> list[str]:
        """Extract sample IDs from the #CHROM header line of a VCF/VCF.GZ file."""
        opener = gzip.open if str(vcf_file_path).endswith(".gz") else open

        with opener(vcf_file_path, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    columns = line.rstrip("\n").split("\t")
                    sample_ids = columns[9:]
                    if not sample_ids:
                        raise ValueError(f"No sample IDs found in VCF header: {vcf_file_path}")
                    return sample_ids

        raise ValueError(f"No #CHROM header found in VCF file: {vcf_file_path}")
    
    def _get_imputation_id_from_vcf(self, vcf_file_path: str) -> str:
        sample_ids = self._get_sample_ids_from_vcf(vcf_file_path)
        if not sample_ids:
            raise ValueError(f"No sample IDs found in VCF file: {vcf_file_path}")
        return self._get_imputation_id_from_sample_id(sample_ids[0])


    def _get_chromosome_from_vcf(self, vcf_file_path: str) -> str:
        """Extract chromosome value from the first variant row in a VCF/VCF.GZ file."""
        opener = gzip.open if str(vcf_file_path).endswith(".gz") else open

        with opener(vcf_file_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if not fields or not fields[0]:
                    continue
                return fields[0]

        raise ValueError(f"No variant rows found in VCF file: {vcf_file_path}")

    def add_imputation_id_to_vcf(self, vcf_file: str, imputation_id: list) -> None:
        """
        Adds imputation ID metadata to the VCF header in place.

        If the VCF already has:
            ##IMPUTATION_ID=12,93

        and imputation_id is:
            [45, 60]

        it becomes:
            ##IMPUTATION_ID=12,93,45,60
        """

        vcf_path = Path(vcf_file)
        new_ids = [int(value) for value in imputation_id]
        existing_ids = self.read_imputation_id_from_vcf(vcf_file)

        merged_ids = existing_ids + [
            value for value in new_ids
            if value not in existing_ids
        ]

        metadata_line = (
            f"##{self.IMPUTATION_ID_KEY}="
            f"{','.join(str(value) for value in merged_ids)}\n"
        )

        def line_rewriter(fin, fout):
            inserted = False
            for line in fin:
                if line.startswith(f"##{self.IMPUTATION_ID_KEY}="):
                    if not inserted:
                        fout.write(metadata_line)
                        inserted = True
                    continue

                if line.startswith("#CHROM") and not inserted:
                    fout.write(metadata_line)
                    inserted = True

                fout.write(line)

        self._rewrite_vcf_preserving_bgzf(vcf_path, line_rewriter)
    
    
    def read_imputation_id_from_vcf(self, vcf_file: str) -> list[int]:
        """
        Reads imputation IDs from the VCF header.

        """

        vcf_path = Path(vcf_file)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file does not exist: {vcf_file}")
        
        opener = gzip.open if self._is_gzipped(vcf_path) else open

        with opener(vcf_path, "rt") as f:
            for line in f:
                line = line.rstrip("\n")

                if line.startswith(f"##{self.IMPUTATION_ID_KEY}="):
                    value = line.split("=", 1)[1]

                    return [
                        int(item)
                        for item in value.split(",")
                        if item
                    ]

                if line.startswith("#CHROM"):
                    break

        return []
    
    
    def get_chromosome_from_vcf(self, vcf_file: str) -> str | None:
        """Reads the chromosome from the first data line in the VCF."""
        vcf_path = Path(vcf_file)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file does not exist: {vcf_file}")
        
        opener = gzip.open if self._is_gzipped(vcf_path) else open

        with opener(vcf_path, "rt") as f:
            for line in f:
                line = line.rstrip("\n")
                
                # Skip header lines
                if line.startswith("#"):
                    continue
                
                # First data line - chromosome is first column
                if line:
                    chrom = line.split("\t")[0]
                    return chrom

        return None

    @staticmethod
    def _is_gzipped(path: Path) -> bool:
        return path.suffix == ".gz"
    
    
    def make_vcf_merging_samplesheet(self, vcf_file_dir: str, output_dir: str) -> str:
        """ Take a dir of vcf files and output a csv file in that dir with the columns full_vcf_path and chrom, which can be used as input for the vcf merging process in the workflow
            This function will produce 22 csv files, one for each chromosome
        """
        
        vcf_files = [f for f in os.listdir(vcf_file_dir) if f.endswith(".vcf.gz")]
        full_vcf_paths = [os.path.join(vcf_file_dir, f) for f in vcf_files]
        chromosomes = [self.get_chromosome_from_vcf(path) for path in full_vcf_paths]
        

        samplesheet_df = pd.DataFrame({
            "full_vcf_path": full_vcf_paths,
            "chrom": chromosomes
        })
        
        for chrom in samplesheet_df["chrom"].unique():
            chrom_df = samplesheet_df[samplesheet_df["chrom"] == chrom]
            output_path = os.path.join(output_dir, f"vcf_merge_sheet_chr{chrom}.csv")
            chrom_df.to_csv(output_path, index=False)
    

    