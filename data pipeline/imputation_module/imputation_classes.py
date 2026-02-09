import os
import pandas as pd
from dataclasses import dataclass
import gzip
import subprocess
import numpy as np
import re
import gc
import polars as pl


@dataclass(frozen=True)
class QCThresholds:
    gp_min: float
    ds_tol: float
    snps_only: bool
    biallelic_only: bool



@dataclass
class DataContainer:
    qc_thresholds: QCThresholds = QCThresholds(None, None, False, False)
    imputed_data: pd.DataFrame = None
    qc_status: bool = False
    user_id: str = None
    
    

    def qc_imputed_data(self) -> None:
        print("starting qc of imputed data...")
        
        x = self.imputed_data
        n = len(x)

        if n == 0:
            print("0.0% of input variants passed qc requirements (no variants)")
            self.qc_status = False
            return

        # ---- build COMPLETE mask for all filters ----
        keep_mask = np.ones(n, dtype=bool)

        # SNP-only: ref and alt length == 1
        if self.qc_thresholds.snps_only:
            ref_len = x["REF"].str.len().to_numpy()
            alt_len = x["ALT"].str.len().to_numpy()
            keep_mask &= (ref_len == 1) & (alt_len == 1)

        # biallelic-only: no ',' in alt
        if self.qc_thresholds.biallelic_only:
            has_comma = x["ALT"].str.contains(",").to_numpy()
            keep_mask &= ~has_comma

        # GP and genotype filters - compute on original data
        gp = x[["GP_00", "GP_01", "GP_11"]].to_numpy()
        gpmax = gp.max(axis=1)
        
        gt = x["GT"].astype(str).to_numpy()
        is_hom_ref = (gt == "0|0") | (gt == "0/0")
        is_het     = (gt == "0|1") | (gt == "1|0") | (gt == "0/1") | (gt == "1/0")
        is_hom_alt = (gt == "1|1") | (gt == "1/1")

        ds = x["DS"].to_numpy(dtype=float)
        gp01 = x["GP_01"].to_numpy(dtype=float)

        hom_ref_ok = is_hom_ref & (ds <= self.qc_thresholds.ds_tol)
        het_ok     = is_het & (gp01 >= self.qc_thresholds.gp_min) & (ds >= 1 - self.qc_thresholds.ds_tol) & (ds <= 1 + self.qc_thresholds.ds_tol)
        hom_alt_ok = is_hom_alt & (ds >= 2 - self.qc_thresholds.ds_tol)

        # Combine ALL filters into one mask
        keep_mask &= (gpmax >= self.qc_thresholds.gp_min) & (hom_ref_ok | het_ok | hom_alt_ok)

        # Apply filter ONCE
        kept = keep_mask.sum()
        pct = (kept / n * 100.0) if n else 0.0
        
        if kept == 0:
            print("0.0% of input variants passed qc requirements")
            self.qc_status = False
            return
        
        print(f"{pct:.1f}% of input variants passed qc requirements")

        # Single filter operation - only ONE new DataFrame created
        self.imputed_data = x.loc[keep_mask]
        self.qc_status = True

        print("QC of imputed data completed")
    
    
    

    
    

class EnvironmentHandler:
    def __init__(self,
                 working_dir: str,
                 java_exe: str,
                 beagle_jar: str,
                 heap_gb: int,
                 threads: int,
                 vcf_files_dir: str,
                 output_dir: str,
                 beagle_reference_dir: str,
                 plink_map_dir: str,
                 vcf_file_paths: dict[str, str] = None,
                 beagle_references: dict[str, str] = None,
                 plink_map_files: dict[str, str] = None,
                 imputed_dir: str = None,
                 imputed_files: dict[str,str] = None,
                 qc_imputed_files: dict[str,str] = None,
                 vcf_plink_reference_mapping: pd.DataFrame = None
                 ):
        
        self.working_dir = working_dir
        self.java_exe = java_exe
        self.beagle_jar = beagle_jar
        self.vcf_plink_reference_mapping = vcf_plink_reference_mapping
        self.imputed_dir = imputed_dir
        self.heap_gb = heap_gb
        self.threads = threads
        self.imputed_files = imputed_files
        self.qc_imputed_files = qc_imputed_files
        self.output_dir = output_dir
        self.beagle_reference_dir = beagle_reference_dir
        self.plink_map_dir = plink_map_dir
        self.plink_map_files = plink_map_files
        self.beagle_references = beagle_references
        self.vcf_files_dir = vcf_files_dir
        self.vcf_file_paths = vcf_file_paths
        self.validate_paths()
        self.make_imputed_directory()
        self.set_beagle_files()
        self.set_plink_map_files()
        self.set_vcf_files()

        
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
    
    def set_vcf_files(self) -> dict:
        vcf_files = {}
        pattern = re.compile(r"_chr(1?\d|2[0-2])_")

        for file in os.listdir(self.vcf_files_dir):
            if file.endswith(".vcf.gz"):
                match = pattern.search(file)
                if match:
                    chrom = match.group(1)
                    vcf_files[chrom] = os.path.join(self.vcf_files_dir, file)
        
        self.vcf_file_paths = vcf_files
    
    
    def make_imputed_directory(self) -> None:
        imputed_dir = os.path.join(self.working_dir, "imputed")
        os.makedirs(imputed_dir, exist_ok=True)
        self.imputed_dir = imputed_dir
    

    def validate_paths(self) -> None:
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        for path in [
            self.working_dir,
            self.java_exe,
            self.beagle_jar,
            self.beagle_reference_dir,
            self.plink_map_dir,
            self.vcf_files_dir
        ]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")


class WorkflowOrchestrator:
    def __init__(self,
                    environment_handler: EnvironmentHandler,
                    data_container: DataContainer
                    ):
        
        self.environment_handler = environment_handler
        self.data_container = data_container
        self.create_vcf_reference_mapping()
        
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

    def make_result_subdir(self) -> None:
        result_dir = os.path.join(self.environment_handler.output_dir, f"user_id_{self.data_container.user_id}")
        os.makedirs(result_dir, exist_ok=False)
        
        self.environment_handler.output_dir = result_dir
        print(f"Created result directory: {result_dir}")


    def impute_data(self, gt_vcf: str, map_file: str,
                        ref_panel: str, chr_number) -> None:
        

        out = os.path.join(f"{self.environment_handler.imputed_dir}/imputed_chr{chr_number}.risk")
        
        cmd = [
            self.environment_handler.java_exe, f"-Xmx{self.environment_handler.heap_gb}g", "-jar", str(self.environment_handler.beagle_jar),
            f"gt={gt_vcf}",
            f"ref={ref_panel}",
            f"map={map_file}",
            f"out={out}",
            f"nthreads={self.environment_handler.threads}",
            f"gp=true",
            # Beagle imputes by default when ref= is provided
        ]
        
        self.run_command(cmd)
        
        out = out + ".vcf.gz"
        
        return out
        
        
    
    def impute_vcf_files(self) -> None:

        
        for _, row in self.environment_handler.vcf_plink_reference_mapping.iterrows():
            vcf_path = row["vcf_file"]
            map_file = row["plink_map_file"]
            ref_file = row["reference_file"]
            chr_number = row["chromosome_number"]
            
            out = self.impute_data(vcf_path, map_file, ref_file, chr_number)
            
            
            self.environment_handler.imputed_files[chr_number] = out
            
    
    
    def load_vcf_to_df_pandas(self, imputed_file) -> None:
        
        """stats from benchmarking loading and qc'ing 3 chromosomes separately:
                    Wall time: 24.98s
                    Peak RAM:  1332.8 MB
                    Avg CPU:   99.6%
                    Peak CPU:  141.1%
        """
        
        dataframe = pd.read_csv(
            imputed_file,
            sep="\t",
            comment="#",
            compression="gzip",
            header=None
        )
        dataframe.columns = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","FAM001_ID001"]

        # ---- FORMAT: GT / DS / GP ----
        fmt = dataframe["FAM001_ID001"].str.split(":", n=2, expand=True)
        dataframe["GT"] = fmt[0]

        # DS: keep first value only
        dataframe["DS"] = pd.to_numeric(fmt[1].str.split(",", n=1).str[0], errors="coerce")

        gplist = fmt[2].astype(str).str.split(",")
        dataframe["GP_00"] = pd.to_numeric(gplist.str[0], errors="coerce")
        dataframe["GP_01"] = pd.to_numeric(gplist.str[1], errors="coerce")
        dataframe["GP_11"] = pd.to_numeric(gplist.str[2], errors="coerce")

        dataframe.drop(columns=["FAM001_ID001", "FORMAT"], inplace=True)

        # ---- INFO: DR2 / AF / IMP ----
        info = dataframe["INFO"].astype(str)

        # DR2: extract value after DR2=, then take first comma-separated entry
        dataframe["DR2"] = pd.to_numeric(
            info.str.extract(r"(?:^|;)DR2=([^;]+)", expand=False).str.split(",", n=1).str[0],
            errors="coerce",
        )

        # AF: extract value after AF=, then take first comma-separated entry
        dataframe["AF"] = pd.to_numeric(
            info.str.extract(r"(?:^|;)AF=([^;]+)", expand=False).str.split(",", n=1).str[0],
            errors="coerce",
        )

        # IMP flag
        dataframe["IMP"] = info.str.contains(r"(?:^|;)IMP(?:;|$)", regex=True, na=False)

        dataframe.drop(columns=["INFO"], inplace=True)
        
        self.data_container.imputed_data = dataframe
    
    
    
    def load_vcf_to_df_polars(self, imputed_file) -> None:
        
        """stats from benchmarking loading and qc'ing 3 chromosomes separately:
                    Wall time: 5.48s
                    Peak RAM:  1856.5 MB
                    Avg CPU:   98.0%
                    Peak CPU:  140.8%
        """
        
        polars = (
            pl.read_csv(
                gzip.open(imputed_file, "rt", encoding="utf-8", newline=""),
                separator="\t",
                comment_prefix="#",
                has_header=False,
                new_columns=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","FAM001_ID001"],
            )

            # -------- FORMAT parsing --------
            .with_columns(
                pl.col("FAM001_ID001").str.split_exact(":", 2).alias("s"),
            )
            .with_columns(
                pl.col("s").struct.field("field_0").alias("GT"),

                # DS: keep FIRST value only
                pl.col("s")
                .struct.field("field_1")
                .str.split(",")
                .list.first()
                .cast(pl.Float64)
                .alias("DS"),

                pl.col("s").struct.field("field_2").alias("GP"),
            )
            .drop(["FAM001_ID001", "FORMAT", "s"])

            # -------- INFO parsing --------
            .with_columns(
                # DR2: extract value, keep FIRST
                pl.col("INFO")
                .str.extract(r"(?:^|;)DR2=([^;]+)", 1)
                .str.split(",")
                .list.first()
                .cast(pl.Float64)
                .alias("DR2"),

                # AF: extract value, keep FIRST
                pl.col("INFO")
                .str.extract(r"(?:^|;)AF=([^;]+)", 1)
                .str.split(",")
                .list.first()
                .cast(pl.Float64)
                .alias("AF"),

                # IMP flag
                pl.col("INFO")
                .str.contains(r"(?:^|;)IMP(?:;|$)")
                .alias("IMP"),
            )
            .drop(["INFO"])

            # -------- GP parsing --------
            .with_columns(
                pl.col("GP").str.split(",").alias("g"),
            )
            .with_columns(
                pl.col("g").list.get(0).cast(pl.Float64).alias("GP_00"),
                pl.col("g").list.get(1).cast(pl.Float64).alias("GP_01"),
                pl.col("g").list.get(2).cast(pl.Float64).alias("GP_11"),
            )
            .drop(["GP", "g"])
        ).to_pandas()
        
        self.data_container.imputed_data = polars
        
        
        
    def write_pandas_to_vcf(self, chunk_size: int = 200_000) -> None:
        """
        Stream-write imputed_data to a single 1-sample VCF file in chunks.
        Much lower peak RAM than building a full string 'body' DataFrame.
        """

        if self.data_container.qc_status is False:
            raise ValueError("Data has not passed QC; cannot write to VCF.")

        df = self.data_container.imputed_data
        if df is None or df.empty:
            raise ValueError("imputed data is empty; run qc_imputed_data() first.")

        sample_id = "imputed_sample"
        chrom_key = str(df["CHROM"].iloc[0])
        out_path = os.path.join(f"{self.environment_handler.output_dir}", f"chr{chrom_key}_imputed_qced_{self.data_container.user_id}.vcf")

        with open(out_path, "w", encoding="utf-8", newline="") as f:
            # --- header ---
            f.write("##fileformat=VCFv4.2\n")
            f.write('##INFO=<ID=DR2,Number=1,Type=Float,Description="Imputation quality (Beagle DR2)">\n')
            f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency of ALT allele">\n')
            f.write('##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed variant">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage of ALT allele">\n')
            f.write('##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype probabilities for 0/0,0/1,1/1">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")

            n = len(df)
            for start in range(0, n, chunk_size):
                chunk = df.iloc[start:start + chunk_size]

                # core fields
                chrom_col = chunk["CHROM"].astype(str)
                pos_col   = chunk["POS"].astype(int)
                id_col    = chunk["ID"].fillna(".").astype(str)
                ref_col   = chunk["REF"].astype(str)
                alt_col   = chunk["ALT"].astype(str)

                qual_raw = chunk.get("QUAL", ".")
                qual_col = qual_raw.where(qual_raw.notna(), ".").astype(str)

                filter_raw = chunk.get("FILTER", "PASS")
                filter_col = filter_raw.where(filter_raw.notna(), "PASS").astype(str)

                # INFO: build minimally
                info = pd.Series(".", index=chunk.index, dtype="object")

                if "DR2" in chunk.columns:
                    dr2 = chunk["DR2"]
                    m = dr2.notna()
                    if m.any():
                        info.loc[m] = "DR2=" + dr2.round(4).astype(str).loc[m]

                if "AF" in chunk.columns:
                    af = chunk["AF"]
                    m = af.notna()
                    if m.any():
                        add = "AF=" + af.round(6).astype(str).loc[m]
                        info.loc[m] = info.loc[m].where(info.loc[m] != ".", "")  # "." -> ""
                        sep = np.where(info.loc[m] != "", ";", "")
                        info.loc[m] = info.loc[m] + sep + add
                        info.loc[m] = info.loc[m].replace("", ".")  # restore if empty

                if "IMP" in chunk.columns:
                    m = chunk["IMP"].fillna(False).astype(bool)
                    if m.any():
                        info.loc[m] = info.loc[m].where(info.loc[m] != ".", "")  # "." -> ""
                        sep = np.where(info.loc[m] != "", ";", "")
                        info.loc[m] = info.loc[m] + sep + "IMP"
                        info.loc[m] = info.loc[m].replace("", ".")

                # FORMAT + sample
                gt_col = chunk["GT"].fillna("./.").astype(str)

                ds = chunk["DS"] if "DS" in chunk.columns else pd.Series(index=chunk.index, dtype=float)
                ds_str = ds.round(4).astype(str).where(ds.notna(), ".")

                gp_str = pd.Series(".", index=chunk.index, dtype="object")
                if all(k in chunk.columns for k in ("GP_00", "GP_01", "GP_11")):
                    gp00, gp01, gp11 = chunk["GP_00"], chunk["GP_01"], chunk["GP_11"]
                    m = gp00.notna() & gp01.notna() & gp11.notna()
                    if m.any():
                        gp_str.loc[m] = (
                            gp00.round(4).astype(str).loc[m] + "," +
                            gp01.round(4).astype(str).loc[m] + "," +
                            gp11.round(4).astype(str).loc[m]
                        )

                format_col = "GT:DS:GP"
                sample_col = gt_col + ":" + ds_str + ":" + gp_str

                # write chunk (no huge full-body dataframe stored beyond chunk)
                out_df = pd.DataFrame({
                    "#CHROM": chrom_col,
                    "POS": pos_col,
                    "ID": id_col,
                    "REF": ref_col,
                    "ALT": alt_col,
                    "QUAL": qual_col,
                    "FILTER": filter_col,
                    "INFO": info,
                    "FORMAT": format_col,
                    sample_id: sample_col,
                })

                out_df.to_csv(f, sep="\t", index=False, header=False)

        self.environment_handler.qc_imputed_files[chrom_key] = out_path
        print(f"Wrote QC'ed imputed data to {out_path}")

        
        


    
    def run_qc_on_imputed_data(self, engine) -> None:
        
        #Per chromsome approach to save memory
        
        for chromosome, file in self.environment_handler.imputed_files.items():
            print(f"Converting chromosome {chromosome} to pandas DataFrame...")
            
            if engine == "pandas":
                self.load_vcf_to_df_pandas(file)
            else:
                self.load_vcf_to_df_polars(file)
            
            print(f"Running QC on chromosome {chromosome} imputed data...")
            self.data_container.qc_imputed_data()
            print(f"Writing QC'ed imputed data for chromosome {chromosome} to VCF...")
            self.write_pandas_to_vcf()
            
            # Force garbage collection to free memory, otherwise it will ramp up
            gc.collect()
            print(f"Completed chromosome {chromosome}. Memory freed.")

    
    def create_vcf_reference_mapping(self) -> None:
        
        mapping_df = pd.DataFrame()
        
        
        vcf_files_df=pd.DataFrame(
            list(self.environment_handler.vcf_file_paths.items()),
            columns=["chromosome_number", "vcf_file"]
        )
        
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
    

        
        self.environment_handler.vcf_plink_reference_mapping = mapping_df
        
        
    def set_user_id_from_vcf(self) -> None:
        vcf_file = self.environment_handler.vcf_file_paths["22"]  # smallest chr
        opener = gzip.open if vcf_file.endswith(".gz") else open

        with opener(vcf_file, "rt") as f:  # <-- text mode
            for line in f:
                if line.startswith("#CHROM"):
                    sample_id = line.rstrip("\n").split("\t")[-1]
                    self.data_container.user_id = sample_id.split("_", 1)[1]
                    break
        
        
    
    def create_samplesheet(self) -> None:
        
        vcf_paths= list(self.environment_handler.qc_imputed_files.values())
        
        vcf_filenames = [os.path.basename(path) for path in vcf_paths]
        
        vcf_base_names = [name.split(".", 1)[0] for name in vcf_filenames]
        
        
        chr_numbers = list(self.environment_handler.qc_imputed_files.keys())

        samplesheet_df = pd.DataFrame({
            "sampleset": self.data_container.user_id,
            "path_prefix": vcf_base_names,
            "chrom": chr_numbers,
            "format": "vcf",
        })
        
        samplesheet_df.to_csv(
            os.path.join(self.environment_handler.output_dir, "samplesheet.csv"),
            index=False
        )
        
    


    

