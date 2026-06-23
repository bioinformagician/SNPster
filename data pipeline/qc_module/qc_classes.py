from dataclasses import dataclass
import os
import pandas as pd
import numpy as np
import polars as pl
import subprocess
import gzip
import re


@dataclass(frozen=True)
class QCThresholds:
    gp_min: float
    ds_tol: float
    snps_only: bool
    biallelic_only: bool



@dataclass
class ImputedDataContainer:
    file_path: str = None # path to the split imputed file
    qc_thresholds: QCThresholds = QCThresholds(None, None, False, False)
    imputed_data: pd.DataFrame = None
    qc_status: bool = False
    imputation_id: int = None
    chromosome: str = None
    
    def set_chromosome_from_data(self):
        self.chromosome = str(self.imputed_data["CHROM"].iloc[0])
    
    
    def load_vcf_to_df_pandas(self) -> None:
        
        """stats from benchmarking loading and qc'ing 3 chromosomes separately:
                    Wall time: 24.98s
                    Peak RAM:  1332.8 MB
                    Avg CPU:   99.6%
                    Peak CPU:  141.1%
        """
        
        # Use the Python CSV engine here to avoid rare native-parser crashes on
        # malformed user-derived rows while keeping QC deterministic.
        dataframe = pd.read_csv(
            self.file_path,
            sep="\t",
            comment="#",
            compression="gzip",
            header=None,
            dtype=str,
            engine="python",
            on_bad_lines="skip",
        )

        if dataframe.empty or dataframe.shape[1] < 10:
            raise ValueError(f"VCF appears malformed or empty after parsing: {self.file_path}")

        # Keep the expected VCF columns and a single-sample field.
        dataframe = dataframe.iloc[:, :10].copy()
        dataframe.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

        # ---- FORMAT: GT / DS / GP ----
        fmt = dataframe["SAMPLE"].str.split(":", n=2, expand=True)
        dataframe["GT"] = fmt[0]

        # DS: keep first value only
        dataframe["DS"] = pd.to_numeric(fmt[1].str.split(",", n=1).str[0], errors="coerce")

        gplist = fmt[2].astype(str).str.split(",")
        dataframe["GP_00"] = pd.to_numeric(gplist.str[0], errors="coerce")
        dataframe["GP_01"] = pd.to_numeric(gplist.str[1], errors="coerce")
        dataframe["GP_11"] = pd.to_numeric(gplist.str[2], errors="coerce")

        dataframe.drop(columns=["SAMPLE", "FORMAT"], inplace=True)

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
        
        self.imputed_data = dataframe
    
    
    
    def load_vcf_to_df_polars(self) -> None:
        
        """stats from benchmarking loading and qc'ing 3 chromosomes separately:
                    Wall time: 5.48s
                    Peak RAM:  1856.5 MB
                    Avg CPU:   98.0%
                    Peak CPU:  140.8%
        """
        
        polars = (
            pl.read_csv(
                gzip.open(self.file_path, "rt", encoding="utf-8", newline=""),
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
        
        self.imputed_data = polars
    
    

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
    
    
    
    
    
    def write_pandas_to_vcf(self, output_dir, chunk_size: int = 500000) -> None:
        """
        Stream-write imputed_data to a single 1-sample VCF file in chunks.
        Much lower peak RAM than building a full string 'body' DataFrame.
        """

        if self.qc_status is False:
            raise ValueError("Data has not passed QC; cannot write to VCF.")

        df = self.imputed_data
        if df is None or df.empty:
            raise ValueError("imputed data is empty; run qc_imputed_data() first.")

        #add qc to the filename
        out_path = os.path.join(f"{output_dir}", os.path.basename(self.file_path).replace(".vcf.gz", "_QC.vcf"))
        
        with open(out_path, "w", encoding="utf-8", newline="") as f:
            # --- header ---
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##contig=<ID={self.chromosome}>\n")
            f.write('##INFO=<ID=DR2,Number=1,Type=Float,Description="Imputation quality (Beagle DR2)">\n')
            f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency of ALT allele">\n')
            f.write('##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed variant">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage of ALT allele">\n')
            f.write('##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype probabilities for 0/0,0/1,1/1">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + f"IMPID_{str(self.imputation_id)}" + "\n")

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
                    f"IMPID_{str(self.imputation_id)}": sample_col,
                })

                out_df.to_csv(f, sep="\t", index=False, header=False)
        

        print(f"Wrote QC'ed imputed data to {out_path}")
        
        return out_path
    
    def zip_vcf(self, vcf_path) -> str:
        """Compress VCF using bgzip for block compression (allows tabix indexing)"""
        gz_path = vcf_path + ".gz"
        subprocess.run(["bgzip", "-c", vcf_path], stdout=open(gz_path, "wb"), check=True)
        os.remove(vcf_path)
        return gz_path
    
    def run_qc_pipeline(self, input_dir, output_dir) -> None:
        
        
        command = [
            "nextflow", "run", "vcf_combiner_pipeline.nf",
            "--vcf_file_dir", input_dir,
            "--output_dir", output_dir
        ]
        
        subprocess.run(command, check=True)