import os
import pandas as pd
from dataclasses import dataclass
import gzip
import subprocess
import numpy as np
import re

@dataclass(frozen=True)
class QCThresholds:
    gp_min: float
    ds_tol: float
    snps_only: bool
    biallelic_only: bool



@dataclass
class DataContainer:
    qc_thresholds: QCThresholds = QCThresholds(None, None, False, False)
    imputed_data: pd.DataFrame = None,
    qced_imputed_data: pd.DataFrame = None
    
    

    def qc_imputed_data(self) -> None:
        
        print("starting qc of imputed data...")
        start_time = pd.Timestamp.now()
        
        x = self.imputed_data
        n = len(x)

        if n == 0:
            print("0.0% of input variants passed qc requirements (no variants)")
            self.qced_imputed_data = x.iloc[0:0][["rsid", "ref", "alt", "GT"]].rename(columns={"GT": "gt"})
            return

        # ---- build a single mask instead of chaining slices ----
        keep_mask = np.ones(n, dtype=bool)

        # SNP-only: ref and alt length == 1
        if self.qc_thresholds.snps_only:
            ref_len = x["ref"].str.len().to_numpy()
            alt_len = x["alt"].str.len().to_numpy()
            keep_mask &= (ref_len == 1) & (alt_len == 1)

        # biallelic-only: no ',' in alt
        if self.qc_thresholds.biallelic_only:
            has_comma = x["alt"].str.contains(",").to_numpy()
            keep_mask &= ~has_comma

        # apply SNP/biallelic filters once
        x = x.loc[keep_mask]
        n = len(x)
        if n == 0:
            print("0.0% of input variants passed qc requirements (after SNP/biallelic filter)")
            self.qced_imputed_data = self.imputed_data.iloc[0:0][["rsid", "ref", "alt", "GT"]].rename(columns={"GT": "gt"})
            return

        # ---- GP certainty (already numeric) ----
        gp = x[["GP_00", "GP_01", "GP_11"]].to_numpy()
        gpmax = gp.max(axis=1)

        # ---- genotype classes ----
        gt = x["GT"].astype(str).to_numpy()

        is_hom_ref = (gt == "0|0") | (gt == "0/0")
        is_het     = (gt == "0|1") | (gt == "1|0") | (gt == "0/1") | (gt == "1/0")
        is_hom_alt = (gt == "1|1") | (gt == "1/1")

        ds = x["DS"].to_numpy(dtype=float)
        gp01 = x["GP_01"].to_numpy(dtype=float)

        hom_ref_ok = is_hom_ref & (ds <= self.qc_thresholds.ds_tol)
        het_ok     = is_het & (gp01 >= self.qc_thresholds.gp_min) & (ds >= 1 - self.qc_thresholds.ds_tol) & (ds <= 1 + self.qc_thresholds.ds_tol)
        hom_alt_ok = is_hom_alt & (ds >= 2 - self.qc_thresholds.ds_tol)

        keep = (gpmax >= self.qc_thresholds.gp_min) & (hom_ref_ok | het_ok | hom_alt_ok)

        qc = x.loc[keep]

        kept = len(qc)
        total = len(x)
        pct = (kept / total * 100.0) if total else 0.0
        print(f"{pct:.1f}% of input variants passed qc requirements")

        # keep full QC'ed dataframe (all columns)
        self.qced_imputed_data = qc

        end_time = pd.Timestamp.now()
        duration = end_time - start_time
        print(f"QC of imputed data completed in {duration}")
    
    
    

    
    

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
        self.read_parquet()
        self.make_imputed_directory()
        self.set_beagle_files()
        self.set_plink_map_files()
        self.set_vcf_files()
        
        
    def set_beagle_files(self) -> dict:
        beagle_files = {}
        for file in os.listdir(self.beagle_reference_dir):
            if file.endswith(".bref3"):
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
        for file in os.listdir(self.vcf_files_dir):
            if file.endswith(".vcf.gz"):
                chrom = re.search(r"chr(\d+)(?=\.vcf(?:\.gz)?$)", file)
                if chrom:
                    vcf_files[chrom.group(1)] = os.path.join(self.vcf_files_dir, file)
        
        self.vcf_file_paths = vcf_files
    

    def read_parquet(self) -> pd.DataFrame:
        self.vcf_plink_reference_mapping = pd.read_parquet(self.vcf_plink_reference_mapping)
    
    def make_imputed_directory(self) -> None:
        imputed_dir = os.path.join(self.working_dir, "imputed")
        os.makedirs(imputed_dir, exist_ok=True)
        self.imputed_dir = imputed_dir
    

    def validate_paths(self) -> None:
        
        for path in [
            self.working_dir,
            self.java_exe,
            self.beagle_jar,
            self.vcf_plink_reference_mapping,
            self.output_dir,
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
            
            
    
    
    def _parse_info(self, info_str: str) -> dict:
            # INFO like: "DR2=0.98;AF=0.12;IMP"
            d = {}
            for kv in info_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = v
                else:
                    # flag field (e.g. IMP)
                    d[kv] = True
            return d

    def _parse_format(self, fmt: str, sample: str) -> dict:
        # e.g. fmt: "GT:DS:GP", sample: "0/1:1.02:0.01,0.95,0.04"
        keys = fmt.split(":")
        vals = sample.split(":")
        return dict(zip(keys, vals))

    def load_vcf_to_df(self) -> None:
        keep_info = ("DR2", "AF")
        keep_format = ("GT", "DS", "GP")

        all_rows = []

        for imputed_file in self.environment_handler.imputed_files.values():

            with gzip.open(imputed_file, "rt", encoding="utf-8", newline="") as f:
                for line in f:
                    # skip all header lines
                    if line.startswith("#"):
                        continue

                    toks = line.rstrip("\n").split("\t")
                    chrom, pos, rsid, ref, alt, qual, flt, info, fmt, sample_str = toks[:10]

                    info_map = self._parse_info(info)
                    fmt_map = self._parse_format(fmt, sample_str)

                    rec = {
                        "chrom": chrom,
                        "pos": int(pos),
                        "rsid": rsid,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual,
                        "filter": flt,
                    }

                    # selected INFO fields (assumed always present + numeric)
                    for k in keep_info:
                        rec[k] = float(info_map[k].split(",")[0])

                    # selected FORMAT fields (assumed always present)
                    for k in keep_format:
                        rec[k] = fmt_map[k]

                    # split GP into three floats (assumed always valid)
                    gp_vals = [float(v) for v in fmt_map["GP"].split(",")[:3]]
                    gp0, gp1, gp2 = gp_vals
                    rec["GP_00"] = gp0
                    rec["GP_01"] = gp1
                    rec["GP_11"] = gp2 
                    rec["DS"] = float(rec["DS"].split(",")[0])

                    # IMP flag: present in INFO => True, else False
                    rec["IMP"] = bool(info_map.get("IMP", False))

                    all_rows.append(rec)

        dataframe = pd.DataFrame(all_rows)
        self.data_container.imputed_data = dataframe
    
    
    
    def write_pandas_to_vcf(self) -> None:
        """
        Write qced_imputed_data to one VCF file per chromosome (1-sample).
        Produces files like: {working_dir}/qc_imputed/imputed_chr{chrom}.vcf

        Expects columns:
        chrom, pos, rsid, ref, alt, qual, filter,
        DR2, AF, GT, DS, GP_00, GP_01, GP_11, IMP
        """

        df = self.data_container.qced_imputed_data
        if df is None or df.empty:
            raise ValueError("qced_imputed_data is empty; run qc_imputed_data() first.")

        sample_id = "imputed_sample"

        # group by chromosome
        for chrom, df_chr in df.groupby("chrom"):
            print(f"Writing QC'ed imputed data for chromosome {chrom}...")

            # sort within chromosome
            df_chr = df_chr.sort_values("pos")

            # ---- build columns needed for VCF body (vectorized) ----

            # CHROM, POS, ID, REF, ALT
            chrom_col = df_chr["chrom"].astype(str)
            pos_col   = df_chr["pos"].astype(int)
            id_col    = df_chr["rsid"].fillna(".").astype(str)
            ref_col   = df_chr["ref"].astype(str)
            alt_col   = df_chr["alt"].astype(str)

            # QUAL & FILTER
            qual_raw = df_chr.get("qual", ".")
            qual_col = qual_raw.where(~pd.isna(qual_raw), ".").astype(str)

            filter_raw = df_chr.get("filter", "PASS")
            filter_col = filter_raw.where(~pd.isna(filter_raw), "PASS").astype(str)

            # ---- INFO field: DR2, AF, IMP ----
            info = pd.Series("", index=df_chr.index, dtype="object")

            # DR2
            if "DR2" in df_chr.columns:
                dr2 = df_chr["DR2"]
                mask_dr2 = dr2.notna()
                if mask_dr2.any():
                    dr2_str = dr2.round(4).astype(str)
                    info[mask_dr2] = "DR2=" + dr2_str[mask_dr2]

            # AF
            if "AF" in df_chr.columns:
                af = df_chr["AF"]
                mask_af = af.notna()
                if mask_af.any():
                    af_str = af.round(6).astype(str)
                    sep = np.where(info[mask_af] != "", ";", "")
                    info.loc[mask_af] = info[mask_af] + sep + "AF=" + af_str[mask_af]

            # IMP flag
            if "IMP" in df_chr.columns:
                mask_imp = df_chr["IMP"].fillna(False).astype(bool)
                if mask_imp.any():
                    sep = np.where(info[mask_imp] != "", ";", "")
                    info.loc[mask_imp] = info[mask_imp] + sep + "IMP"

            info_col = info.replace("", ".")

            # ---- FORMAT + sample field ----
            # GT
            gt_col = df_chr["GT"].fillna("./.").astype(str)

            # DS
            ds = df_chr["DS"] if "DS" in df_chr.columns else pd.Series(index=df_chr.index, dtype=float)
            ds_str = ds.round(4).astype(str)
            ds_str = ds_str.where(~ds.isna(), ".")

            # GP
            gp00 = df_chr.get("GP_00")
            gp01 = df_chr.get("GP_01")
            gp11 = df_chr.get("GP_11")

            # default all GP to "."
            gp_str = pd.Series(".", index=df_chr.index, dtype="object")
            if gp00 is not None and gp01 is not None and gp11 is not None:
                mask_gp = gp00.notna() & gp01.notna() & gp11.notna()
                if mask_gp.any():
                    gp00_str = gp00.round(4).astype(str)
                    gp01_str = gp01.round(4).astype(str)
                    gp11_str = gp11.round(4).astype(str)
                    gp_str.loc[mask_gp] = (
                        gp00_str[mask_gp] + "," +
                        gp01_str[mask_gp] + "," +
                        gp11_str[mask_gp]
                    )

            format_col = pd.Series("GT:DS:GP", index=df_chr.index)
            sample_col = gt_col + ":" + ds_str + ":" + gp_str

            # assemble final VCF body dataframe
            body = pd.DataFrame({
                "#CHROM": chrom_col,
                "POS": pos_col,
                "ID": id_col,
                "REF": ref_col,
                "ALT": alt_col,
                "QUAL": qual_col,
                "FILTER": filter_col,
                "INFO": info_col,
                "FORMAT": format_col,
                sample_id: sample_col,
            })

            out_path = os.path.join(self.environment_handler.output_dir, f"imputed_chr{chrom}.vcf")
            with open(out_path, "w", encoding="utf-8", newline="") as f:
                # --- minimal header ---
                f.write("##fileformat=VCFv4.2\n")
                f.write('##INFO=<ID=DR2,Number=1,Type=Float,Description="Imputation quality (Beagle DR2)">\n')
                f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency of ALT allele">\n')
                f.write('##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed variant">\n')
                f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                f.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage of ALT allele">\n')
                f.write('##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype probabilities for 0/0,0/1,1/1">\n')
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")

                # write body in one go (no Python row loop)
                body.to_csv(f, sep="\t", index=False, header=False)

            self.environment_handler.qc_imputed_files[chrom] = out_path
            print(f"Wrote QC'ed imputed data for chromosome {chrom} to {out_path}")
    
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
        

        
        self.environment_handler.vcf_plink_reference_mapping = mapping_df

    


    

