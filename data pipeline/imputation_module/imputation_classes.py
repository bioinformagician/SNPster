import os
import pandas as pd
from dataclasses import dataclass
import gzip
import subprocess
import numpy as np


@dataclass(frozen=True)
class QCThresholds:
    gp_min: float = 0.90
    ds_tol: float = 0.10
    snps_only: bool = True
    biallelic_only: bool = True



@dataclass
class DataContainer:
    qc_thresholds: QCThresholds = QCThresholds()
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

        out = qc[["rsid", "ref", "alt", "GT"]].rename(columns={"GT": "gt"})
        self.qced_imputed_data = out
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
                 vcf_plink_reference_mapping: pd.DataFrame,
                 imputed_dir: str = None,
                 imputed_files: dict[str,str] = None
                 ):
        
        self.working_dir = working_dir
        self.java_exe = java_exe
        self.beagle_jar = beagle_jar
        self.vcf_plink_reference_mapping = vcf_plink_reference_mapping
        self.imputed_dir = imputed_dir
        self.heap_gb = heap_gb
        self.threads = threads
        self.imputed_files = imputed_files
        self.validate_paths()
        self.read_parquet()
        self.make_imputed_directory()
    

    def read_parquet(self) -> pd.DataFrame:
        self.vcf_plink_reference_mapping = pd.read_parquet(self.vcf_plink_reference_mapping)
    
    def make_imputed_directory(self) -> None:
        imputed_dir = os.path.join(self.working_dir, "imputed")
        os.makedirs(imputed_dir, exist_ok=True)
        self.working_dir = imputed_dir

    def validate_paths(self) -> None:
        
        for path in [
            self.working_dir,
            self.java_exe,
            self.beagle_jar,
            self.vcf_plink_reference_mapping
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
        

        out = os.path.join(f"{self.environment_handler.working_dir}/imputed_chr{chr_number}.risk")
        
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
    


if __name__ == "__main__":
    from config import BEAGLE_JAR, JAVA_EXE, HEAP_GB, THREADS, GP_MIN, DS_TOL, SNPS_ONLY, BIALLELIC_ONLY
    
    environment_handler = EnvironmentHandler(
        working_dir=r"C:\Users\frezz\pipeline_testing\temp_working_dir_7213",
        java_exe=JAVA_EXE,
        beagle_jar=BEAGLE_JAR,
        heap_gb=HEAP_GB,
        threads=THREADS,
        vcf_plink_reference_mapping = r"C:\Users\frezz\pipeline_testing\temp_working_dir_7213\results\vcf_reference_mapping.parquet",
        imputed_files = {}
    )
    
    qc_thresholds = QCThresholds(
        gp_min=GP_MIN,
        ds_tol=DS_TOL,
        snps_only=SNPS_ONLY,
        biallelic_only=BIALLELIC_ONLY
    )
    
    data_container = DataContainer(
        qc_thresholds=qc_thresholds
    )
    
    orchestrator = WorkflowOrchestrator(
        environment_handler=environment_handler,
        data_container=data_container
    )
    
    orchestrator.impute_vcf_files()
    
    orchestrator.load_vcf_to_df()
    
    orchestrator.data_container.qc_imputed_data()
    
    print(orchestrator.data_container.qced_imputed_data.head())
    
    

