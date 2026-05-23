
import os

ENGINE = os.getenv("ENGINE", "polars").lower()
THREADS = int(os.getenv("THREADS", "2"))
GP_MIN = float(os.getenv("GP_MIN", "0.90"))
DS_TOL = float(os.getenv("DS_TOL", "0.05"))
SNPS_ONLY = os.getenv("SNPS_ONLY", "True").lower() == "true"
BIALLELIC_ONLY = os.getenv("BIALLELIC_ONLY", "True").lower() == "true"
ACCEPTED_DF_ENGINES = ["pandas", "polars"]