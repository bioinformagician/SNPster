import os

BEAGLE_JAR = os.getenv("BEAGLE_JAR", "/app/dependencies/beagle.27Feb25.75f.jar")
JAVA_EXE = os.getenv("JAVA_EXE", "/usr/bin/java")
VCF_REFERENCE_MAPPING = os.getenv("VCF_REFERENCE_MAPPING", "/work/harmonization_results/vcf_reference_mapping.parquet")
OUTPUT_DIR = os.getenv("OUTPUT_DIR", "/output")

# use docker environment variables or default values
HEAP_GB = int(os.getenv("HEAP_GB", "8"))
THREADS = int(os.getenv("THREADS", "2"))
GP_MIN = float(os.getenv("GP_MIN", "0.90"))
DS_TOL = float(os.getenv("DS_TOL", "0.05"))
SNPS_ONLY = os.getenv("SNPS_ONLY", "True").lower() == "true"
BIALLELIC_ONLY = os.getenv("BIALLELIC_ONLY", "True").lower() == "true"