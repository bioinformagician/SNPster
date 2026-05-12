import os

BEAGLE_JAR = os.getenv("BEAGLE_JAR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\dependencies\beagle.27Feb25.75f.jar")
JAVA_EXE = os.getenv("JAVA_EXE", r"C:\Program Files\Java\jre1.8.0_471\bin\java.exe")
OUTPUT_DIR = os.getenv("OUTPUT_DIR", r"C:\Users\frezz\Desktop\imputer_output\testing")
VCF_FILES_DIR = os.getenv("VCF_FILES_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\dependencies\test_data")
PLINK_MAP_DIR = os.getenv("PLINK_MAP_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\big_dependencies\plink.GRCh38.map")
DF_ENGINE = os.getenv("DF_ENGINE", "polars")

BEAGLE_REFERENCE_DIR = os.getenv("BEAGLE_REFERENCE_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\big_dependencies\beagle_references")
# use docker environment variables or default values
HEAP_GB = int(os.getenv("HEAP_GB", "8"))
THREADS = int(os.getenv("THREADS", "2"))
GP_MIN = float(os.getenv("GP_MIN", "0.90"))
DS_TOL = float(os.getenv("DS_TOL", "0.05"))
SNPS_ONLY = os.getenv("SNPS_ONLY", "True").lower() == "true"
BIALLELIC_ONLY = os.getenv("BIALLELIC_ONLY", "True").lower() == "true"
ACCEPTED_DF_ENGINES = ["pandas", "polars"]
MERGED_SAMPLESHEET_DIR = os.getenv("MERGED_SAMPLESHEET_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\merged_samplesheets")