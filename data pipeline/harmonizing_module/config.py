import os

#use docker environment variables or default values

PLINK_PREFIX = os.getenv("PLINK_PREFIX", "plink_temp")


PLINK_REFERENCE_FASTA = os.getenv(
    "PLINK_REFERENCE_FASTA",
    r"C:\Users\frezz\Desktop\harmonizer_dependencies\Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
)

TEST_FILE = os.getenv(
    "TEST_FILE",
    r"C:\Users\frezz\Downloads\snpster\data pipeline\harmonizing_module\dependencies\standardized_microarray_data.parquet",
)
PVAR_REF_FILE = os.getenv(
    "PVAR_REF_FILE",
    r"C:\Users\frezz\Desktop\harmonizer_dependencies\all_phase3.pvar.zst",
)

PLINK_1_9_PATH = os.getenv(
    "PLINK_1_9_PATH",
    r"C:\Users\frezz\Desktop\harmonizer_dependencies\plink.exe",
)

PLINK_1_9_MEMORY_MB = os.getenv(
    "PLINK_1_9_MEMORY_MB",
    "1000",
)

PLINK_1_9_THREADS = os.getenv(
    "PLINK_1_9_THREADS",
    "1",
)


PLINK_2_0_PATH = os.getenv(
    "PLINK_2_0_PATH",
    r"C:\Users\frezz\Desktop\harmonizer_dependencies\plink2.exe",
)

