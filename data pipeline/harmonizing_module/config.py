import os

#use docker environment variables or default values

PLINK_PREFIX = os.getenv("PLINK_PREFIX", "plink_temp")


PLINK_REFERENCE_FASTA = os.getenv(
    "PLINK_REFERENCE_FASTA",
    "/data/hs37d5.fa.zst",
)

TEST_FILE = os.getenv(
    "TEST_FILE",
    "/app/dependencies/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt",
)
PVAR_REF_FILE = os.getenv(
    "PVAR_REF_FILE",
    "/data/all_phase3.pvar.zst",
)
PLINK_1_9_PATH = os.getenv(
    "PLINK_1_9_PATH",
    "/app/dependencies/plink",
)

PLINK_1_9_MEMORY_MB = os.getenv(
    "PLINK_1_9_MEMORY_MB",
    "4000",
)

PLINK_2_0_PATH = os.getenv(
    "PLINK_2_0_PATH",
    "/app/dependencies/plink2",
)
