import os

#use docker environment variables or default values

PLINK_PREFIX = os.getenv("PLINK_PREFIX", "plink_temp")

PLINK_MAP_DIR = os.getenv(
    "PLINK_MAP_DIR",
    "/data/plink.GRCh37.map",
)
PLINK_REFERENCE_FASTA = os.getenv(
    "PLINK_REFERENCE_FASTA",
    "/data/hs37d5.fa.zst",
)
BEAGLE_REFERENCE_DIR = os.getenv(
    "BEAGLE_REFERENCE_DIR",
    "/data/beagle_references",
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
PLINK_2_0_PATH = os.getenv(
    "PLINK_2_0_PATH",
    "/app/dependencies/plink2",
)
