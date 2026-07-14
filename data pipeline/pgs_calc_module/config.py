import os

OUTPUT_DIR = os.getenv("OUTPUT_DIR")
REFERENCE_DATA_PATH = os.getenv("REFERENCE_DATA_PATH")
SCORING_FILE_SOURCE_DIR = os.getenv("SCORING_FILE_SOURCE_DIR")
SCORING_FILE_TARGET_DIR = os.getenv("SCORING_FILE_TARGET_DIR")
NF_WORK_DIR = os.getenv("NF_WORK_DIR", "/app/work")
BCFTOOLS_THREADS = int(os.getenv("BCFTOOLS_THREADS", "8"))
SAMPLESET_NAME = os.getenv("SAMPLESET_NAME", "MergedUserData")
VCF_MERGE_SHEET_DIR = os.getenv("VCF_MERGE_SHEET_DIR", "/app/vcf_merge_sheets")
PGS_CHUNK_SIZE = int(os.getenv("PGS_CHUNK_SIZE", "50"))
PGS_RESULT_DIR = os.getenv("PGS_RESULT_DIR", "/app/pgs_results")