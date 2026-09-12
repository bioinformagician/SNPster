import os

OUTPUT_DIR = os.getenv("OUTPUT_DIR")
REFERENCE_DATA_PATH = os.getenv("REFERENCE_DATA_PATH")
SCORING_FILE_SOURCE_DIR = os.getenv("SCORING_FILE_SOURCE_DIR")
SCORING_FILE_TARGET_DIR = os.getenv("SCORING_FILE_TARGET_DIR")
NF_WORK_DIR = os.getenv("NF_WORK_DIR", "/app/work")
BCFTOOLS_THREADS = int(os.getenv("BCFTOOLS_THREADS", "8"))
SAMPLESET_NAME = os.getenv("SAMPLESET_NAME", "MergedUserData")
VCF_MERGE_SHEET_DIR = os.getenv("VCF_MERGE_SHEET_DIR", "/app/vcf_merge_sheets")
N_PRSC_JOBS = int(os.getenv("N_PRSC_JOBS", "303"))
PGS_RESULT_DIR = os.getenv("PGS_RESULT_DIR", "/app/pgs_results")
NEXTFLOW_PGS_CONFIG = os.getenv("NEXTFLOW_PGS_CONFIG", "/app/nextflow_pgs.config")
NEXTFLOW_PGS_RESOURCE_CONFIG = os.getenv("NEXTFLOW_PGS_RESOURCE_CONFIG", "/app/nextflow_pgs_resources.config")
NEXTFLOW_VCF_MERGING_CONFIG = os.getenv("NEXTFLOW_VCF_MERGING_CONFIG", "/app/nextflow_vcf_merging.config")

# Resource envelope the pgsc_calc executor is allowed to use; chunk sizes are derived from it.
PGS_TOTAL_CPUS = int(os.getenv("PGS_TOTAL_CPUS", "11"))
PGS_TOTAL_MEM_GB = float(os.getenv("PGS_TOTAL_MEM_GB", "110"))
PGS_MEM_SAFETY_FRACTION = float(os.getenv("PGS_MEM_SAFETY_FRACTION", "0.8"))
PGS_MATCH_VARIANTS_CPUS = int(os.getenv("PGS_MATCH_VARIANTS_CPUS", "2"))
# MATCH_COMBINE runs alone and passes task.cpus straight to POLARS_MAX_THREADS; gains flatten off
# past ~6 threads because a large part of the step is single-threaded.
PGS_MATCH_COMBINE_CPUS = int(os.getenv("PGS_MATCH_COMBINE_CPUS", "6"))
# COMBINE_SCOREFILES sets no thread limit and still only reaches ~110% CPU, so extra cores are wasted.
PGS_COMBINE_SCOREFILES_CPUS = int(os.getenv("PGS_COMBINE_SCOREFILES_CPUS", "1"))
PGS_TARGET_VARIANTS_PER_CHUNK = int(os.getenv("PGS_TARGET_VARIANTS_PER_CHUNK", "85000000"))
PGS_MAX_FILES_PER_CHUNK = int(os.getenv("PGS_MAX_FILES_PER_CHUNK", "1000"))

# Reference panel samples that pgsc_calc carries in every run, on top of the user cohort.
PGS_REFERENCE_PANEL_SAMPLES = int(os.getenv("PGS_REFERENCE_PANEL_SAMPLES", "4000"))
# Multiplier applied to the fitted peak RSS before it is requested from the executor.
PGS_MEM_HEADROOM_FACTOR = float(os.getenv("PGS_MEM_HEADROOM_FACTOR", "1.10"))
# One MATCH_VARIANTS task as (serial_seconds, parallelisable_cpu_seconds) at the variant count
# below, fitted from 2m20s at 187% CPU on 2 cpus in the nf_reports traces.
PGS_MATCH_VARIANTS_TASK_MODEL_S = (18.0, 244.0)
PGS_MODEL_REFERENCE_VARIANTS = 85_000_000