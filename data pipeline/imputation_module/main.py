import os
import argparse
from imputation_classes import EnvironmentHandler, QCThresholds, DataContainer, WorkflowOrchestrator, vcf_classes
from config import BEAGLE_JAR, JAVA_EXE, HEAP_GB, THREADS, GP_MIN, DS_TOL, SNPS_ONLY, BIALLELIC_ONLY, OUTPUT_DIR, BEAGLE_REFERENCE_DIR, PLINK_MAP_DIR, VCF_FILES_DIR, DF_ENGINE, ACCEPTED_DF_ENGINES, MERGED_SAMPLESHEET_DIR

parser = argparse.ArgumentParser()
parser.add_argument('--vcf_files', type=str, required=False, default=VCF_FILES_DIR) #default is vcf files for testing baked into image
parser.add_argument('--df_engine', type=str, required=False, default=DF_ENGINE)
args = parser.parse_args()

if args.df_engine not in ACCEPTED_DF_ENGINES:
    raise ValueError(f"df_engine must be one of {ACCEPTED_DF_ENGINES}")

environment_handler = EnvironmentHandler(
    working_dir=os.getcwd(),
    java_exe=JAVA_EXE,
    beagle_jar=BEAGLE_JAR,
    heap_gb=HEAP_GB,
    threads=THREADS,
    output_dir=OUTPUT_DIR,
    merged_samplesheet_dir = MERGED_SAMPLESHEET_DIR,
    beagle_reference_dir=BEAGLE_REFERENCE_DIR, #mounted
    plink_map_dir=PLINK_MAP_DIR, #mounted
    vcf_files_dir=args.vcf_files, #mounted from harmonization step
)

qc_thresholds = QCThresholds(
    gp_min=GP_MIN,
    ds_tol=DS_TOL,
    snps_only=SNPS_ONLY,
    biallelic_only=BIALLELIC_ONLY
)

vcf_environment_handler = vcf_classes.VCFEnvironmentHandler(
    vcf_samplesheet_path = None,
    output_dir = None)

vcf_handler = vcf_classes.VCFHandler(
    vcf_environment_handler=vcf_environment_handler
)

vcf_handler.vcf_environment_handler.output_dir = environment_handler.merged_samplesheet_dir

orchestrator = WorkflowOrchestrator(
    environment_handler=environment_handler,
    data_containers=[],
    vcf_handler=vcf_handler,
    qc_thresholds=qc_thresholds
)


orchestrator.create_vcf_samplesheet_for_vcf_merge()

orchestrator.split_vcf_samplesheet_by_chromosome()

vcf_handler.run_vcf_merging(environment_handler.merged_samplesheet_dir) #merging in parallel

orchestrator.update_vcf_plink_reference_mapping()

orchestrator.impute_vcf_files()# parallelize?

orchestrator.create_vcf_samplesheets_for_vcf_split()

vcf_handler.run_vcf_splitting(environment_handler.merged_samplesheet_dir)

orchestrator.delete_merged_vcf()

orchestrator.setup_datacontainers()

orchestrator.run_qc_on_imputed_data(args.df_engine) #also need to do this in parallel with nextflow, rly slow step

orchestrator.create_samplesheet()

print("Imputation and QC pipeline completed successfully.")
