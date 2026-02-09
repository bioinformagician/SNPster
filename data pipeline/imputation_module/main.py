import os
import argparse
from imputation_classes import EnvironmentHandler, QCThresholds, DataContainer, WorkflowOrchestrator
from config import BEAGLE_JAR, JAVA_EXE, HEAP_GB, THREADS, GP_MIN, DS_TOL, SNPS_ONLY, BIALLELIC_ONLY, OUTPUT_DIR, BEAGLE_REFERENCE_DIR, PLINK_MAP_DIR, VCF_FILES_DIR, DF_ENGINE, ACCEPTED_DF_ENGINES

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
    beagle_reference_dir=BEAGLE_REFERENCE_DIR, #mounted
    plink_map_dir=PLINK_MAP_DIR, #mounted
    vcf_files_dir=args.vcf_files, #mounted from harmonization step
    imputed_files = {},
    qc_imputed_files= {}
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


orchestrator.set_user_id_from_vcf()

orchestrator.make_result_subdir()

orchestrator.impute_vcf_files()

orchestrator.run_qc_on_imputed_data(args.df_engine)

orchestrator.create_samplesheet()

print("Imputation and QC pipeline completed successfully.")
