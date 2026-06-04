import os
from config import PLINK_PREFIX, TEST_FILE, PLINK_1_9_PATH, PLINK_2_0_PATH, PLINK_REFERENCE_FASTA, PLINK_1_9_MEMORY_MB, PLINK_1_9_THREADS
from harmonizer_classes import EnvironmentHandler, DataContainer, WorkflowOrchestrator
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=False, default=os.getcwd())
parser.add_argument('--microarray_file', type=str, required=False, default = TEST_FILE, help='Path to the user microarray data file (e.g., 23andMe, Ancestry, Myheritage...).')

args = parser.parse_args()

environment_handler = EnvironmentHandler(
    output_dir=args.output_dir,
    user_upload_file=args.microarray_file,
    plink_1_9_path=PLINK_1_9_PATH,
    plink_2_0_path=PLINK_2_0_PATH,
    PLINK_PREFIX=PLINK_PREFIX,
    plink_reference_fasta=PLINK_REFERENCE_FASTA,
    plink_1_9_memory_mb=PLINK_1_9_MEMORY_MB,
    plink_1_9_threads=PLINK_1_9_THREADS
)


data_container = DataContainer()

workflow_orchestrator = WorkflowOrchestrator(
    environment_handler=environment_handler,
    data_container=data_container
)


#workflow_orchestrator.initiate_working_directory()
print(f"Working directory initialized at: {workflow_orchestrator.environment_handler.output_dir}")


workflow_orchestrator.read_parquet()

print("Converting harmonized file to BED format...")
workflow_orchestrator.convert_23andme_to_bed()
workflow_orchestrator.confirm_path_exist(workflow_orchestrator.environment_handler.bed_file_path)
print("Conversion to BED format complete")

print("Converting BED file to VCF format...")
workflow_orchestrator.convert_bed_to_vcf()
workflow_orchestrator.confirm_path_exist(workflow_orchestrator.environment_handler.vcf_file_path)
workflow_orchestrator.add_imputation_id_to_vcfs()
print("Conversion to VCF format complete.")


