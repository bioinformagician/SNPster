import os
from config import PLINK_PREFIX, TEST_FILE, PLINK_1_9_PATH, PLINK_2_0_PATH, PVAR_REF_FILE, PLINK_REFERENCE_FASTA, PLINK_1_9_MEMORY_MB, ACCEPTED_VENDORS_DICT, GENOME_BUILD_DICT, PLINK_1_9_THREADS
from harmonizer_classes import EnvironmentHandler, DataContainer, FileHandler, WorkflowOrchestrator
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
    pvar_ref_file=PVAR_REF_FILE,
    plink_1_9_memory_mb=PLINK_1_9_MEMORY_MB,
    plink_1_9_threads=PLINK_1_9_THREADS
)


file_handler = FileHandler(
    user_file=args.microarray_file,
    accepted_vendors_dict=ACCEPTED_VENDORS_DICT,
    genome_build_dict=GENOME_BUILD_DICT
)

data_container = DataContainer()

workflow_orchestrator = WorkflowOrchestrator(
    environment_handler=environment_handler,
    data_container=data_container,
    file_handler=file_handler
)


#workflow_orchestrator.initiate_working_directory()
print(f"Working directory initialized at: {workflow_orchestrator.environment_handler.output_dir}")

workflow_orchestrator.set_vendor()
workflow_orchestrator.set_genome_build()
workflow_orchestrator.set_microarray_data()

print("Microarray data set in data container:")
print(workflow_orchestrator.data_container.microarray_data.head())


workflow_orchestrator.run_harmonization_workflow() #harmonize data if is_forward_strand is false (not implemented yet) or skip if true

print("Splitting harmonized data into chromosome-specific files...")
workflow_orchestrator.create_harmonized_chromosome_files()
print("Splitting complete. Chromosome-specific file paths:")

print("Converting harmonized files to BED format...")
workflow_orchestrator.convert_23andme_to_bed()
workflow_orchestrator.confirm_paths_exist(workflow_orchestrator.environment_handler.bed_file_paths)
print("Conversion to BED format complete. BED file paths:")

print("Converting BED files to VCF format...")
workflow_orchestrator.convert_bed_to_vcf()
workflow_orchestrator.confirm_paths_exist(workflow_orchestrator.environment_handler.vcf_file_paths)



