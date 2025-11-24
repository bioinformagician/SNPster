import os
from config import PLINK_PREFIX, PLINK_MAP_DIR, BEAGLE_REFERENCE_DIR, TEST_FILE, PLINK_1_9_PATH, PLINK_2_0_PATH, PVAR_REF_FILE, PLINK_REFERENCE_FASTA
from harmonizer_classes import EnvironmentHandler, DataContainer, WorkflowOrchestrator

environment_handler = EnvironmentHandler(
    working_dir=os.getcwd(),
    user_upload_file=TEST_FILE,
    plink_1_9_path=PLINK_1_9_PATH,
    plink_2_0_path=PLINK_2_0_PATH,
    plink_map_dir=PLINK_MAP_DIR,
    PLINK_PREFIX=PLINK_PREFIX,
    plink_reference_fasta=PLINK_REFERENCE_FASTA,
    pvar_ref_file=PVAR_REF_FILE,
    beagle_references=BEAGLE_REFERENCE_DIR
)

environment_handler.set_beagle_files()

data_container = DataContainer()

workflow_orchestrator = WorkflowOrchestrator(
    
    environment_handler=environment_handler,
    data_container=data_container,
    working_dir=os.getcwd() #remove this later
)


#workflow_orchestrator.initiate_working_directory()
print(f"Working directory initialized at: {workflow_orchestrator.working_dir}")

workflow_orchestrator.set_microarray_data()
print("Microarray data set in data container:")
print(workflow_orchestrator.data_container.microarray_data.head())

workflow_orchestrator.create_user_snp_list()
print("User SNP list created.")

print("Extracting reference data...")
workflow_orchestrator.extract_reference_data()
workflow_orchestrator.data_container.reference_data = workflow_orchestrator.read_vcf_like_to_df()
print("Reference data extracted and set in data container")

print("Harmonizing data...")
workflow_orchestrator.data_container.harmonize_data()
print("Data harmonization complete. Harmonized data:")
print(workflow_orchestrator.data_container.harmonization_stats)
print(workflow_orchestrator.data_container.harmonized_data.head())

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
print("Conversion to VCF format complete. VCF file paths:")
print(workflow_orchestrator.environment_handler.vcf_file_paths)
print(workflow_orchestrator.environment_handler.split_harmonized_file_paths)


print("Creating VCF to reference mapping...")
vcf_reference_mapping_df = workflow_orchestrator.create_vcf_reference_mapping()
print("VCF to reference mapping created:")
print(vcf_reference_mapping_df)

print("Writing VCF to reference mapping to output file...")
os.makedirs(os.path.join(workflow_orchestrator.working_dir, "harmonization_results"), exist_ok=True)
vcf_reference_mapping_df.to_parquet(os.path.join(workflow_orchestrator.working_dir, "harmonization_results/vcf_reference_mapping.parquet"))
print("VCF to reference mapping written to output file.")


