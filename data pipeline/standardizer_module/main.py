from standardizer_classes import *
import argparse
import os
from config import TEST_FILE, ACCEPTED_VENDORS_DICT, GENOME_BUILD_DICT, CHAIN_FILE_DICT, GRCH_TO_HG_IDENTIFIER_DICT, FORWARD_STRAND_VENDORS

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=False, default=os.getcwd())
parser.add_argument('--microarray_file', type=str, required=False, default = TEST_FILE, help='Path to the user microarray data file (e.g., 23andMe, Ancestry, Myheritage...).')
parser.add_argument(
    '--identifier',
    type=str,
    help='A unique identifier for the current standardization run (e.g., user ID or timestamp).',
    default=f"{pd.Timestamp.now().strftime('%Y%m%d%H%M%S')}-{os.urandom(8).hex()}",
)

args = parser.parse_args()


data_container = DataContainer(identifier=args.identifier)

environment_handler = EnvironmentHandler(
    chain_file_dict=CHAIN_FILE_DICT,
    grch_to_hg_identifier_dict=GRCH_TO_HG_IDENTIFIER_DICT,
    user_file = args.microarray_file,
    output_dir = args.output_dir
)


file_handler = FileHandler(
    user_file = environment_handler.user_file,
    accepted_vendors_dict=ACCEPTED_VENDORS_DICT,
    genome_build_dict=GENOME_BUILD_DICT,
    forward_strand_vendors=FORWARD_STRAND_VENDORS
)


workflow_orchestrator = WorkflowOrchestrator(
    data_container=data_container,
    file_handler=file_handler,
    environment_handler=environment_handler
)

workflow_orchestrator.environment_handler.validate_environment()
workflow_orchestrator.check_dictionary_coherence()
workflow_orchestrator.evaluate_zipping()
workflow_orchestrator.set_genome_build()
workflow_orchestrator.set_vendor()
workflow_orchestrator.set_strand_direction()
workflow_orchestrator.set_microarray_data()
workflow_orchestrator.evaluate_liftover()
workflow_orchestrator.write_parquet_output()
workflow_orchestrator.write_meta_data_output()


