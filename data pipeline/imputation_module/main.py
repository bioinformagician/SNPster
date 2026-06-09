import os
import argparse
from imputation_classes import EnvironmentHandler, ImputationRunner
from config import BEAGLE_JAR, JAVA_EXE, HEAP_GB, BEAGLE_REFERENCE_DIR, PLINK_MAP_DIR, OUTPUT_DIR, THREADS

parser = argparse.ArgumentParser()
parser.add_argument('--vcf_file', type=str, required=True) #default is vcf files for testing baked into image
parser.add_argument('--output_dir', type=str, required=False, default=OUTPUT_DIR) #default is output dir for testing baked into image

args = parser.parse_args()

environment_handler = EnvironmentHandler(
    working_dir=os.getcwd(),
    java_exe=JAVA_EXE,
    beagle_jar=BEAGLE_JAR,
    heap_gb=HEAP_GB,
    beagle_threads=THREADS,
    output_dir=args.output_dir,
    beagle_reference_dir=BEAGLE_REFERENCE_DIR, #mounted
    plink_map_dir=PLINK_MAP_DIR, #mounted
    vcf_file=args.vcf_file
)



orchestrator = ImputationRunner(
    environment_handler=environment_handler
)


orchestrator.create_vcf_reference_mapping()


orchestrator.impute_vcf_file()

