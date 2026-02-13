from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator
from config import SAMPLESHEET_PATH, OUTPUT_DIR, REFERENCE_DATA_PATH, PGS_ID_FILE
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--low_memory', type=str, required=False, default="false") #setting it true only saved about 4gb but increases computational time x3 min

args = parser.parse_args()

environment_handler = EnvironmentHandler(
    samplesheet_path=SAMPLESHEET_PATH,
    output_dir=OUTPUT_DIR,
    low_memory = args.low_memory,
    reference_data_path=REFERENCE_DATA_PATH,
    pgs_id_file=PGS_ID_FILE
)

pgs_calculator_config = PGSCalculator_Config(
    environment_handler=environment_handler,
    pgs_id_str=None)

pgs_calculator = PGSCalculator(
    environment_handler=environment_handler,
    pgscalculator_config=pgs_calculator_config
)

pgs_calculator.run_pgs_calculation()

