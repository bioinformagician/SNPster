from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator
from config import OUTPUT_DIR, REFERENCE_DATA_PATH
import argparse

#adjust to new pgs_classes.py file



parser = argparse.ArgumentParser()
parser.add_argument('--low_memory', type=str, required=False, default="false") #setting it true only saved about 4gb but increases computational time x3 min

args = parser.parse_args()

environment_handler = EnvironmentHandler(
    output_dir=OUTPUT_DIR,
    low_memory = args.low_memory,
    reference_data_path=REFERENCE_DATA_PATH,
)

pgs_calculator_config = PGSCalculator_Config(
    environment_handler=environment_handler)

pgs_calculator = PGSCalculator(
    environment_handler=environment_handler,
    pgscalculator_config=pgs_calculator_config
)

environment_handler.connect_to_db()
pgs_calculator.run_pgs_calculation()
environment_handler.close_db_connection()
