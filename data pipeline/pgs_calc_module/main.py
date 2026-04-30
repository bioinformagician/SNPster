from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator
from config import OUTPUT_DIR, REFERENCE_DATA_PATH
import argparse
import time
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


while True:
    environment_handler.connect_to_db()
    
    try:
        pgs_calculator.set_job_parameters()
    except ValueError as e:
        print(f"No pending jobs found: {e}")
        environment_handler.close_db_connection()
        time.sleep(60)
        continue
    #self.environment_handler.set_db_job_status("running") implement this when heartbeat function is implemented
    environment_handler.close_db_connection()
    pgs_calculator.run_pgs_calculation()
    
    try:
        pgs_calculator.validate_results() 
    except Exception as e:
        print(f"Validation failed: {e}")
        environment_handler.connect_to_db()
        environment_handler.set_db_job_status("failed")
        environment_handler.close_db_connection()
        environment_handler.clear_directories()
        continue  # Skip uploading results and move to the next job
    
    environment_handler.connect_to_db()
    pgs_calculator.upload_results()
    environment_handler.set_db_job_status("completed")  #maybe this should be a trigger function in db upon upload of data to prsc_job_results instead --- IGNORE ---
    environment_handler.close_db_connection()
    environment_handler.clear_directories()
    
