from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator
import time


environment_handler = EnvironmentHandler()

pgs_calculator_config = PGSCalculator_Config(
    environment_handler=environment_handler)

pgs_calculator = PGSCalculator(
    environment_handler=environment_handler,
    pgscalculator_config=pgs_calculator_config)


while True:
    environment_handler.connect_to_db()
    
    
    jobs=pgs_calculator.set_job_parameters()
    
    if not jobs:
        print("No pending jobs found.")
        environment_handler.close_db_connection()
        time.sleep(60)
        continue
    
    #self.environment_handler.set_db_job_status("running") implement this when heartbeat function is implemented
    #pgs_calculator.copy_imputed_files()
    #pgs_calculator.create_combined_samplesheet()
    environment_handler.close_db_connection()
    
    
    for idx, samplesheet in enumerate(pgs_calculator.environment_handler.samplesheet_paths):
        
        print(f"Processing job {idx+1}/{len(pgs_calculator.environment_handler.samplesheet_paths)}: {samplesheet}")
        
        pgs_calculator.run_pgs_calculation(samplesheet)
    
        results_valid = pgs_calculator.validate_results(imputation_id=pgs_calculator.environment_handler.imputation_ids[idx]) 
        
        if not results_valid:
            print("Validation failed.")
            environment_handler.connect_to_db()
            environment_handler.set_db_job_status("failed", pgs_calculator.environment_handler.prsc_ids[idx])
            environment_handler.close_db_connection()
            
            continue  # Skip uploading results and move to the next job
        
        environment_handler.connect_to_db()
        pgs_calculator.upload_results(imputation_id=pgs_calculator.environment_handler.imputation_ids[idx],
                                      prsc_id=pgs_calculator.environment_handler.prsc_ids[idx])
        
        environment_handler.set_db_job_status("completed", pgs_calculator.environment_handler.prsc_ids[idx])  #maybe this should be a trigger function in db upon upload of data to prsc_job_results instead --- IGNORE ---
        environment_handler.close_db_connection()
        environment_handler.clear_output_directory()
    
    environment_handler.clear_directories()
    
