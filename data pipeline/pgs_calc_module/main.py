from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator, VCFHandler
import time


environment_handler = EnvironmentHandler()

pgs_calculator_config = PGSCalculator_Config(
    environment_handler=environment_handler)

pgs_calculator = PGSCalculator(
    environment_handler=environment_handler,
    pgscalculator_config=pgs_calculator_config)

vcf_handler = VCFHandler(environment_handler=environment_handler)


while True:
    environment_handler.connect_to_db()
    
    
    jobs=pgs_calculator.set_job_parameters()
    
    if not jobs:
        print("No pending jobs found.")
        environment_handler.close_db_connection()
        time.sleep(60)
        continue
    
    #self.environment_handler.set_db_job_status("running") implement this when heartbeat function is implemented
    
    pgs_calculator.create_combined_samplesheet()
    environment_handler.close_db_connection()
    
    if len(environment_handler.imputation_ids) > 1:
        vcf_handler.set_user_imputation_id_dict()
        vcf_handler.set_vcf_file_dict()
        vcf_handler.merge_vcf_files()
    
        pgs_calculator.run_pgs_calculation(environment_handler.merged_sample_sheet)
        pgs_score_path = f"{environment_handler.output_dir}/{environment_handler.sampleset_name}/score/{environment_handler.sampleset_name}_pgs.txt.gz"

    else:
        pgs_calculator.run_pgs_calculation(environment_handler.samplesheet_paths[0])
        pgs_score_path = f"{environment_handler.output_dir}/{environment_handler.imputation_ids[0]}/score/{environment_handler.imputation_ids[0]}_pgs.txt.gz"
    
    
    results_valid = pgs_calculator.validate_results(pgs_score_path) 
    
    if not results_valid:
        print("Validation failed.")
        environment_handler.connect_to_db()
        environment_handler.set_db_job_status("failed")
        environment_handler.close_db_connection()
        
        continue  # Skip uploading results and move to the next job
    
    environment_handler.connect_to_db()
    pgs_calculator.upload_results(path = pgs_score_path)
    
    environment_handler.set_db_job_status("completed")  #maybe this should be a trigger function in db upon upload of data to prsc_job_results instead --- IGNORE ---
    environment_handler.close_db_connection()
    environment_handler.clear_directories()
    
