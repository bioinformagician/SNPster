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
        print("No eligible PGS jobs found. Requires queued PRSC jobs with completed imputations and matching prsc_job_parameters rows.")
        environment_handler.close_db_connection()
        time.sleep(60)
        continue
    
    #self.environment_handler.set_db_job_status("running") implement this when heartbeat function is implemented
    
    #create samplesheets for each imputed dir
    
    pgs_calculator.create_samplesheets()
    #pgs_calculator.create_chr_specific_samplesheet()
    
    
    environment_handler.close_db_connection()
    print(environment_handler.imputation_ids)
    if len(environment_handler.imputation_ids) > 1:
        
        environment_handler.create_full_path_samplesheets() #this creates one sheet per chrom in the vcf_merge_sheet_dir with columns: sampleset, path_prefix, chrom, format for the nextflow process to merge them fast
        #pgs_calculator.create_chr_specific_samplesheet()
        pgs_calculator.run_vcf_merging()
        
        
        environment_handler.create_merged_vcf_samplesheet()
    
        pgs_calculator.run_pgs_calculation(environment_handler.merged_sample_sheet)
        pgs_score_path = f"{environment_handler.output_dir}/{environment_handler.sampleset_name}/score/{environment_handler.sampleset_name}_pgs.txt.gz"

    else:
        pgs_calculator.run_pgs_calculation(environment_handler.samplesheet_paths[0])
        pgs_score_path = f"{environment_handler.output_dir}/{environment_handler.imputation_ids[0]}/score/{environment_handler.imputation_ids[0]}_pgs.txt.gz"
        pgs_scoring_summary_path = f"{environment_handler.output_dir}/{environment_handler.imputation_ids[0]}/match/{environment_handler.imputation_ids[0]}_summary.csv"
    
    
    results_valid = pgs_calculator.validate_results(pgs_score_path)
    results_valid_summary = pgs_calculator.validate_results(pgs_scoring_summary_path)
    
    if not results_valid or not results_valid_summary:
        print("Validation failed.")
        environment_handler.connect_to_db()
        environment_handler.set_db_job_status("failed")
        environment_handler.close_db_connection()
        environment_handler.clear_directories()
        continue  # Skip uploading results and move to the next job
    
    environment_handler.connect_to_db()
    pgs_calculator.upload_results(path = pgs_score_path)
    
    environment_handler.set_db_job_status("completed")  #maybe this should be a trigger function in db upon upload of data to prsc_job_results instead --- IGNORE ---
    environment_handler.close_db_connection()
    environment_handler.clear_directories()
    
    
