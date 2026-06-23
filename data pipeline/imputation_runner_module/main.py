from helper_classes import EnvironmentHandler, DatabaseQueryHandler, ImputationRunner
import time


environment_handler = EnvironmentHandler()
database_query_handler = DatabaseQueryHandler()
imputation_runner = ImputationRunner(env_handler=environment_handler,
                                   query_handler=database_query_handler)



while True:

    imputation_runner.env_handler.validate_environment()

    with imputation_runner.query_handler.db_utils.db_handler:
        imputation_runner.set_job_df()
    
    if imputation_runner.job_df.empty:
        print("No pending imputation jobs found. Retrying in 5 minutes...")
        imputation_runner.env_handler.clear_environment()
        time.sleep(300)  # Wait for 5 minutes before checking again
        continue


    imputation_ids = imputation_runner.job_df["imputation_id"].tolist()

    imputation_runner.write_samplesheet()

    #imputation_runner.query_handler.db_utils.db_handler.connect()

    #imputation_runner.query_handler.mark_jobs_running(imputation_ids)

    #imputation_runner.query_handler.db_utils.db_handler.close()


    try:
        imputation_runner.run_imputation()

        

        completed_ids, failed_ids = imputation_runner.evaluate_results()

        with imputation_runner.query_handler.db_utils.db_handler:
            if completed_ids:
                imputation_runner.query_handler.mark_jobs_completed(completed_ids)
            if failed_ids:
                imputation_runner.query_handler.mark_jobs_failed(failed_ids)
    except Exception as e:
        print(f"Error occurred during imputation: {e}")

        # Even when Nextflow exits non-zero, some IDs may still have complete outputs.
        completed_ids, failed_ids = imputation_runner.evaluate_results()
        remaining_failed_ids = [
            imp_id for imp_id in imputation_ids
            if imp_id not in completed_ids and imp_id not in failed_ids
        ]

        with imputation_runner.query_handler.db_utils.db_handler:
            if completed_ids:
                imputation_runner.query_handler.mark_jobs_completed(completed_ids)
            to_fail = failed_ids + remaining_failed_ids
            if to_fail:
                imputation_runner.query_handler.mark_jobs_failed(to_fail)

        imputation_runner.env_handler.clear_environment()
        continue

    imputation_runner.env_handler.clear_environment()
    #check file output and if file is missing, set job to failed, imputation id 8 is an example of a failed one (like three files inside compressed folder)

