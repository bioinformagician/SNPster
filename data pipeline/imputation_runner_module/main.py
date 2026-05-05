from helper_classes import EnvironmentHandler, DatabaseQueryHandler, ImputationRunner



environment_handler = EnvironmentHandler()
database_query_handler = DatabaseQueryHandler()
imputation_runner = ImputationRunner(env_handler=environment_handler,
                                   query_handler=database_query_handler)


max = 2
counter = 0

while counter < max:

    imputation_runner.env_handler.validate_environment()

    imputation_runner.query_handler.db_utils.db_handler.connect()

    imputation_runner.set_job_df()

    imputation_runner.query_handler.db_utils.db_handler.close()

    imputation_ids = imputation_runner.job_df["imputation_id"].tolist()

    imputation_runner.write_samplesheet()

    #imputation_runner.query_handler.db_utils.db_handler.connect()

    #imputation_runner.query_handler.mark_jobs_running(imputation_ids)

    #imputation_runner.query_handler.db_utils.db_handler.close()


    try:
        imputation_runner.run_imputation()
        imputation_runner.query_handler.db_utils.db_handler.connect()
        imputation_runner.query_handler.mark_jobs_completed(imputation_ids)
        imputation_runner.query_handler.db_utils.db_handler.close()
    except Exception as e:
        print(f"Error occurred during imputation: {e}")
        #imputation_runner.query_handler.db_utils.db_handler.connect()
        #imputation_runner.query_handler.mark_jobs_failed(imputation_ids)
        #imputation_runner.query_handler.db_utils.db_handler.close()

    imputation_runner.env_handler.clear_environment()
    #check file output and if file is missing, set job to failed, imputation id 8 is an example of a failed one (like three files inside compressed folder)

    counter += 1