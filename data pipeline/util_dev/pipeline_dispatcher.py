import subprocess
import sys
import time 
import os

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))
print(sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module")))
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH









if __name__ == "__main__":

    while True:

        db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
        db_handler.connect()

        query = f"""SELECT ij.user_id, ij.imputation_id, ui.genefile_location 
                    FROM snpster_users.user_information ui
                    JOIN snpster_users.imputation_jobs ij ON ui.user_id = ij.user_id
                    WHERE ij.imputation_status = 'queued'
                    LIMIT 10;"""
        
        results = db_handler.execute_query(query)
        for user_id, imputation_id, genefile_location in results:
            print(f"Processing user_id: {user_id}, imputation_id: {imputation_id}, genefile_location: {genefile_location}")
            # Here you would call your processing function, e.g.:
            # process_genefile(genefile_location)
            # After processing, update the status in the database
            update_query = f"""UPDATE snpster_users.imputation_jobs
                               SET imputation_status = 'processed' 
                               WHERE user_id = '{user_id}' AND imputation_id = '{imputation_id}';"""
            db_handler.execute_query(update_query)