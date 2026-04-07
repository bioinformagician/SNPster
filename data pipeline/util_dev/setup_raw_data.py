import os
import shutil
import zstandard as zstd
from pathlib import Path
import sys
import random

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module"))
print(sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "database_module")))
from db_handler import DbHandler, DbUtils
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH


RAW_DATA_DIR = "/home/frederik/snpster_project/zipped"
TARGET_DIR = "/srv/raw"


db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
db_handler.connect()

def transfer_files(source_dir:str, target_dir:str, db_handler:DbHandler=db_handler) -> None:
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)



    for _, filename in enumerate(os.listdir(source_dir)):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy2(source_file, target_file)
            print(f"Copied: {source_file} to {target_file}")

            user_id = _
            email = f"{_}@example.com"
            created_at = "NOW()"
            genefile_storage_backend = "local"
            genefile_location = target_file

            query = f"""
            INSERT INTO snpster_users.user_information (user_id, email, created_at, genefile_storage_backend, genefile_location)
            VALUES ('{user_id}', '{email}', {created_at}, '{genefile_storage_backend}', '{genefile_location}');
            """
            print(f"Executing query: {query}")
            db_handler.execute_query(query)
            print(f"Updated database with file location for {target_file}")



if __name__ == "__main__":
    transfer_files(RAW_DATA_DIR, TARGET_DIR)


"""CREATE TABLE snpster_users.user_information (
    user_id varchar(100) PRIMARY KEY,
    email VARCHAR(100) NOT NULL UNIQUE,
    password_hash VARCHAR(255) NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    genefile_location TEXT -- stored on linux server in folder /srv/raw
);"""