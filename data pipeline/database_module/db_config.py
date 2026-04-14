import os


USERNAME = os.getenv("DB_USER", os.getenv("POSTGRES_USER", "postgres"))
PASSWORD = os.getenv("DB_PASSWORD", os.getenv("POSTGRES_PASSWORD", "zod50902"))
DATABASE_NAME = os.getenv("DB_NAME", os.getenv("POSTGRES_DB", "snpster_db"))
HOST = os.getenv("DB_HOST", "127.0.0.1")
PORT = int(os.getenv("DB_PORT", "5433"))
PGS_EXCEL_FILEPATH = os.getenv("PGS_EXCEL_FILEPATH", r"/home/frederik/snpster_project/pgs_all_metadata.xlsx")