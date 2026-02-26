import os

QMD_PATH = os.getenv("QMD_PATH", "data pipeline\\reporting_module\\prs_report.qmd")
OUT_DIR = os.getenv("OUT_DIR", "test/path/to/_rendered")
DATA_PATH = os.getenv("DATA_PATH", "test/path/to/data.parquet")
