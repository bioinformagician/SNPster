import os

BASE_OUTPUT_DIR = os.getenv("BASE_OUTPUT_DIR", "/srv/imputed")
IMPUTATION_DEPENDENCIES = os.getenv("IMPUTATION_DEPENDENCIES", "/home/frederik/shared_drive/snpster_dependencies/imputer_dependencies")
HARMONIZER_DEPENDENCIES = os.getenv("HARMONIZER_DEPENDENCIES", "/home/frederik/shared_drive/snpster_dependencies/harmonizer_dependencies")
NEXTFLOW_BIN = os.getenv("NEXTFLOW_BIN", "/home/frederik/.local/bin/nextflow")
SAMPLESHEET_PATH = os.getenv("SAMPLESHEET_PATH", "/home/frederik/snpster_project/samplesheets/samplesheet.csv")