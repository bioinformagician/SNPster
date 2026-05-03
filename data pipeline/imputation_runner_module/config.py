import os

BASE_OUTPUT_DIR = os.getenv("BASE_OUTPUT_DIR", "/srv/imputed")
IMPUTATION_DEPENDENCIES = os.getenv("IMPUTATION_DEPENDENCIES", "/home/frederik/shared_drive/snpster_dependencies/imputer_dependencies")
HARMONIZER_DEPENDENCIES = os.getenv("HARMONIZER_DEPENDENCIES", "/home/frederik/shared_drive/snpster_dependencies/harmonizer_dependencies")
NEXTFLOW_BIN = os.getenv("NEXTFLOW_BIN", "/home/frederik/.local/bin/nextflow")
NEXTFLOW_CONFIG = os.getenv("NEXTFLOW_CONFIG", "/app/nextflow/nextflow.config")
SAMPLESHEET_PATH = os.getenv("SAMPLESHEET_PATH", "/app/job_runner_module/samplesheet.csv")
PIPELINE_PATH = os.getenv("PIPELINE_PATH", "/app/nextflow/pipeline_orchestrator.nf")