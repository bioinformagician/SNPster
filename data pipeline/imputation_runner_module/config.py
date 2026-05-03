import os

BASE_OUTPUT_DIR = os.getenv("BASE_OUTPUT_DIR", "/srv/imputed")
IMPUTATION_DEPENDENCIES = os.getenv("IMPUTATION_DEPENDENCIES", "/srv/dependencies/imputation_runner/imputer") #mounted
HARMONIZER_DEPENDENCIES = os.getenv("HARMONIZER_DEPENDENCIES", "/srv/dependencies/imputation_runner/harmonizer") #mounted
NEXTFLOW_BIN = os.getenv("NEXTFLOW_BIN", "/home/frederik/.local/bin/nextflow")
NEXTFLOW_CONFIG = os.getenv("NEXTFLOW_CONFIG", "/app/nextflow/nextflow.config")
NEXTFLOW_WORK_DIR = os.getenv("NXF_WORK", "/srv/imputed/nf_runtime/nf_work")
SAMPLESHEET_PATH = os.getenv("SAMPLESHEET_PATH", "/app/imputation_runner_module/samplesheet.csv")
PIPELINE_PATH = os.getenv("PIPELINE_PATH", "/app/nextflow/pipeline_orchestrator.nf")
