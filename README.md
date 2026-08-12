# SNPster (active codebase)

This repository contains the active SNPster implementation for genetic data ingestion, imputation, ancestry estimation, polygenic score calculation, and report generation.

The `Old_POC/` directory is intentionally excluded from this README; the sections below describe the current pipeline modules and runtime behavior under `data pipeline/` and `webapp/`.

## 1) System structure (current + planned cloud architecture)

SNPster is split into two main planes:

- **Application plane (web app)**: a Django application under `webapp/demo/` that represents the user-facing layer.
- **Compute plane (data pipeline)**: containerized, CPU-heavy jobs under `data pipeline/` executed locally via Docker + Nextflow.

Planned production deployment direction:

- **Django web app on Azure** (future): user signup, file upload orchestration, report ordering, and status display.
- **Azure-hosted PostgreSQL** (future): central state store for jobs, metadata, and results.
- **S3 object storage for user files** (future): raw uploads and derived artifacts moved from local mounted paths to object storage.
- **Local or dedicated compute workers for heavy genomics tasks**: Nextflow fan-out across containerized stages remains the main strategy for expensive processing.

Current development/runtime reality in this repo:

- Pipeline workers use a PostgreSQL container from `data pipeline/docker/docker-compose.yaml`.
- The Django demo project currently uses SQLite by default (`webapp/demo/demo/settings.py`).
- Imputation and PGS workers poll job tables and execute independently as long-running services.

## 2) End-to-end job flow

At a high level, the lifecycle is:

1. A user signs up and uploads one or more genotype files.
2. The application records files and creates an `imputation_job` (status `queued`).
3. The imputation runner polls for queued jobs and executes a Nextflow orchestration.
4. Per-chromosome outputs are quality-checked and stored.
5. When a user orders a report, the app creates one or more `prsc_jobs` linked to that imputation.
6. The PGS worker polls for eligible queued `prsc_jobs` (only where imputation is complete), runs scoring, and uploads final report inputs/results.

In this codebase, development helpers for this flow are in `data pipeline/util_dev/setup_user.py`, which can create:

- `user_information`
- `user_files`
- `imputation_jobs` + `imputation_job_parameters`
- one queued `prsc_job` per report type + `prsc_job_parameters`

## 3) Imputation pipeline details

### 3.1 Queueing and orchestration

The long-running worker in `data pipeline/imputation_runner_module/main.py` repeatedly:

- validates runtime dependencies (Docker CLI and mounted Docker socket),
- queries queued jobs from `snpster_users.imputation_jobs`,
- writes a Nextflow samplesheet,
- runs the pipeline, and
- marks jobs `completed` or `failed` based on output validation.

The orchestration entrypoint is:

- `data pipeline/nextflow/pipeline_orhcestrator.nf`

Each major stage runs in its own container image:

- `standardizer`
- `file_combiner`
- `harmonizer`
- `vcf_merger`
- `ancestry`
- `imputer`
- `vcf_splitter`
- `imputation_qc`

### 3.2 Per-chromosome processing model

The workflow standardizes and harmonizes microarray inputs, then groups VCFs by chromosome, merges within chromosome, imputes, and splits/validates outputs for downstream use.

This chromosome-wise decomposition is central to scaling: it keeps each task unit smaller, parallelizable, and easier to retry on failure.

### 3.3 Ancestry estimation (ADMIXTURE)

The ancestry module (`data pipeline/ancestry_module/`) runs on merged chromosome data (chr22 in the current workflow for speed) and performs supervised ADMIXTURE-based inference against 1000 Genomes panel data.

Scientific reference for ADMIXTURE:

- Alexander DH, Novembre J, Lange K. **Fast model-based estimation of ancestry in unrelated individuals**. *Genome Research*. 2009;19(9):1655-1664. doi:10.1101/gr.094052.109

## 4) Polygenic score pipeline (PRSC jobs)

After imputation finishes, report ordering creates queued `prsc_jobs` linked to the completed imputation(s). The PGS worker in `data pipeline/pgs_calc_module/main.py` then:

1. Selects eligible queued PRSC jobs where imputation status is `completed`.
2. Batches compatible jobs (same score set) for compute efficiency.
3. Builds chromosome-aware merge samplesheets from imputed QC outputs.
4. Runs a Nextflow VCF merging step.
5. Runs polygenic scoring with the PGS Catalog calculator pipeline (`pgsc-calc`).
6. Validates outputs, uploads `prsc_job_results`, and marks job state.

The scoring stage delegates to the `pgsc-calc` Nextflow pipeline (`/opt/pgsc_calc/main.nf` inside the pgs container), with implementation details documented here:

- https://pgsc-calc.readthedocs.io/

In other words, SNPster wraps the pgsc-calc workflow inside its own job orchestration and database-driven scheduling model.

## 5) Folder guide (active parts)

- `webapp/demo/`: Django app (user-facing layer, currently minimal demo scaffolding)
- `data pipeline/database_module/`: shared PostgreSQL utilities and data I/O helpers
- `data pipeline/imputation_runner_module/`: DB polling + Nextflow launcher for imputation
- `data pipeline/nextflow/`: orchestration workflows/configs for containerized genomics stages
- `data pipeline/ancestry_module/`: ancestry inference and DB upload
- `data pipeline/pgs_calc_module/`: PRSC queue handling + pgsc-calc orchestration
- `data pipeline/docker/`: Dockerfiles and compose stack used by workers
- `data pipeline/reporting_module/`: report construction templates and logic

## 6) Notes on evolution toward cloud

As SNPster moves to Azure + object storage, the core architecture can remain stable if:

- job state continues to be DB-driven,
- heavy compute remains isolated behind worker services,
- storage paths are abstracted from local mounts to object URIs, and
- pipeline stages stay containerized with explicit inputs/outputs.

That keeps the scientific compute workflow reproducible while letting the web platform scale independently.
