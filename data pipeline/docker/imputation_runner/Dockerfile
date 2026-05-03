FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"
ENV BASE_OUTPUT_DIR=/srv/imputed
ENV NEXTFLOW_BIN=/usr/local/bin/nextflow
ENV SAMPLESHEET_PATH=/app/imputation_runner_module/samplesheet.csv
ENV PIPELINE_PATH=/app/nextflow/pipeline_orchestrator.nf
ENV NEXTFLOW_CONFIG=/app/nextflow/nextflow.config
ENV NXF_HOME=/srv/imputed/nf_runtime/.nextflow
ENV NXF_WORK=/srv/imputed/nf_runtime/nf_work

RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre \
    curl \
    bash \
    ca-certificates \
    docker-cli \
    git \
    procps \
    unzip \
    gcc \
    libpq-dev \
 && rm -rf /var/lib/apt/lists/*

RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/nextflow \
 && chmod +x /usr/local/bin/nextflow

WORKDIR /app/imputation_runner_module

COPY database_module/requirements.txt /tmp/database_requirements.txt
COPY imputation_runner_module/requirements.txt /tmp/imputation_runner_requirements.txt

RUN pip install --no-cache-dir -r /tmp/database_requirements.txt \
 && pip install --no-cache-dir -r /tmp/imputation_runner_requirements.txt

COPY database_module/ /app/database_module/
COPY imputation_runner_module/ /app/imputation_runner_module/
COPY nextflow/ /app/nextflow/

RUN mkdir -p /srv/imputed/nf_runtime

CMD ["python", "/app/imputation_runner_module/main.py"]