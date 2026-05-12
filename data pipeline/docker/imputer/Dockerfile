# Use Python base image
FROM python:3.10

# cant use multiple FROM statements in a single Dockerfile, so instaall java on top of python base image
RUN apt-get update && \
    apt-get install -y default-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    tabix \
    && rm -rf /var/lib/apt/lists/*


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


# Set working directory
WORKDIR /app


# Copy and install Python requirements
COPY imputation_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy Python files (exclude large dependencies)
COPY imputation_module/*.py .
COPY vcf_merging_module/vcf* .

# Copy only needed dependencies (Beagle JAR and test data)
COPY imputation_module/dependencies/beagle.27Feb25.75f.jar ./dependencies/
COPY imputation_module/dependencies/test_data ./dependencies/test_data



ENV BEAGLE_JAR=/app/dependencies/beagle.27Feb25.75f.jar
ENV JAVA_EXE=/usr/bin/java
ENV HEAP_GB=6
ENV THREADS=2
ENV GP_MIN=0.90
ENV DS_TOL=0.05
ENV SNPS_ONLY=True
ENV BIALLELIC_ONLY=True
ENV POLARS_MAX_THREADS=1
ENV DF_ENGINE=polars
ENV MERGED_SAMPLESHEET_DIR=/merged_samplesheets

# Create directories for volume mounts
RUN mkdir -p /data /work /output ${MERGED_SAMPLESHEET_DIR}
# ---- Paths *inside the container* ----
# Beagle JAR is baked into the image at /app/dependencies/



#otherwise polars will use all cores available

# These big files are NOT in the image.
#expect them to appear under /data via a volume mount at runtime.
ENV PLINK_MAP_DIR=/data/plink.GRCh38.map
ENV BEAGLE_REFERENCE_DIR=/data/beagle_references
ENV OUTPUT_DIR=/output


#default test data for testing
ENV VCF_FILES_DIR=/app/dependencies/test_data


# Run your script
CMD ["python", "main.py"]
