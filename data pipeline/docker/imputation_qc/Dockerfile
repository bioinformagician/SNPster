# Use Python base image
FROM python:3.10


RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
 && wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 \
 && tar -xjf htslib-1.17.tar.bz2 \
 && cd htslib-1.17 \
 && ./configure \
 && make \
 && make install \
 && cd .. \
 && rm -rf htslib-1.17 htslib-1.17.tar.bz2


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
COPY qc_module/requirements.txt .
RUN pip install -r requirements.txt
COPY vcf_utilities_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy Python files (exclude large dependencies)
COPY qc_module/ .
COPY vcf_utilities_module/vcf_classes.py .



ENV ENGINE=polars
ENV THREADS=1
ENV HEAP_GB=6
ENV THREADS=2
ENV GP_MIN=0.90
ENV DS_TOL=0.05
ENV SNPS_ONLY=True
ENV BIALLELIC_ONLY=True
ENV POLARS_MAX_THREADS=1
ENV OUTPUT_DIR=/output



# Create directories for volume mounts
RUN mkdir -p /data /work ${OUTPUT_DIR} ${MERGED_SAMPLESHEET_DIR}
# ---- Paths *inside the container* ----
# Beagle JAR is baked into the image at /app/dependencies/



#default test data for testing, should be mounted
ENV VCF_FILES_TEST_DIR=/app/test_data


CMD ["python", "main.py"]
