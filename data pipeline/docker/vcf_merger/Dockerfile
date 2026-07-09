# Use Python base image
FROM python:3.10

# cant use multiple FROM statements in a single Dockerfile, so instaall java on top of python base image
RUN apt-get update && \
    apt-get install -y default-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre \
    beagle \
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


RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*



# Set working directory
WORKDIR /app


# Copy and install Python requirements
COPY vcf_utilities_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy Python files (exclude large dependencies)
COPY vcf_utilities_module/vcf_classes.py .
COPY vcf_utilities_module/vcf_combiner.py .
COPY vcf_utilities_module/vcf_combiner_pipeline.nf .
COPY vcf_utilities_module/merging_test_data /app/merging_test_data


# Run your script
CMD ["python", "vcf_combiner.py"]
