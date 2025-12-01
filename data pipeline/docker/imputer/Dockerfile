# Use Python base image
FROM python:3.9-slim

# cant use multiple FROM statements in a single Dockerfile, so instaall java on top of python base image
RUN apt-get update && \
    apt-get install -y default-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy and install Python requirements
COPY imputation_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy all project files from imputation_module
COPY imputation_module/ .

# Create directories for volume mounts
RUN mkdir -p /data /work /output

# ---- Paths *inside the container* ----
# Beagle JAR is baked into the image at /app/dependencies/
ENV BEAGLE_JAR=/app/dependencies/beagle.27Feb25.75f.jar
ENV JAVA_EXE=/usr/bin/java
ENV HEAP_GB=8
ENV THREADS=2
ENV GP_MIN=0.90
ENV DS_TOL=0.05
ENV SNPS_ONLY=True
ENV BIALLELIC_ONLY=True

# These big files are NOT in the image.
# We expect them to appear under /data via a volume mount at runtime.
ENV PLINK_MAP_DIR=/data/plink.GRCh37.map
ENV PLINK_REFERENCE_FASTA=/data/hs37d5.fa.zst
ENV BEAGLE_REFERENCE_DIR=/data/beagle_references
ENV PVAR_REF_FILE=/data/all_phase3.pvar.zst

# Output directory for imputed results
ENV OUTPUT_DIR=/output

# Run your script
CMD ["python", "main.py"]
