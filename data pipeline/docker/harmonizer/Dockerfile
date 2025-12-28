# Use Python base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app


# ---- Paths *inside the container* ----
# These three are baked into the image
ENV PLINK_1_9_PATH=/app/dependencies/plink
ENV PLINK_2_0_PATH=/app/dependencies/plink2
ENV TEST_FILE=/app/dependencies/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt
ENV PLINK_1_9_MEMORY_MB=1000
ENV PLINK_1_9_THREADS=1

# These big files are NOT in the image.
# We expect them to appear under /data via a volume mount at runtime.

ENV PLINK_PREFIX=plink_temp
ENV PLINK_REFERENCE_FASTA=/data/hs37d5.fa.zst
ENV PVAR_REF_FILE=/data/all_phase3.pvar.zst

# Copy and install Python requirements
COPY harmonizing_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy all project files from harmonizing_module
COPY harmonizing_module/ .

# Run your script
CMD ["python", "main.py"]
