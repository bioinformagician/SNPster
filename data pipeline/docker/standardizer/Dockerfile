# Use Python base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app


# ---- Paths *inside the container* ----
# These three are baked into the image
ENV HG18_TO_HG38_CHAIN_FILE=/app/dependencies/hg18ToHg38.over.chain.gz
ENV HG19_TO_HG38_CHAIN_FILE=/app/dependencies/hg19ToHg38.over.chain.gz
ENV TEST_FILE=/app/dependencies/test_file/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt
ENV PVAR_REF_FILE=/data/all_hg38.pvar.zst

# Install system dependencies + PLINK2
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    ca-certificates \
    && wget -O /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip \
    && unzip /tmp/plink2.zip -d /tmp/plink2 \
    && mv /tmp/plink2/plink2 /usr/local/bin/plink2 \
    && chmod +x /usr/local/bin/plink2 \
    && rm -rf /tmp/plink2.zip /tmp/plink2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python requirements
COPY standardizer_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy all project files from harmonizing_module
COPY standardizer_module/ .
# Run your script
CMD ["python", "main.py"]
