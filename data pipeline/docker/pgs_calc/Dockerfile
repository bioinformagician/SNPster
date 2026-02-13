# Use Python base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install Java (required for Nextflow), Singularity and other dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-21-jre-headless \
    curl \
    procps \
    wget \
    squashfs-tools \
    cryptsetup \
    uidmap \
    libseccomp-dev && \
    rm -rf /var/lib/apt/lists/*

# singularity is used to avoid docker in docker issues  
RUN apt-get update && apt-get install -y --no-install-recommends \
    apptainer \
 && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow


RUN mkdir -p /app/pgs_output

#set environment variables, need to mount volume as :/data
ENV REFERENCE_DATA_PATH=/data/pgsc_HGDP+1kGP_v1.tar.zst
ENV PGS_ID_FILE=/data/pgs_ids.txt
ENV SAMPLESHEET_PATH=/data/samplesheet.csv
ENV OUTPUT_DIR=/app/pgs_output



#copy module
COPY pgs_calc_module/ .

CMD ["python", "main.py"]
