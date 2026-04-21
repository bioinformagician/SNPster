FROM python:3.9-slim-bookworm

WORKDIR /app

# ---- system deps (Java for Nextflow, plus tools Apptainer commonly needs) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    ca-certificates curl wget procps \
    squashfs-tools cryptsetup uidmap libseccomp-dev \
    fuse3 libfuse3-3 \
 && rm -rf /var/lib/apt/lists/*

# apptainer/singularity so i can use pgs pipeline from prs catalog inside container (to avoid docker in docker issues)
ARG APPTAINER_VER=1.4.5
RUN set -eux; \
    cd /tmp; \
    wget -q https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VER}/apptainer_${APPTAINER_VER}_amd64.deb; \
    wget -q https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VER}/apptainer-suid_${APPTAINER_VER}_amd64.deb; \
    apt-get update; \
    apt-get install -y --no-install-recommends ./apptainer_${APPTAINER_VER}_amd64.deb ./apptainer-suid_${APPTAINER_VER}_amd64.deb; \
    rm -f /tmp/apptainer_*.deb; \
    rm -rf /var/lib/apt/lists/*

# nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

RUN mkdir -p /app/pgs_output

ENV REFERENCE_DATA_PATH=/data/pgsc_HGDP+1kGP_v1.tar.zst
ENV OUTPUT_DIR=/app/pgs_output

# Install both requirement sets with distinct file names to avoid overwrite.
COPY pgs_calc_module/requirements.txt /tmp/pgs_calc_requirements.txt
COPY database_module/requirements.txt /tmp/database_requirements.txt
RUN pip install --no-cache-dir -r /tmp/pgs_calc_requirements.txt \
 && pip install --no-cache-dir -r /tmp/database_requirements.txt

COPY pgs_calc_module/ .
COPY database_module/ .

CMD ["python", "main.py"]