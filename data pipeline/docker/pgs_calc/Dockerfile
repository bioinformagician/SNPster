FROM python:3.9-slim-bookworm

WORKDIR /app

# ---- system deps (Java for Nextflow, plus tools Apptainer commonly needs) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    ca-certificates curl wget procps \
    squashfs-tools cryptsetup uidmap libseccomp-dev \
    fuse3 libfuse3-3 \
 && rm -rf /var/lib/apt/lists/*

# ---- install Apptainer from official release .deb ----
ARG APPTAINER_VER=1.4.5
RUN set -eux; \
    cd /tmp; \
    wget -q https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VER}/apptainer_${APPTAINER_VER}_amd64.deb; \
    wget -q https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VER}/apptainer-suid_${APPTAINER_VER}_amd64.deb; \
    apt-get update; \
    apt-get install -y --no-install-recommends ./apptainer_${APPTAINER_VER}_amd64.deb ./apptainer-suid_${APPTAINER_VER}_amd64.deb; \
    rm -f /tmp/apptainer_*.deb; \
    rm -rf /var/lib/apt/lists/*

# ---- install Nextflow ----
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

RUN mkdir -p /app/pgs_output

ENV REFERENCE_DATA_PATH=/data/pgsc_HGDP+1kGP_v1.tar.zst
ENV PGS_ID_FILE=/data/pgs_ids.txt
ENV SAMPLESHEET_PATH=/data/samplesheet.csv
ENV OUTPUT_DIR=/app/pgs_output

COPY pgs_calc_module/ .

CMD ["python", "main.py"]