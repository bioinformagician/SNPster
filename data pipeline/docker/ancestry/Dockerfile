FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------------------------------------------------------
# System dependencies
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    bzip2 \
    git \
    build-essential \
    python3 \
    python3-pip \
    bcftools \
    tabix \
    openjdk-17-jre-headless \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/python3 /usr/bin/python


# -----------------------------------------------------------------------------
# PLINK 2.0
# -----------------------------------------------------------------------------
RUN wget -q \
    https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip \
    -O /tmp/plink2.zip \
    && unzip /tmp/plink2.zip -d /tmp/plink2 \
    && mv /tmp/plink2/plink2 /usr/local/bin/plink2 \
    && chmod +x /usr/local/bin/plink2 \
    && rm -rf /tmp/plink2 /tmp/plink2.zip


# -----------------------------------------------------------------------------
# ADMIXTURE 1.3.0
# -----------------------------------------------------------------------------
RUN wget -q \
    http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
    -O /tmp/admixture.tar.gz \
    && tar -xzf /tmp/admixture.tar.gz -C /tmp \
    && mv /tmp/dist/admixture_linux-1.3.0/admixture /usr/local/bin/admixture \
    && chmod +x /usr/local/bin/admixture \
    && rm -rf /tmp/dist /tmp/admixture.tar.gz


# -----------------------------------------------------------------------------
# FLARE
#
# Current Browning Lab version: 0.6.0
# Requires Java >= 11
# -----------------------------------------------------------------------------
ENV FLARE_HOME=/opt/flare
ENV FLARE_JAR=/opt/flare/flare.jar

# Default Java memory allocation.
# Can be overridden from Docker/Nextflow:
#
#   -e FLARE_JAVA_OPTS="-Xmx64g"
#
ENV FLARE_JAVA_OPTS="-Xmx32g"

RUN mkdir -p ${FLARE_HOME} \
    && wget -q \
        https://faculty.washington.edu/browning/flare.jar \
        -O ${FLARE_JAR} \
    && test -s ${FLARE_JAR}


# Convenience wrapper:
#
# Instead of:
#   java -Xmx32g -jar /opt/flare/flare.jar ...
#
# you can run:
#   flare ...
#
RUN printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -e' \
    'exec java ${FLARE_JAVA_OPTS} -jar ${FLARE_JAR} "$@"' \
    > /usr/local/bin/flare \
    && chmod +x /usr/local/bin/flare


# -----------------------------------------------------------------------------
# Python dependencies
# -----------------------------------------------------------------------------
COPY ancestry_module/requirements.txt /tmp/requirements.txt

RUN pip3 install \
    --no-cache-dir \
    -r /tmp/requirements.txt

    RUN pip3 install --no-cache-dir scikit-learn

COPY database_module/requirements.txt /tmp/database_requirements.txt

RUN pip3 install \
    --no-cache-dir \
    -r /tmp/database_requirements.txt


# -----------------------------------------------------------------------------
# Application
# -----------------------------------------------------------------------------
WORKDIR /app

COPY ancestry_module/ /app/
COPY database_module/ /app/


# -----------------------------------------------------------------------------
# Default environment
# -----------------------------------------------------------------------------

# Mounted HGDP+1kGP reference directory
ENV REFERENCE_VCF_DIR=/data

ENV POPULATION_PANEL_FILE=/data/1000G_panel.txt

ENV ANCESTRY_CHROMOSOME=22
ENV K_POPULATIONS=10
ENV ANCESTRY_METHOD=ADMIXTURE
ENV REFERENCE_PANEL=HGDP+1kGP


CMD ["python3", "main.py"]