FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    bzip2 \
    git \
    build-essential \
    python3 \
    python3-pip \
    bcftools \
    tabix \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/python3 /usr/bin/python

# Install PLINK 2.0 (latest stable)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && \
    unzip plink2_linux_x86_64_latest.zip && \
    mv plink2 /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm plink2_linux_x86_64_latest.zip

# Install ADMIXTURE
RUN wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz && \
    tar -xvf admixture_linux-1.3.0.tar.gz && \
    mv dist/admixture_linux-1.3.0/admixture /usr/local/bin/ && \
    chmod +x /usr/local/bin/admixture && \
    rm -rf dist admixture_linux-1.3.0.tar.gz

# Install Python dependencies
COPY ancestry_module/requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

COPY database_module/requirements.txt /tmp/database_requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/database_requirements.txt



# Copy module code
WORKDIR /app
COPY ancestry_module/ /app/
COPY database_module/ /app/


#mount /srv/dependencies/imputation_runner/imputer/beagle_references$ as data
ENV REFERENCE_VCF_DIR=/data
ENV POPULATION_PANEL_FILE=/data/1000G_panel.txt 
ENV ANCESTRY_CHROMOSOME=22
ENV K_POPULATIONS=10
ENV ANCESTRY_METHOD=ADMIXTURE
ENV REFERENCE_PANEL=HGDP+1kGP

CMD ["python3", "main.py"]