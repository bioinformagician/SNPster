FROM pgscatalog/pgsc_calc:v2-blob

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

ARG BCFTOOLS_VERSION=1.22
ARG HTSLIB_VERSION=1.22

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    bash \
    ca-certificates \
    gcc \
    make \
    autoconf \
    automake \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libpq-dev \
    tabix \
 && rm -rf /var/lib/apt/lists/*

# Install htslib
RUN curl -L "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
      -o /tmp/htslib.tar.bz2 \
 && tar -xjf /tmp/htslib.tar.bz2 -C /tmp \
 && cd /tmp/htslib-${HTSLIB_VERSION} \
 && ./configure --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make install \
 && rm -rf /tmp/htslib*

# Install bcftools
RUN curl -L "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
      -o /tmp/bcftools.tar.bz2 \
 && tar -xjf /tmp/bcftools.tar.bz2 -C /tmp \
 && cd /tmp/bcftools-${BCFTOOLS_VERSION} \
 && ./configure --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make install \
 && rm -rf /tmp/bcftools*

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
ENV PYTHONUNBUFFERED=1

RUN bcftools --version \
 && bcftools merge --help 2>&1 | grep -q -- "--no-index" \
 && nextflow -version

RUN MATCH_INTERSECT=/opt/pgscatalog_utils/lib/python3.12/site-packages/pgscatalog/match/cli/intersect_cli.py \
 && grep -q 'if "bim" in pathlib.Path(path).name:' "$MATCH_INTERSECT" \
 && sed -i 's/if "bim" in pathlib.Path(path).name:/if pathlib.Path(path).name.endswith((".bim", ".bim.gz", ".bim.zst")):/' "$MATCH_INTERSECT" \
 && grep -q 'name.endswith((".bim", ".bim.gz", ".bim.zst"))' "$MATCH_INTERSECT"

WORKDIR /app

ENV REFERENCE_DATA_PATH=/data/pgsc_HGDP+1kGP_v1.tar.zst
ENV OUTPUT_DIR=/app/pgs_output
ENV SCORING_FILE_SOURCE_DIR=/data/scoring_files
ENV SCORING_FILE_TARGET_DIR=/app/scoring_files
ENV NF_WORK_DIR=/app/work
ENV BCFTOOLS_THREADS=8
ENV SAMPLESET_NAME=MergedUserData
ENV VCF_MERGE_SHEET_DIR=/app/vcf_merge_sheets
ENV PGS_TOTAL_CPUS=11
ENV PGS_TOTAL_MEM_GB=110
ENV N_PRSC_JOBS=303
ENV PGS_RESULT_DIR=/app/pgs_results
ENV NEXTFLOW_PGS_CONFIG=/app/nextflow_pgs.config
ENV NEXTFLOW_PGS_RESOURCE_CONFIG=/app/nextflow_pgs_resources.config
ENV NEXTFLOW_VCF_MERGING_CONFIG=/app/nextflow_vcf_merging.config
ENV NF_REPORT_DIR=/app/nf_reports

RUN mkdir -p "$OUTPUT_DIR" "$SCORING_FILE_TARGET_DIR" "$NF_WORK_DIR" "$VCF_MERGE_SHEET_DIR" "$PGS_RESULT_DIR" "$NF_REPORT_DIR"

COPY pgs_calc_module/requirements.txt /tmp/pgs_calc_requirements.txt
COPY database_module/requirements.txt /tmp/database_requirements.txt
COPY vcf_utilities_module/requirements.txt /tmp/vcf_utilities_requirements.txt

RUN pip install --no-cache-dir -r /tmp/pgs_calc_requirements.txt \
 && pip install --no-cache-dir -r /tmp/database_requirements.txt \
 && pip install --no-cache-dir -r /tmp/vcf_utilities_requirements.txt

COPY pgs_calc_module/ .
COPY database_module/ .
COPY vcf_utilities_module/ .

CMD ["python", "main.py"]