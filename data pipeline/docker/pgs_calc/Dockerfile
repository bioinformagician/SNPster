FROM pgscatalog/pgsc_calc:v2-blob
#this is the official base image for the pgs calculator pipeline
WORKDIR /app

ENV REFERENCE_DATA_PATH=/data/pgsc_HGDP+1kGP_v1.tar.zst
ENV OUTPUT_DIR=/app/pgs_output
ENV SCORING_FILE_SOURCE_DIR=/data/scoring_files
ENV SCORING_FILE_TARGET_DIR=/app/scoring_files
ENV NF_WORK_DIR=/app/work

RUN mkdir -p "$OUTPUT_DIR" "$SCORING_FILE_TARGET_DIR" "$NF_WORK_DIR"

COPY pgs_calc_module/requirements.txt /tmp/pgs_calc_requirements.txt
COPY database_module/requirements.txt /tmp/database_requirements.txt

RUN pip install --no-cache-dir -r /tmp/pgs_calc_requirements.txt \
 && pip install --no-cache-dir -r /tmp/database_requirements.txt

COPY pgs_calc_module/ .
COPY database_module/ .

CMD ["python", "main.py"]