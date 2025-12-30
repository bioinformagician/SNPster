# Use Python base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app


# ---- Paths *inside the container* ----
# These three are baked into the image
ENV HG18_TO_HG38_CHAIN_FILE=/app/dependencies/hg18ToHg38.over.chain.gz
ENV HG19_TO_HG38_CHAIN_FILE=/app/dependencies/hg19ToHg38.over.chain.gz
ENV TEST_FILE=/app/dependencies/test_file/23andMe.zip



# Copy and install Python requirements
COPY standardizer_module/requirements.txt .
RUN pip install -r requirements.txt

# Copy all project files from harmonizing_module
COPY standardizer_module/ .
# Run your script
CMD ["python", "main.py"]
