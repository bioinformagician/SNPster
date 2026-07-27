# Use Python base image
FROM python:3.10-slim

# Set working directory
WORKDIR /app


# Copy all project files from harmonizing_module
COPY file_combiner_module/ .
RUN pip install -r requirements.txt

ENV MAXIMUM_OVERLAP_PERCENTAGE=5

# Run your script
CMD ["python", "main.py"]
