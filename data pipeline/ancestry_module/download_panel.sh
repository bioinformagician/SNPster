#!/bin/bash

# Download 1000 Genomes Phase 3 population panel file
# This file maps sample IDs to populations and superpopulations

PANEL_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR="/srv/dependencies/imputation_runner"
OUTPUT_FILE="${OUTPUT_DIR}/1000G_panel.txt"

echo "Downloading 1000 Genomes population panel..."
echo "(Using sudo to write to ${OUTPUT_DIR})"

# Download to temp location first
TEMP_FILE="/tmp/1000G_panel.txt"
wget -q -O "${TEMP_FILE}" "${PANEL_URL}"

if [ -f "${TEMP_FILE}" ]; then
    # Move to final location with sudo
    sudo mv "${TEMP_FILE}" "${OUTPUT_FILE}"
    sudo chmod 644 "${OUTPUT_FILE}"
    
    echo "✓ Population panel downloaded successfully"
    echo "  Location: ${OUTPUT_FILE}"
    echo ""
    echo "Sample format:"
    head -n 3 "${OUTPUT_FILE}"
    echo ""
    echo "Total samples: $(tail -n +2 "${OUTPUT_FILE}" | wc -l)"
else
    echo "✗ Failed to download population panel"
    exit 1
fi
