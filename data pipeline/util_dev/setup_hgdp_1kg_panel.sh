#!/bin/bash

set -euo pipefail

# Setup HGDP+1KG metadata and panel file
# This script downloads the metadata and converts it to panel format

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"
METADATA_DIR="/srv/dependencies/imputation_runner/imputer/beagle_references_hgdp"
PANEL_DIR="/srv/dependencies/imputation_runner"
TEMP_METADATA="/tmp/gnomad_meta_updated.tsv"
TEMP_PANEL="/tmp/hgdp_1kg_panel.txt"

echo "=================================================="
echo "HGDP+1KG Metadata and Panel Setup"
echo "=================================================="
echo ""

# Check if Python script exists
if [ ! -f "${SCRIPT_DIR}/convert_hgdp_1kg_metadata_to_panel.py" ]; then
    echo "Error: convert_hgdp_1kg_metadata_to_panel.py not found!"
    echo "Expected location: ${SCRIPT_DIR}/convert_hgdp_1kg_metadata_to_panel.py"
    exit 1
fi

# Download metadata
echo "Step 1: Downloading metadata file..."
echo "  URL: ${METADATA_URL}"
echo "  Temp location: ${TEMP_METADATA}"

if [ -f "${TEMP_METADATA}" ]; then
    echo "  Metadata already exists in /tmp/, skipping download"
else
    curl -L --retry 5 --retry-all-errors --progress-bar \
        -o "${TEMP_METADATA}" \
        "${METADATA_URL}"
    echo "  ✓ Download complete"
fi

echo ""

# Convert metadata to panel format
echo "Step 2: Converting metadata to panel format..."

# Activate virtual environment if it exists
if [ -d "/home/frederik/github_projects/SNPster/.venv-1" ]; then
    echo "  Activating virtual environment..."
    source /home/frederik/github_projects/SNPster/.venv-1/bin/activate
fi

# Check if pandas is available
if ! python3 -c "import pandas" 2>/dev/null; then
    echo "Error: pandas is not installed!"
    echo "Install it with: pip install pandas"
    exit 1
fi

python3 "${SCRIPT_DIR}/convert_hgdp_1kg_metadata_to_panel.py" \
    "${TEMP_METADATA}" \
    "${TEMP_PANEL}"

echo ""
echo "  ✓ Panel file created: ${TEMP_PANEL}"
echo ""

# Show preview
echo "Preview of panel file:"
head -5 "${TEMP_PANEL}"
echo ""

# Install files (requires sudo)
echo "Step 3: Installing files to system directories..."
echo ""
echo "The following files need to be copied (requires sudo):"
echo "  1. ${TEMP_METADATA}"
echo "     → ${METADATA_DIR}/gnomad_meta_updated.tsv"
echo "  2. ${TEMP_PANEL}"
echo "     → ${PANEL_DIR}/hgdp_1kg_panel.txt"
echo ""

read -p "Copy files now? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    sudo mkdir -p "${METADATA_DIR}"
    sudo mkdir -p "${PANEL_DIR}"
    
    sudo cp "${TEMP_METADATA}" "${METADATA_DIR}/gnomad_meta_updated.tsv"
    sudo cp "${TEMP_PANEL}" "${PANEL_DIR}/hgdp_1kg_panel.txt"
    
    sudo chmod 644 "${METADATA_DIR}/gnomad_meta_updated.tsv"
    sudo chmod 644 "${PANEL_DIR}/hgdp_1kg_panel.txt"
    
    echo ""
    echo "✓ Files installed successfully!"
    echo ""
    echo "Installed files:"
    ls -lh "${METADATA_DIR}/gnomad_meta_updated.tsv"
    ls -lh "${PANEL_DIR}/hgdp_1kg_panel.txt"
else
    echo ""
    echo "Files not copied. To copy manually, run:"
    echo "  sudo cp ${TEMP_METADATA} ${METADATA_DIR}/"
    echo "  sudo cp ${TEMP_PANEL} ${PANEL_DIR}/"
fi

echo ""
echo "=================================================="
echo "Setup complete!"
echo "=================================================="
echo ""
echo "Next steps:"
echo "  1. Update ancestry_module/config.py to use the new panel:"
echo "     POPULATION_PANEL_FILE = '${PANEL_DIR}/hgdp_1kg_panel.txt'"
echo "     REFERENCE_PANEL = 'HGDP+1KG'"
echo ""
echo "  2. Ensure reference VCF files are prepared:"
echo "     nextflow run prepare_hgdp_1kg_reference_panel.nf"
echo ""
