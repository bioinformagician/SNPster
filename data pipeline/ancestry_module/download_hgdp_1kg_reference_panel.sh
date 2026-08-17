#!/usr/bin/env bash

set -euo pipefail

# Downloads the autosomal HGDP+1kGP phased GRCh38 reference panel and its indexes.
# Source: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/

REFERENCE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
METADATA_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"
OUTPUT_DIR="${1:-/srv/dependencies/imputation_runner/imputer/beagle_references_hgdp}"

if ! command -v curl >/dev/null 2>&1; then
    echo "curl is required but was not found on PATH." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    if [[ -f "${destination}" ]]; then
        echo "Already downloaded: $(basename "${destination}")"
        return
    fi

    echo "Downloading: $(basename "${destination}")"
    curl \
        --fail \
        --location \
        --retry 5 \
        --retry-all-errors \
        --continue-at - \
        --output "${partial}" \
        "${url}"
    mv "${partial}" "${destination}"
}

for chromosome in {1..22}; do
    filename="hgdp1kgp_chr${chromosome}.filtered.SNV_INDEL.phased.shapeit5.bcf"
    download_file "${REFERENCE_URL}/${filename}" "${OUTPUT_DIR}/${filename}"
    download_file "${REFERENCE_URL}/${filename}.csi" "${OUTPUT_DIR}/${filename}.csi"
done

download_file "${METADATA_URL}" "${OUTPUT_DIR}/gnomad_meta_updated.tsv"

echo "HGDP+1kGP reference panel download complete: ${OUTPUT_DIR}"