#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

MAX_MODELS="${MAX_MODELS:-0}"
MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-0}"
SEARCH_SPACE="${SEARCH_SPACE:-config/search_space.json}"

mkdir -p config/generated/model_configs results_v3
rm -f config/generated/panel_configs.tsv config/generated/model_manifest.tsv config/generated/model_configs.json
rm -f config/generated/model_configs/*.json

python3 scripts/generate_model_configs.py \
  --search-space "$SEARCH_SPACE" \
  --outdir config/generated \
  --max-models "$MAX_MODELS"

python3 scripts/preflight_v3.py \
  --params config/params.json \
  --search-space "$SEARCH_SPACE"

# Changing config or helper scripts should invalidate downstream tasks,
# while identical FLARE artifacts remain cacheable under -resume.
CONFIG_FINGERPRINT="$(bash scripts/config_fingerprint.sh "$SEARCH_SPACE")"

echo "Search space:       ${SEARCH_SPACE}"
echo "Config fingerprint: ${CONFIG_FINGERPRINT}"

nextflow run main.nf \
  -resume \
  -params-file config/params.json \
  --search_space "/workspace/${SEARCH_SPACE}" \
  --generated_model_configs '/workspace/config/generated/model_configs/*.json' \
  --outdir /workspace/results_v3 \
  --config_fingerprint "$CONFIG_FINGERPRINT" \
  --max_stage2_candidates "$MAX_STAGE2_CANDIDATES"
