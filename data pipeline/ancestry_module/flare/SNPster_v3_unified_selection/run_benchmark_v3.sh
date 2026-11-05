#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

MAX_MODELS="${MAX_MODELS:-0}"
MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-0}"
SEARCH_SPACE="${SEARCH_SPACE:-config/search_space.json}"

mkdir -p config/generated/model_configs results_v3
bash scripts/generate_model_configs_cached.sh "$SEARCH_SPACE" config/generated "$MAX_MODELS"

python3 scripts/preflight_v3.py \
  --params config/params.json \
  --search-space "$SEARCH_SPACE"

# Changing config or helper scripts should invalidate downstream tasks,
# while identical FLARE artifacts remain cacheable under -resume.
CONFIG_FINGERPRINT="$(bash scripts/config_fingerprint.sh "$SEARCH_SPACE")"

echo "Search space:       ${SEARCH_SPACE}"
echo "Config fingerprint: ${CONFIG_FINGERPRINT}"

PEDIGREE_ARGS=(--pedigree_enabled "${PEDIGREE_ENABLED:-true}")
if [[ "${PEDIGREE_ENABLED:-true}" == "true" ]]; then
  PEDIGREE_FINGERPRINT="$(find pedigree_benchmark/scripts -type f -name '*.py' -print0 | sort -z | xargs -0 -r sha256sum | sha256sum | cut -d' ' -f1)"
  PEDIGREE_ARGS+=(
    --pedigree_fingerprint "$PEDIGREE_FINGERPRINT"
    --unified_fingerprint "$(sha256sum config/unified_selection.json | cut -d' ' -f1)"
    --pedigree_families "${PEDIGREE_FAMILIES:-90}"
    --pedigree_batch_size "${PEDIGREE_BATCH_SIZE:-10}"
    --pedigree_seed "${PEDIGREE_SEED:-20261008}"
    --pedigree_allow_reuse "${PEDIGREE_ALLOW_REUSE:-true}"
    --pedigree_keep_high_density "${PEDIGREE_KEEP_HIGH_DENSITY:-false}"
    --pedigree_chromosomes "${PEDIGREE_CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
  )
  echo "Unified ancestry selection: enabled; ${PEDIGREE_FAMILIES:-90} families"
fi

nextflow run main.nf \
  -resume \
  -params-file config/params.json \
  --search_space "/workspace/${SEARCH_SPACE}" \
  --generated_model_configs '/workspace/config/generated/model_configs/*.json' \
  --outdir /workspace/results_v3 \
  --config_fingerprint "$CONFIG_FINGERPRINT" \
  --max_stage2_candidates "$MAX_STAGE2_CANDIDATES" \
  "${PEDIGREE_ARGS[@]}"
