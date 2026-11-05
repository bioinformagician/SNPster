#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
MAX_MODELS="${MAX_MODELS:-0}"
MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-6000}"
SEARCH_SPACE="${SEARCH_SPACE:-config/search_space_100.json}"
PEDIGREE_DESIGN="${PEDIGREE_DESIGN:-reduced}"
RUN_OUTDIR="${RUN_OUTDIR:-/workspace/results_v3_sweep_100_disjoint_${PEDIGREE_DESIGN}}"
RECOVERY_STAGE2_MANIFEST="${RECOVERY_STAGE2_MANIFEST:-}"
RECOVERY_SOURCE_OUTDIR="${RECOVERY_SOURCE_OUTDIR:-$RUN_OUTDIR}"
RESUME_SESSION="${RESUME_SESSION:-}"
FROZEN_CANDIDATES_MANIFEST="${FROZEN_CANDIDATES_MANIFEST:-}"
RUN_WORKDIR="${RUN_WORKDIR:-/workspace/work}"
if [[ -n "$FROZEN_CANDIDATES_MANIFEST" && -n "$RECOVERY_STAGE2_MANIFEST" ]]; then
  echo 'Choose either frozen pilot mode or stage-2 recovery, not both' >&2
  exit 2
fi
if [[ -n "$FROZEN_CANDIDATES_MANIFEST" && "${PEDIGREE_ENABLED:-true}" != true ]]; then
  echo 'Frozen pilot requires pedigree validation' >&2
  exit 2
fi
if [[ "${PEDIGREE_ENABLED:-true}" == true ]]; then
  python3 pedigree_benchmark/scripts/preflight_founders.py \
    --params config/params.json --outdir "$RUN_OUTDIR" \
    --design "$PEDIGREE_DESIGN" --families "${PEDIGREE_FAMILIES:-90}" \
    --batch-size "${PEDIGREE_BATCH_SIZE:-10}" --seed "${PEDIGREE_SEED:-20261008}" \
    --allow-reuse "${PEDIGREE_ALLOW_REUSE:-true}"
fi
mkdir -p config/generated/model_configs "$RUN_OUTDIR" "$RUN_WORKDIR" .nextflow
if [[ -n "$FROZEN_CANDIDATES_MANIFEST" ]]; then
  [[ -s "$FROZEN_CANDIDATES_MANIFEST" ]] || { echo 'Missing frozen candidate manifest' >&2; exit 2; }
  echo 'Frozen pilot: skipping representation fitting, decoder fitting and reference inference'
elif [[ -z "$RECOVERY_STAGE2_MANIFEST" ]]; then
  bash scripts/generate_model_configs_cached.sh "$SEARCH_SPACE" config/generated "$MAX_MODELS"
  python3 scripts/preflight_v3.py --params config/params.json --search-space "$SEARCH_SPACE"
else
  [[ -s "$RECOVERY_STAGE2_MANIFEST" ]] || { echo "Missing recovery stage2 manifest: $RECOVERY_STAGE2_MANIFEST" >&2; exit 2; }
  [[ -s "$RECOVERY_SOURCE_OUTDIR/selection/selected_representations.tsv" ]] || { echo 'Missing selected_representations.tsv for recovery' >&2; exit 2; }
  echo 'Recovery mode: bypassing all upstream training and stage2 fitting'
fi
CONFIG_FINGERPRINT="$(bash scripts/config_fingerprint.sh "$SEARCH_SPACE")"
PDB_ARGS=(--pedigree_enabled "${PEDIGREE_ENABLED:-true}")
if [[ "${PEDIGREE_ENABLED:-true}" == true ]]; then
  PEDIGREE_FINGERPRINT="$(find pedigree_benchmark/scripts -name '*.py' -type f -print0 | sort -z | xargs -0 -r sha256sum | sha256sum | awk '{print $1}')"
  PDB_ARGS+=(--pedigree_fingerprint "$PEDIGREE_FINGERPRINT" \
    --unified_fingerprint "$(sha256sum config/unified_selection.json | awk '{print $1}')" \
    --pedigree_families "${PEDIGREE_FAMILIES:-90}" \
    --pedigree_batch_size "${PEDIGREE_BATCH_SIZE:-10}" \
    --pedigree_seed "${PEDIGREE_SEED:-20261008}" \
    --pedigree_design "$PEDIGREE_DESIGN" \
    --pedigree_allow_reuse "${PEDIGREE_ALLOW_REUSE:-true}" \
    --pedigree_keep_high_density "${PEDIGREE_KEEP_HIGH_DENSITY:-false}" \
    --pedigree_chromosomes "${PEDIGREE_CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}")
fi
[[ -z "$RECOVERY_STAGE2_MANIFEST" ]] || PDB_ARGS+=(--resume_stage2_manifest "$RECOVERY_STAGE2_MANIFEST" --resume_source_outdir "$RECOVERY_SOURCE_OUTDIR")
ENTRYPOINT=main.nf
REPORT_ARGS=()
if [[ -n "$FROZEN_CANDIDATES_MANIFEST" ]]; then
  ENTRYPOINT=pedigree_pilot.nf
  PDB_ARGS+=(--frozen_candidates_manifest "$FROZEN_CANDIDATES_MANIFEST")
  [[ -z "${PEDIGREE_REUSE_DIR:-}" ]] || PDB_ARGS+=(--pedigree_reuse_dir "$PEDIGREE_REUSE_DIR")
  REPORT_ARGS=(-with-trace "$RUN_OUTDIR/trace.txt" -with-report "$RUN_OUTDIR/report.html"
               -with-timeline "$RUN_OUTDIR/timeline.html")
fi
RESUME=(-resume)
[[ -z "$RESUME_SESSION" ]] || RESUME=(-resume "$RESUME_SESSION")
echo "Search: $SEARCH_SPACE | output: $RUN_OUTDIR | fingerprint: $CONFIG_FINGERPRINT"
nextflow run "$ENTRYPOINT" "${RESUME[@]}" "${REPORT_ARGS[@]}" -work-dir "$RUN_WORKDIR" -params-file config/params.json \
  --search_space "/workspace/${SEARCH_SPACE}" \
  --generated_model_configs '/workspace/config/generated/model_configs/*.json' \
  --outdir "$RUN_OUTDIR" \
  --config_fingerprint "$CONFIG_FINGERPRINT" \
  --max_stage2_candidates "$MAX_STAGE2_CANDIDATES" \
  "${PDB_ARGS[@]}"
