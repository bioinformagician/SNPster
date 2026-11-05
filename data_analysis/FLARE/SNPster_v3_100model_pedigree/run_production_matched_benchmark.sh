#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cohort="${SNPSTER_MATCHED_COHORT:-$project/results_v3_observed_truth_support_20261022_114223/frozen/frozen/cohort.json}"
test -s "$cohort"
suffix="$(date +%Y%m%d_%H%M%S)"
run="${SNPSTER_MATCHED_RUN:-$project/results_v3_production_matched_$suffix}"
if [[ -e "$run" ]]; then
    echo "Experiment already exists; use resume_production_matched_benchmark.sh instead." >&2; exit 1
fi
mkdir -p "$run/smoke/launch" "$run/full/launch" "$run/controller_home"
sha256sum "$project/production_matched_benchmark.nf" "$project/production_matched.config" \
    "$project/config/production_matched_protocol.json" "$project/run_production_matched_benchmark.sh" \
    "$project/pedigree_benchmark/scripts/production_matched_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/prepare_production_truth_support.py" \
    "$project/pedigree_benchmark/scripts/production_matched_models.py" \
    "$project/pedigree_benchmark/scripts/verify_production_matched.py" \
    "$project/pedigree_benchmark/scripts/finalize_production_matched.py" \
    "$project/pedigree_benchmark/scripts/model_family_gnomix.py" \
    "$project/pedigree_benchmark/scripts/screen_model_families.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/results_v3_global_precision_01/precision.jar" "$cohort" > "$run/source_sha256.txt"
cp "$cohort" "$run/frozen_cohort.json"
SNPSTER_MATCHED_COHORT="$cohort" bash "$project/resume_production_matched_benchmark.sh" "$run"
