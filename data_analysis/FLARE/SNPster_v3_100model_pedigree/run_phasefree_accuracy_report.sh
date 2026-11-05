#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_ACCURACY_OUTDIR:-results_v3_phasefree_accuracy_targets_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_accuracy_targets_[a-zA-Z0-9_]+$ || $# -ne 0 ]]; then
    echo "Usage: bash run_phasefree_accuracy_report.sh" >&2
    exit 1
fi
output="$project_dir/$output_name"
mkdir "$output"
image_id="$(cat "$project_dir/results_v3_phasefree_array_validation_01/image.txt")"
docker run --rm --network none --cpus 1 --memory 3g --user "$(id -u):$(id -g)" \
    -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_phasefree_stack_01:/stack:ro" \
    -v "$project_dir/results_v3_phasefree_array_validation_01:/validation:ro" \
    -v "$output:/accuracy:rw" -w /tmp "$image_id" \
    python3 /workspace/pedigree_benchmark/scripts/report_phasefree_accuracy.py
