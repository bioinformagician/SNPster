#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_ERROR_AUDIT_OUTDIR:-results_v3_phasefree_error_audit_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_error_audit_[a-zA-Z0-9_]+$ || $# -ne 0 ]]; then
    echo "Usage: bash run_phasefree_error_audit.sh" >&2
    exit 1
fi
output="$project_dir/$output_name"
mkdir "$output"
image_id="$(cat "$project_dir/results_v3_phasefree_array_validation_01/image.txt")"
docker run --rm --network none --cpus 2 --memory 5g --user "$(id -u):$(id -g)" \
    -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=2 -e OMP_NUM_THREADS=2 \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_phasefree_stack_01:/stack:ro" \
    -v "$output:/error-audit:rw" -w /tmp "$image_id" \
    python3 /workspace/pedigree_benchmark/scripts/audit_phasefree_errors.py \
    --output /error-audit
