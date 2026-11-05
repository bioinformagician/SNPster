#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_READINESS_OUTDIR:-results_v3_phasefree_readiness_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_readiness_[a-zA-Z0-9_]+$ || $# -ne 0 ]]; then
    echo "Usage: bash run_phasefree_readiness_audit.sh" >&2
    exit 1
fi
output="$project_dir/$output_name"
mkdir "$output"
cp "$project_dir/results_v3_dosage_screen_01/image.txt" "$output/image.txt"
image_id="$(cat "$output/image.txt")"
repo_dir="$(git -C "$project_dir" rev-parse --show-toplevel)"
name="snpster-${output_name}"
cleanup() {
    local status="$?"
    trap - EXIT INT TERM
    if [[ "$status" != 0 ]]; then
        if [[ "$(docker inspect --format '{{.State.Running}}' "$name" 2>/dev/null || true)" == true ]]; then
            docker stop --timeout 10 "$name" >/dev/null
        fi
    fi
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
docker run --rm --network none --name "$name" \
    --cpus=2 --memory=8g --user "$(id -u):$(id -g)" \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OPENBLAS_NUM_THREADS=2 -e OMP_NUM_THREADS=2 \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
    -v "$project_dir/results_v3_matched_global_01:/pilot:ro" \
    -v "$output:/readiness:rw" \
    -v "$repo_dir/data pipeline/qc_module:/production-qc:ro" \
    -v /srv/dependencies/ancestry:/data:ro \
    -v /srv/dependencies/imputation_runner/imputer:/imputer:ro \
    -w /tmp "$image_id" \
    python3 /workspace/pedigree_benchmark/scripts/audit_phasefree_readiness.py \
    --screen /results --matched /pilot \
    --source-metadata /data/flare2_test/source_model \
    --dense-vcf /imputer/beagle_references/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz \
    --qc-module /production-qc --output /readiness \
    2>&1 | tee "$output/audit.log"
