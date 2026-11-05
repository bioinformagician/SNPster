#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_VALIDATION_OUTDIR:-results_v3_phasefree_array_validation_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_array_validation_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run && "$1" != --resume) ]]; then
    echo "Usage: bash run_phasefree_array_validation.sh --prepare|--smoke|--run|--resume" >&2
    exit 1
fi
output="$project_dir/$output_name"
prefix="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    mkdir "$output"
    printf '%s\n' sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d > "$output/image.txt"
fi
image_id="$(cat "$output/image.txt")"
mkdir "$output/scheduler.lock"
active=()
cleanup() {
    local status="$?"
    trap - EXIT INT TERM
    if [[ "$status" != 0 ]]; then
        for name in "${active[@]}"; do
            if [[ "$(docker inspect --format '{{.State.Running}}' "$name" 2>/dev/null || true)" == true ]]; then
                docker stop --timeout 10 "$name" >/dev/null
            fi
        done
    fi
    rmdir "$output/scheduler.lock"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
retry=()
if [[ "$1" == --resume ]]; then retry=(--retry-incomplete); fi
stage() {
    local action="$1" chrom="${2:-}"
    local args=(--stage "$action")
    if [[ -n "$chrom" ]]; then args+=(--chromosome "$chrom" "${retry[@]}"); fi
    docker run --rm --network none --name "${prefix}-${action}${chrom:+-$chrom}" \
        --cpus=6 --memory=30g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
        -v "$output:/validation:rw" \
        -v /srv/dependencies/ancestry:/data:ro \
        -v /srv/dependencies/imputation_runner/imputer:/imputer:ro \
        -w /tmp "$image_id" \
        python3 /workspace/pedigree_benchmark/scripts/validate_phasefree_arrays.py \
        "${args[@]}" 2>&1 | tee -a "$output/${action}${chrom:+_chr$chrom}.log"
}
if [[ "$1" == --prepare ]]; then
    active=("${prefix}-prepare")
    stage prepare
elif [[ "$1" == --smoke ]]; then
    active=("${prefix}-chromosome-22")
    stage chromosome 22
else
    if [[ ! -f "$output/chr22/completion.json" ]]; then
        echo "Complete the chromosome22 --smoke preflight first" >&2
        exit 1
    fi
    for ((chrom=1; chrom<=22; chrom+=2)); do
        active=("${prefix}-chromosome-$chrom" "${prefix}-chromosome-$((chrom+1))")
        stage chromosome "$chrom" & first="$!"
        stage chromosome "$((chrom+1))" & second="$!"
        first_status=0; second_status=0
        wait "$first" || first_status="$?"
        wait "$second" || second_status="$?"
        if [[ "$first_status" != 0 || "$second_status" != 0 ]]; then
            echo "Chromosome jobs failed: $chrom=$first_status, $((chrom+1))=$second_status" >&2
            exit 1
        fi
    done
    active=("${prefix}-summarize")
    stage summarize
fi
