#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${MATCHED_GLOBAL_OUTDIR:-results_v3_matched_global_01}"
if [[ ! "$output_name" =~ ^results_v3_matched_global_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_matched_global_training.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/$output_name"
prefix="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    if [[ -e "$output" ]]; then echo "Preparation requires a fresh output" >&2; exit 1; fi
    mkdir "$output"
    docker image inspect snpster-ancestry-benchmark:matched-imputation \
        --format '{{.Id}}' > "$output/image.txt"
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
stage() {
    local action="$1" chrom="${2:-}"
    local args=(--stage "$action")
    if [[ -n "$chrom" ]]; then args+=(--chromosome "$chrom"); fi
    docker run --rm --network none --name "${prefix}-${action}${chrom:+-$chrom}" \
        --cpus=6 --memory=40g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -e NUMBA_CACHE_DIR=/tmp/numba-cache -e MPLCONFIGDIR=/tmp/matplotlib \
        -e XDG_CACHE_HOME=/tmp/cache -e TMPDIR=/tmp \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
        -v "$project_dir/results_v3_matched_imputation_02:/previous:ro" \
        -v "$project_dir/results_v3_phase_robustness_01:/robust:ro" \
        -v "$output:/pilot:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -v /srv/dependencies/imputation_runner/imputer:/imputer:ro \
        -w /tmp "$image_id" \
        python3 /workspace/pedigree_benchmark/scripts/train_matched_global.py "${args[@]}" \
        2>&1 | tee -a "$output/${action}${chrom:+_chr$chrom}.log"
}
chromosome() {
    stage impute "$1"
    stage features "$1"
}
names() {
    active+=("${prefix}-impute-$1" "${prefix}-features-$1")
}
if [[ "$1" == --prepare ]]; then
    active=("${prefix}-prepare")
    stage prepare
elif [[ "$1" == --smoke ]]; then
    names 22
    chromosome 22
else
    if [[ ! -f "$output/features/chr22/audit.json" ]]; then
        echo "Complete chromosome22 preflight with --smoke first" >&2
        exit 1
    fi
    for ((chrom=1; chrom<=22; chrom+=2)); do
        active=()
        names "$chrom"; names "$((chrom+1))"
        chromosome "$chrom" & first="$!"
        chromosome "$((chrom+1))" & second="$!"
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
