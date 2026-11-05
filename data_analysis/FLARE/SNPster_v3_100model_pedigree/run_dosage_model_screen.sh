#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${DOSAGE_SCREEN_OUTDIR:-results_v3_dosage_screen_01}"
if [[ ! "$output_name" =~ ^results_v3_dosage_screen_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_dosage_model_screen.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/$output_name"
prefix="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    if [[ -e "$output" ]]; then echo "Preparation requires a fresh output" >&2; exit 1; fi
    mkdir "$output"
    docker image inspect snpster-ancestry-benchmark:gnomix-screen \
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
        -e OPENBLAS_NUM_THREADS=6 -e OMP_NUM_THREADS=6 -e MPLCONFIGDIR=/tmp/matplotlib \
        -e XDG_CACHE_HOME=/tmp/cache -e TMPDIR=/tmp \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
        -v "$project_dir/results_v3_matched_global_01:/pilot:ro" \
        -v "$output:/dosage:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image_id" \
        python3 /workspace/pedigree_benchmark/scripts/screen_dosage_models.py "${args[@]}" \
        2>&1 | tee -a "$output/${action}${chrom:+_chr$chrom}.log"
}
if [[ "$1" == --prepare ]]; then
    active=("${prefix}-prepare")
    stage prepare
elif [[ "$1" == --smoke ]]; then
    active=("${prefix}-extract-22")
    stage extract 22
else
    if [[ ! -f "$output/chr22/audit.json" ]]; then
        echo "Complete chromosome22 compatibility check with --smoke first" >&2
        exit 1
    fi
    for ((chrom=1; chrom<=22; chrom+=2)); do
        active=("${prefix}-extract-$chrom" "${prefix}-extract-$((chrom+1))")
        stage extract "$chrom" & first="$!"
        stage extract "$((chrom+1))" & second="$!"
        first_status=0; second_status=0
        wait "$first" || first_status="$?"
        wait "$second" || second_status="$?"
        if [[ "$first_status" != 0 || "$second_status" != 0 ]]; then
            echo "Dosage extraction failed: $chrom=$first_status, $((chrom+1))=$second_status" >&2
            exit 1
        fi
    done
    active=("${prefix}-fit")
    stage fit
fi
