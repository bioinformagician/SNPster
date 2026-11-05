#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${FINE_PROFILE_OUTDIR:-results_v3_fine_profiles_01}"
if [[ ! "$output_name" =~ ^results_v3_fine_profiles_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_fine_profile_screen.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/$output_name"
if [[ "$1" == --prepare ]]; then
    if [[ -e "$output" ]]; then echo "Preparation requires a fresh output" >&2; exit 1; fi
    mkdir "$output"
    cp "$project_dir/results_v3_dosage_screen_01/image.txt" "$output/image.txt"
fi
image_id="$(cat "$output/image.txt")"
mkdir "$output/scheduler.lock"
action="${1#--}"
name="snpster-${output_name}-${action}"
cleanup() {
    local status="$?"
    trap - EXIT INT TERM
    if [[ "$status" != 0 ]]; then
        if [[ "$(docker inspect --format '{{.State.Running}}' "$name" 2>/dev/null || true)" == true ]]; then
            docker stop --timeout 10 "$name" >/dev/null
        fi
    fi
    rmdir "$output/scheduler.lock"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
if [[ "$action" == run && ! -f "$output/smoke.json" ]]; then
    echo "Complete --smoke before the full fine-profile comparison" >&2
    exit 1
fi
docker run --rm --network none --name "$name" \
    --cpus=6 --memory=30g --user "$(id -u):$(id -g)" \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OPENBLAS_NUM_THREADS=6 -e OMP_NUM_THREADS=6 -e MPLCONFIGDIR=/tmp/matplotlib \
    -e XDG_CACHE_HOME=/tmp/cache -e TMPDIR=/tmp \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
    -v "$project_dir/results_v3_dosage_screen_01:/dosage:ro" \
    -v "$output:/fine:rw" \
    -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
    -w /tmp "$image_id" \
    python3 /workspace/pedigree_benchmark/scripts/screen_fine_profiles.py --stage "$action" \
    2>&1 | tee -a "$output/${action}.log"
