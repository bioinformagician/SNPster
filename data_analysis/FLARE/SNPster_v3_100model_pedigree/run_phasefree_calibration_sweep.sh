#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_CALIBRATION_OUTDIR:-results_v3_phasefree_calibration_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_calibration_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_phasefree_calibration_sweep.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/$output_name"
name="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    mkdir "$output"
    cp "$project_dir/results_v3_phasefree_array_validation_01/image.txt" "$output/image.txt"
fi
image_id="$(cat "$output/image.txt")"
args=(--init --network none --cpus=4 --memory=16g --user "$(id -u):$(id -g)"
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1
    -e OPENBLAS_NUM_THREADS=2 -e OMP_NUM_THREADS=2
    -v "$project_dir:/workspace:ro"
    -v "$project_dir/results_v3_phasefree_stack_01:/stack:ro"
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro"
    -v "$project_dir/results_v3_phasefree_array_validation_01:/validation:ro"
    -v "$output:/sweep:rw" -w /tmp)
stage="${1#--}"
if [[ "$1" == --run ]]; then
    if [[ ! -f "$output/smoke_passed.json" ]]; then
        echo "Complete the ten-model --smoke fold first." >&2
        exit 1
    fi
    if docker inspect "$name" >/dev/null 2>&1; then
        echo "Existing durable sweep container: inspect it; use docker start only after reviewing a stopped run." >&2
        exit 1
    fi
    docker run --detach --name "$name" --log-driver local --log-opt max-size=10m --log-opt max-file=3 \
        "${args[@]}" "$image_id" \
        python3 /workspace/pedigree_benchmark/scripts/sweep_phasefree_calibration.py --stage run \
        > "$output/container.id"
    docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
else
    if [[ -n "$(docker ps --filter "name=${name}" --format '{{.Names}}')" ]]; then
        echo "A sweep container is already active." >&2
        exit 1
    fi
    docker run --rm --name "${name}-${stage}" "${args[@]}" "$image_id" \
        python3 /workspace/pedigree_benchmark/scripts/sweep_phasefree_calibration.py --stage "$stage" \
        2>&1 | tee "$output/${stage}.log"
fi
