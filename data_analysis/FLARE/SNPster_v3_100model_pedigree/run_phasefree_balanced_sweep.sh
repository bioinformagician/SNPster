#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_BALANCED_OUTDIR:-results_v3_phasefree_balanced_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_balanced_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_phasefree_balanced_sweep.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/$output_name"
name="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    mkdir "$output"
    cp "$project_dir/results_v3_phasefree_rich_01/image.txt" "$output/image.txt"
fi
args=(--init --network none --cpus=4 --memory=16g --user "$(id -u):$(id -g)"
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1
    -e OPENBLAS_NUM_THREADS=2 -e OMP_NUM_THREADS=2
    -v "$project_dir:/workspace:ro"
    -v "$project_dir/results_v3_phasefree_stack_01:/stack:ro"
    -v "$project_dir/results_v3_phasefree_calibration_01:/calibration:ro"
    -v "$project_dir/results_v3_phasefree_rich_01:/rich:ro"
    -v "$project_dir/results_v3_dosage_screen_01:/dosage:ro"
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro"
    -v "$project_dir/results_v3_phasefree_array_validation_01:/validation:ro"
    -v "$output:/sweep:rw" -w /tmp)
image="$(cat "$output/image.txt")"
if [[ "$1" == --run ]]; then
    if [[ ! -f "$output/smoke_passed.json" || ! -f "$output/augmentation_verification.json" ]]; then
        echo "Verify the augmentation bank and complete the five-model smoke fold first." >&2
        exit 1
    fi
    if docker inspect "$name" >/dev/null 2>&1; then
        echo "Inspect the retained balanced-sweep container before explicitly restarting it." >&2
        exit 1
    fi
    docker run --detach --name "$name" --log-driver local --log-opt max-size=10m --log-opt max-file=3 \
        "${args[@]}" "$image" python3 /workspace/pedigree_benchmark/scripts/sweep_phasefree_balanced.py \
        --stage run > "$output/container.id"
    docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
else
    if [[ -n "$(docker ps --filter "name=${name}" --format '{{.Names}}')" ]]; then
        echo "A balanced-mixture sweep is already running." >&2
        exit 1
    fi
    docker run --rm --name "${name}-${1#--}" "${args[@]}" "$image" \
        python3 /workspace/pedigree_benchmark/scripts/sweep_phasefree_balanced.py --stage "${1#--}" \
        2>&1 | tee "$output/${1#--}.log"
fi
