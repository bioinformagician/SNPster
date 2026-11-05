#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${FLARE_NATIVE_OUTDIR:-results_v3_flare_native_01}"
if [[ ! "$output_name" =~ ^results_v3_flare_native_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --run) ]]; then
    echo "Usage: bash run_flare_native_diagnostic.sh --prepare|--run" >&2
    exit 1
fi
output="$project/$output_name"
name="snpster-${output_name}"
if [[ "$1" == --prepare ]]; then
    mkdir "$output"
    cp "$project/results_v3_observed_dense_01/image.txt" "$output/image.txt"
fi
args=(--init --network none --cpus=12 --memory=32g --user "$(id -u):$(id -g)"
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1
    -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1
    -v "$project:/workspace:ro"
    -v "$project/results_v3_observed_dense_01:/benchmark:ro"
    -v "$project/results_v3_model_family_screen_05:/results:ro"
    -v "$project/results_v3_phasefree_array_validation_01:/legacy:ro"
    -v "$project/results_v3_phasefree_array_validation_01:/validation:ro"
    -v /srv/dependencies/ancestry:/data:ro
    -v "$output:/diagnostic:rw" -w /tmp)
image="$(cat "$output/image.txt")"
if [[ "$1" == --run ]]; then
    if [[ ! -f "$output/experiment.json" ]]; then
        echo "Prepare immutable inputs and training-only selection first." >&2
        exit 1
    fi
    if docker inspect "$name" >/dev/null 2>&1; then
        echo "Inspect existing diagnostic container before explicit recovery." >&2
        exit 1
    fi
    docker run --detach --name "$name" --log-driver local --log-opt max-size=10m --log-opt max-file=3 \
        "${args[@]}" "$image" python3 /workspace/pedigree_benchmark/scripts/diagnose_flare_native.py \
        --stage run > "$output/container.id"
    docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
else
    docker run --rm --name "${name}-prepare" "${args[@]}" "$image" \
        python3 /workspace/pedigree_benchmark/scripts/diagnose_flare_native.py --stage prepare \
        2>&1 | tee "$output/prepare.log"
fi
