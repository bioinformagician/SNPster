#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run) ]]; then
    echo "Usage: bash run_phasefree_allele_recovery.sh --prepare|--smoke|--run" >&2
    exit 1
fi
output="$project_dir/results_v3_phasefree_array_validation_01"
prefix="snpster-results_v3_phasefree_array_validation_01"
name="${prefix}-alleles-02"
if [[ -n "$(docker ps --filter "name=${prefix}" --format '{{.Names}}')" ]]; then
    echo "An array-validation container is already running." >&2
    exit 1
fi
args=(--init --network none --cpus=12 --memory=64g --user "$(id -u):$(id -g)"
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1
    -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1
    -v "$project_dir:/workspace:ro"
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro"
    -v "$output:/validation:rw" -v /srv/dependencies/ancestry:/data:ro
    -v /srv/dependencies/imputation_runner/imputer:/imputer:ro -w /tmp)
image="$(cat "$output/image.txt")"
if [[ "$1" == --run ]]; then
    if [[ ! -f "$output/allele_smoke_passed.json" ]]; then
        echo "Complete exact-allele smoke verification before durable recovery." >&2
        exit 1
    fi
    if docker inspect "$name" >/dev/null 2>&1; then
        echo "Inspect the existing retained recovery container before explicitly restarting it." >&2
        exit 1
    fi
    docker run --detach --name "$name" --log-driver local --log-opt max-size=10m --log-opt max-file=3 \
        "${args[@]}" "$image" python3 /workspace/pedigree_benchmark/scripts/recover_phasefree_alleles.py \
        --stage run > "$output/allele_container.id"
    docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
else
    docker run --rm --name "${name}-${1#--}" "${args[@]}" "$image" \
        python3 /workspace/pedigree_benchmark/scripts/recover_phasefree_alleles.py --stage "${1#--}" \
        2>&1 | tee "$output/allele_${1#--}.log"
fi
