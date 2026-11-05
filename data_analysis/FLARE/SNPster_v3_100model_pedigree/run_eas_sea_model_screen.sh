#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export EAS_SEA_SCREEN_OUTDIR="${EAS_SEA_SCREEN_OUTDIR:-results_v3_eas_sea_screen_01}"
output_name="$EAS_SEA_SCREEN_OUTDIR"
if [[ ! "$output_name" =~ ^results_v3_eas_sea_screen_[a-zA-Z0-9_]+$ ]]; then
    echo "EAS_SEA_SCREEN_OUTDIR must name an isolated results_v3_eas_sea_screen_* directory" >&2
    exit 1
fi
image="snpster-ancestry-benchmark:pedigree-100"

run_container() {
    local stage="$1" task="${2:-}" cpus=2 memory=8g name="$1"
    local args=(--stage "$stage" --outdir /results)
    if [[ "$stage" == prepare && -n "${EAS_SEA_REUSE_FROM:-}" ]]; then
        if [[ ! "$EAS_SEA_REUSE_FROM" =~ ^results_v3_eas_sea_screen_[a-zA-Z0-9_]+$ ]]; then
            echo "EAS_SEA_REUSE_FROM must name an existing isolated screen directory" >&2
            return 1
        fi
        args+=(--reuse-dir "/workspace/$EAS_SEA_REUSE_FROM")
    fi
    if [[ "$stage" == infer ]]; then
        cpus=6; memory=40g; name="infer-$task"
        args+=(--task "$task")
    fi
    docker run --rm --network none --name "snpster-${output_name}-${name}" \
        --cpus="$cpus" --memory="$memory" --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e SKLEARN_N_JOBS=2 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image" \
        python3 /workspace/pedigree_benchmark/scripts/screen_eas_sea_models.py \
            "${args[@]}" \
        2>&1 | tee "$project_dir/$output_name/$name.log"
}

if [[ "${1:-}" == --infer-task ]]; then
    run_container infer "$2"
    exit 0
fi
if [[ $# != 0 ]]; then
    echo "Unexpected launcher argument" >&2
    exit 1
fi
if [[ -e "$project_dir/$output_name" ]]; then
    echo "Output already exists; choose a new EAS/SEA screen directory" >&2
    exit 1
fi
docker image inspect "$image" --format '{{.Id}}' >/dev/null
mkdir "$project_dir/$output_name"
docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/container_image.txt"
run_container prepare
tail -n +2 "$project_dir/$output_name/pending_tasks.tsv" | cut -f1 |
    xargs -r -P3 -I{} bash "$project_dir/run_eas_sea_model_screen.sh" --infer-task '{}'
for stage in diag_k30_ridge diag_k30_extratrees diag_k36_ridge diag_k36_extratrees summarize; do
    run_container "$stage"
done
