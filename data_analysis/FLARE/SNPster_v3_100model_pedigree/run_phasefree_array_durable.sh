#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${PHASEFREE_VALIDATION_OUTDIR:-results_v3_phasefree_array_validation_01}"
if [[ ! "$output_name" =~ ^results_v3_phasefree_array_validation_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --run && "$1" != --recover) ]]; then
    echo "Usage: bash run_phasefree_array_durable.sh --run|--recover" >&2
    exit 1
fi
output="$project_dir/$output_name"
prefix="snpster-${output_name}"
name="${prefix}-durable-01"
if [[ -n "$(docker ps --filter "name=${prefix}" --format '{{.Names}}')" ]]; then
    echo "A validation container is already running; inspect it before launching another." >&2
    exit 1
fi
if docker inspect "$name" >/dev/null 2>&1; then
    echo "Container $name already exists; inspect its exit status and logs rather than overwrite it." >&2
    exit 1
fi
args=()
if [[ "$1" == --recover ]]; then
    for pid_file in scheduler.pid verification_scheduler.pid; do
        if [[ -f "$output/$pid_file" ]]; then
            pid="$(cat "$output/$pid_file")"
            if [[ ! "$pid" =~ ^[0-9]+$ ]]; then echo "Invalid PID in $pid_file" >&2; exit 1; fi
            if kill -0 "$pid" 2>/dev/null; then
                echo "Recorded scheduler PID $pid still exists; inspect it before recovery." >&2
                exit 1
            fi
        fi
    done
    if [[ -d "$output/scheduler.lock" ]]; then
        archive="$output/scheduler.lock.interrupted-$(date +%Y%m%dT%H%M%S)"
        mv -T "$output/scheduler.lock" "$archive"
    fi
    args+=(--recover-preparation)
elif [[ -e "$output/scheduler.lock" ]]; then
    echo "Existing scheduler lock; inspect ownership and use --recover only after all workers stop." >&2
    exit 1
fi
image_id="$(cat "$output/image.txt")"
docker run --detach --init --name "$name" --network none \
    --cpus=12 --memory=64g --user "$(id -u):$(id -g)" \
    --log-driver local --log-opt max-size=10m --log-opt max-file=3 \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
    -v "$output:/validation:rw" \
    -v /srv/dependencies/ancestry:/data:ro \
    -v /srv/dependencies/imputation_runner/imputer:/imputer:ro \
    -w /tmp "$image_id" \
    python3 /workspace/pedigree_benchmark/scripts/schedule_phasefree_arrays.py \
    --output /validation "${args[@]}" > "$output/durable_container.id"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
