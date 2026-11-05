#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 || ( $# -eq 3 && "$3" != --resume ) ]]; then
    echo "Usage: bash run_observed.sh ABSOLUTE_VERIFIED_SOURCE ABSOLUTE_NEW_RUN [--resume]" >&2
    exit 2
fi
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$source_dir/../.." && pwd)"
source="$(realpath "$1")"
experiment="$(realpath -m "$2")"
cpu="${SNPSTER_CPUS:-22}"
memory="${SNPSTER_MEMORY_GB:-178}"
task_memory="${SNPSTER_TASK_MEMORY_GB:-8}"
for path in "$source" "$experiment" "$repository"; do
    [[ "$path" =~ ^/[A-Za-z0-9_./-]+$ && "$path" != / && "$path" != "$HOME" ]] || {
        echo "Invalid absolute mount path: $path" >&2; exit 2;
    }
done
[[ -f "$source/completion.json" && "$source" == "$repository/"* ]] || {
    echo "The verified observed source must be retained inside this repository." >&2; exit 2;
}
[[ "$experiment" != "$repository" && "$experiment" != "$source" && "$experiment" != "$source_dir" ]] || {
    echo "Choose a new dedicated experiment directory." >&2; exit 2;
}
[[ "$cpu" =~ ^[1-9][0-9]*$ && "$memory" =~ ^[1-9][0-9]*$ && "$task_memory" =~ ^[1-9][0-9]*$ ]] || {
    echo "CPU/GiB resource budgets must be positive integers." >&2; exit 2;
}
(( memory >= 16 && task_memory >= 2 && task_memory <= memory )) || {
    echo "Require at least 16 GiB total (PCA base), and valid per-task memory." >&2; exit 2;
}
resume=false
if [[ "${3:-}" == --resume ]]; then
    resume=true
    [[ -f "$experiment/source.sha256" ]] || { echo "Missing exact-run snapshot." >&2; exit 2; }
    [[ "$(cat "$experiment/input-source.path")" == "$source" ]] || {
        echo "Resume source differs from the frozen run." >&2; exit 2;
    }
else
    [[ ! -e "$experiment" ]] || { echo "Run directory exists; use --resume." >&2; exit 2; }
    mkdir -p "$experiment"/{source,reports,controller_home,prepare/launch,prepare/work,prepare/results,smoke/launch,smoke/work,smoke/results,full/launch,full/work,full/results}
    cp "$source_dir/"*.py "$source_dir/"*.nf "$source_dir/nextflow.config" "$experiment/source/"
    cp "$source_dir/post_sweep_v1.draft.json" "$experiment/template.json"
    printf '%s\n' "$source" > "$experiment/input-source.path"
    (cd "$experiment/source"; sha256sum *.py *.nf nextflow.config) > "$experiment/source.sha256"
fi
(cd "$experiment/source"; sha256sum --check "$experiment/source.sha256")
suffix="$(date +%Y%m%d_%H%M%S_%N)"
name="snpster-ancestry-observed-$suffix"
printf '%s\n' "$name" > "$experiment/reports/$suffix.controller.name"
cp "$source_dir/run_observed.sh" "$experiment/reports/$suffix.launcher.sh"
sha256sum "$experiment/template.json" "$source_dir/run_observed.sh" > "$experiment/reports/$suffix.launch-inputs.sha256"
image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v $experiment/prepare/results:$experiment/prepare/results:ro -v /srv/dependencies:/srv/dependencies:ro"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Fresh run: $experiment"
echo "Budget: $cpu one-core task slots, $memory GiB tasks plus 2 GiB controller."
echo "Stages: verified input preparation -> real chr22 smoke -> independent audit -> all-autosome run -> independent audit."
echo "No old model fits imported. Imputed routes disabled; independent_test=false."
set +e
docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e SOURCE="$source" \
    -e SUFFIX="$suffix" -e BUDGET_CPUS="$cpu" -e BUDGET_MEMORY="$memory" \
    -e TASK_MEMORY="$task_memory" -e CONTAINER_OPTIONS="$options" -e RESUME="$resume" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment" \
    --entrypoint /bin/bash "$image" -euo pipefail -c '
        stty cols 180 rows 55
        trap '\''status=$?; printf "%s\n" "$status" > "$EXPERIMENT/reports/$SUFFIX.pipeline.exit-status"; echo "Controller pipeline exited: $status (work/history retained)"; exit "$status"'\'' EXIT
        docker run --rm --network none --cpus=1 --memory=128m \
            --entrypoint /bin/true sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
        echo "Task-container engine probe passed."
        extra=()
        [[ "$RESUME" == true ]] && extra=(-resume)
        pipeline() {
            local stage=$1 workflow=$2 config=$3
            shift 3
            cd "$EXPERIMENT/$stage/launch"
            echo "=== $stage: $workflow ==="
            nextflow -log "$EXPERIMENT/reports/$SUFFIX.$stage.nextflow.log" \
                -C "$EXPERIMENT/source/nextflow.config" run "$EXPERIMENT/source/$workflow" \
                "${extra[@]}" -ansi-log true -work-dir "$EXPERIMENT/$stage/work" \
                -with-trace "$EXPERIMENT/reports/$SUFFIX.$stage.trace.tsv" \
                -with-report "$EXPERIMENT/reports/$SUFFIX.$stage.report.html" \
                -with-timeline "$EXPERIMENT/reports/$SUFFIX.$stage.timeline.html" \
                -with-dag "$EXPERIMENT/reports/$SUFFIX.$stage.pipeline.html" \
                --experiment "$config" --output "$EXPERIMENT/$stage/results" \
                --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" \
                --task_memory "$TASK_MEMORY GB" --pca_memory "16 GB" --container_options "$CONTAINER_OPTIONS" "$@"
        }
        audit() {
            local stage=$1
            echo "=== Independent durable-output audit: $stage ==="
            docker run --rm --network none --cpus=1 --memory="${TASK_MEMORY}g" \
                --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -v "$EXPERIMENT:$EXPERIMENT:rw" -v /srv/dependencies:/srv/dependencies:ro \
                -v "'"$repository:$repository"':ro" -w "$EXPERIMENT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                audit.py "$EXPERIMENT/$stage/results" 2>&1 | tee "$EXPERIMENT/reports/$SUFFIX.$stage.independent-audit.log"
        }
        pipeline prepare prepare_observed.nf "$EXPERIMENT/template.json" --source "$SOURCE"
        pipeline smoke main.nf "$EXPERIMENT/prepare/results/ready/smoke.json"
        audit smoke
        echo "Real-cohort chr22 smoke and independent audit passed; starting full run."
        pipeline full main.nf "$EXPERIMENT/prepare/results/ready/full.json"
        audit full
        echo "FULL WORKFLOW AND INDEPENDENT AUDIT PASSED."
        echo "Durable results: $EXPERIMENT/full/results"
    '
status=$?
set -e
docker logs "$name" > "$experiment/reports/$suffix.controller.log" 2>&1
docker inspect "$name" > "$experiment/reports/$suffix.controller.json"
printf '%s\n' "$status" > "$experiment/reports/$suffix.exit-status"
echo "Controller $name exited $status; reports and exact-run history retained."
exit "$status"
