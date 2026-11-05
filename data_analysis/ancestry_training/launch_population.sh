#!/usr/bin/env bash
set -euo pipefail
[[ $# == 2 || ( $# == 3 && "$3" == --resume ) ]] || {
    echo "Usage: launch_population.sh specification.json new_run [--resume]" >&2; exit 2;
}
source_dir="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$source_dir" rev-parse --show-toplevel)"
specification="$(realpath "$1")"
root="$(realpath -m "$2")"
resume=false
[[ $# == 3 ]] && resume=true
[[ "$root" != / && "$root" != "$HOME" && "$root" != "$repo" && "$root" != "$source_dir" &&
   "$root" =~ ^/[A-Za-z0-9_./-]+$ && "$repo" =~ ^/[A-Za-z0-9_./-]+$ ]] || {
    echo "Use a dedicated run directory with safe container mount paths." >&2; exit 2;
}
cpus="${SNPSTER_CPUS:-22}"
memory="${SNPSTER_MEMORY_GB:-160}"
chunk_size="${SNPSTER_INFERENCE_CHUNK_SIZE:-1}"
[[ "$cpus" =~ ^[1-9][0-9]*$ && "$memory" =~ ^[1-9][0-9]*$ ]] && (( memory >= 6 )) || {
    echo "CPU budget must be positive and memory at least 6 GiB." >&2; exit 2;
}
[[ "$chunk_size" =~ ^[1-9][0-9]*$ ]] || {
    echo "SNPSTER_INFERENCE_CHUNK_SIZE must be a positive integer." >&2; exit 2;
}
check_source="$source_dir"
check_spec="$specification"
if [[ "$resume" == true ]]; then
    [[ -f "$root/population-followup" && -d "$root/planning/launch/.nextflow" &&
       ! -f "$root/reports/scoped-cache-cleanup.json" ]] || {
        echo "Population run is completed/cleaned or lacks retained planning history." >&2; exit 2;
    }
    [[ "$(cat "$root/population-followup")" == true && "$(cat "$root/java-profile")" == c1 ]] || {
        echo "Resume requires this exact dedicated C1 population experiment." >&2; exit 2;
    }
    (cd "$root"; sha256sum --check runtime.sha256)
    cmp --silent "$specification" "$root/experiment.json" || {
        echo "Do not change the frozen population protocol on resume." >&2; exit 2;
    }
    check_source="$root/source"
    check_spec="$root/experiment.json"
else
    [[ ! -e "$root" ]] || { echo "Use --resume for this exact existing population run." >&2; exit 2; }
fi
# Refuse early without creating a run; never inspect predecessor test scores.
docker run --rm --network none --cpus=1 --memory=6g --user "$(id -u):$(id -g)" \
    -e OPENBLAS_NUM_THREADS=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -v "$repo:$repo:ro" -v "$check_source:$check_source:ro" -v "$check_spec:$check_spec:ro" \
    -w "$check_source" --entrypoint python3 \
    sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
    population_followup.py check-ready --specification "$check_spec"
if [[ "$resume" == false ]]; then
    mkdir -p "$root"/{source,runtime,controller_home/framework,reports}
    cp "$source_dir/"*.py "$source_dir/"*.nf "$source_dir/nextflow.config" "$root/source/"
    cp "$source_dir/java_c1.config" "$root/runtime/java-c1-pinned.config"
    (cd "$root"; sha256sum runtime/java-c1-pinned.config) > "$root/runtime.sha256"
    cp "$specification" "$root/experiment.json"
    printf 'c1\n' > "$root/java-profile"
    printf 'true\n' > "$root/population-followup"
    (cd "$root/source"; sha256sum *.py *.nf nextflow.config) > "$root/source.sha256"
fi
(cd "$root/source"; sha256sum --check "$root/source.sha256")
[[ -z "$(docker ps -q --filter "label=snpster.population.root=$root")" ]] || {
    echo "A controller for this population run is already active." >&2; exit 2;
}
suffix="$(date +%Y%m%d_%H%M%S_%N)"
name="snpster-population-followup-$suffix"
printf '%s\n' "$name" > "$root/reports/$suffix.controller.name"
cp "$source_dir/launch_population.sh" "$root/reports/$suffix.launcher.sh"
options="--user $(id -u):$(id -g) -v $repo:$repo:ro -v $root/source:$root/source:ro -v /srv/dependencies:/srv/dependencies:ro"
recovery="${SNPSTER_NATIVE_RECOVERY:-}"
if [[ "$resume" == true && -f "$root/runtime/native-recovery.path" ]]; then
    saved_recovery="$(cat "$root/runtime/native-recovery.path")"
    [[ -z "$recovery" || "$recovery" == "$saved_recovery" ]] || {
        echo "Do not change the persisted native recovery manifest on resume." >&2; exit 2;
    }
    recovery="$saved_recovery"
    (cd "$root"; sha256sum --check runtime/native-recovery.sha256)
fi
[[ -z "$recovery" || ( "$resume" == true && -f "$recovery" ) ]] || {
    echo "Native recovery requires --resume and an existing verified manifest." >&2; exit 2;
}
if [[ -n "$recovery" && ! -f "$root/runtime/native-recovery.path" ]]; then
    recovery="$(realpath "$recovery")"
    [[ "$recovery" == "$root/reports/"* ]] || {
        echo "Native recovery must be a durable snapshot in this exact run's reports." >&2; exit 2;
    }
    printf '%s\n' "$recovery" > "$root/runtime/native-recovery.path"
    sha256sum "$recovery" > "$root/runtime/native-recovery.sha256"
fi
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Inference scheduling: $chunk_size person(s) per chunk; up to $cpus single-CPU tasks"
if docker run -t --name "$name" --label "snpster.population.root=$root" --init --network none --cpus=2 --memory=2g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c %g /var/run/docker.sock)" \
    -e NXF_OFFLINE=true -e NXF_VER=26.04.4 -e NXF_HOME="$root/controller_home" \
    -e ROOT="$root" -e SUFFIX="$suffix" -e RESUME="$resume" -e TASK_OPTIONS="$options" \
    -e BUDGET_CPUS="$cpus" -e BUDGET_MEMORY="$memory" \
    -e INFERENCE_CHUNK_SIZE="$chunk_size" \
    -e NATIVE_RECOVERY="$recovery" \
    -v "$repo:$repo:ro" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$root:$root:rw" -v "$root/source:$root/source:ro" \
    -v "$HOME/.nextflow/framework/26.04.4:$root/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$root" --entrypoint /bin/bash \
    sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 55
        docker run --rm --network none --cpus=1 --memory=256m --entrypoint /bin/true \
            sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
        echo "Task-container engine probe passed; population development comparison; no test or automatic expansion"
        runtime="$ROOT/source/nextflow.config,$ROOT/runtime/java-c1-pinned.config"
        extra=()
        [[ "$RESUME" == true ]] && extra=(-resume)
        nextflow_stage() {
            local stage=$1 workflow=$2
            shift 2
            mkdir -p "$ROOT/$stage"/{launch,work}
            cd "$ROOT/$stage/launch"
            echo "=== Population comparison: $stage ==="
            nextflow -log "$ROOT/reports/$SUFFIX.$stage.nextflow.log" -C "$runtime" run "$ROOT/source/$workflow" \
                "${extra[@]}" -pool-size "$BUDGET_CPUS" -ansi-log true -work-dir "$ROOT/$stage/work" \
                -with-trace "$ROOT/reports/$SUFFIX.$stage.trace.tsv" \
                -with-report "$ROOT/reports/$SUFFIX.$stage.report.html" \
                -with-timeline "$ROOT/reports/$SUFFIX.$stage.timeline.html" \
                -with-dag "$ROOT/reports/$SUFFIX.$stage.pipeline.html" \
                --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" --task_memory "6 GB" \
                --container_options "$TASK_OPTIONS" "$@"
        }
        python_stage() {
            docker run --rm --network none --cpus=1 --memory=6g $TASK_OPTIONS -v "$ROOT:$ROOT:rw" \
                -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -w "$ROOT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d "$@"
        }
        nextflow_stage planning population_followup.nf --population_stage plan --root "$ROOT" \
            --experiment "$ROOT/experiment.json" --output "$ROOT/planning"
        cat "$ROOT/planning/design/scope.json"
        native_manifest="${NATIVE_RECOVERY:-$ROOT/planning/design/native-reuse.json}"
        nextflow_stage screen guided_followup.nf --guided_stage run \
            --experiment "$ROOT/planning/design/config.json" --output "$ROOT/screen/results" \
            --phase_reuse "$ROOT/planning/design/phase-reuse.json" \
            --inference_chunk_size "$INFERENCE_CHUNK_SIZE" \
            --reuse_manifest "$native_manifest"
        python_stage audit.py "$ROOT/screen/results" | tee "$ROOT/reports/$SUFFIX.screen.independent-audit.log"
        nextflow_stage reporting population_followup.nf --population_stage report --root "$ROOT" --output "$ROOT/reporting"
        cat "$ROOT/reporting/report/population-comparison.json"
        python_stage population_followup.py cleanup --root "$ROOT" --output "$ROOT/reports/scoped-cache-cleanup.json"
        echo "POPULATION DEVELOPMENT COMPARISON AND INDEPENDENT AUDIT PASSED; CACHES CLEANED; NO TEST OR PROMOTION"
    '; then
    status=0
else
    status=$?
fi
docker inspect "$name" > "$root/reports/$suffix.controller.inspect.json"
exit "$status"
