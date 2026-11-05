#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: bash launch.sh ABSOLUTE_CONFIG ABSOLUTE_RUN_DIR [--smoke-first|--staged-development|--guided-followup] [--resume]" >&2
    exit 2
fi
resume=false
smoke_first=false
staged=false
guided=false
for argument in "${@:3}"; do
    case "$argument" in
        --resume) resume=true ;;
        --smoke-first) smoke_first=true ;;
        --staged-development) staged=true ;;
        --guided-followup) guided=true ;;
        *) echo "Unknown launch option: $argument" >&2; exit 2 ;;
    esac
done
[[ "$staged" == false || "$smoke_first" == false ]] || {
    echo "Staged development replaces, rather than extends, the old smoke/full sweep." >&2; exit 2;
}
[[ "$guided" == false || ( "$staged" == false && "$smoke_first" == false ) ]] || {
    echo "Guided follow-up cannot be combined with another expansion policy." >&2; exit 2;
}
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$source_dir/../.." && pwd)"
config="$(realpath "$1")"
experiment="$(realpath -m "$2")"
[[ -f "$config" ]] || { echo "Config file does not exist: $config" >&2; exit 2; }
[[ "$experiment" == /* && "$experiment" != "$repository" && "$experiment" != / && "$experiment" != "$HOME" ]] || {
    echo "Use a new dedicated run directory, never a repository/home/filesystem root." >&2
    exit 2
}
[[ "$repository" =~ ^/[A-Za-z0-9_./-]+$ && "$experiment" =~ ^/[A-Za-z0-9_./-]+$ ]] || {
    echo "Container mount paths must contain only letters, digits, dots, slashes, underscores and hyphens." >&2
    exit 2
}
cpu="${SNPSTER_CPUS:-1}"
budget="${SNPSTER_MEMORY_GB:-8}"
task_memory="${SNPSTER_TASK_MEMORY_GB:-8}"
pca_memory="${SNPSTER_PCA_MEMORY_GB:-$task_memory}"
[[ "$cpu" =~ ^[1-9][0-9]*$ && "$budget" =~ ^[1-9][0-9]*$ && "$task_memory" =~ ^[1-9][0-9]*$ &&
   "$pca_memory" =~ ^[1-9][0-9]*$ ]] || {
    echo "CPU and GiB memory budgets must be positive integers." >&2
    exit 2
}
(( task_memory <= budget && pca_memory <= budget )) || { echo "Task memory exceeds the whole task budget." >&2; exit 2; }
previous_java_profile=default
[[ -f "$experiment/java-profile" ]] && previous_java_profile="$(cat "$experiment/java-profile")"
java_profile="${SNPSTER_JAVA_PROFILE:-$previous_java_profile}"
[[ "$java_profile" == default || "$java_profile" == c1 ]] || {
    echo "SNPSTER_JAVA_PROFILE must be default or c1." >&2; exit 2;
}
[[ "$previous_java_profile" != c1 || "$java_profile" == c1 ]] || {
    echo "Do not mix an existing C1 recovery with default-JIT caches; use a new run to change back." >&2; exit 2;
}
[[ "$staged" == false || "$java_profile" == c1 ]] || {
    echo "Staged recovery requires SNPSTER_JAVA_PROFILE=c1 to match the retained native evidence." >&2; exit 2;
}
[[ "$guided" == false || "$java_profile" == c1 ]] || {
    echo "Guided follow-up requires SNPSTER_JAVA_PROFILE=c1 to match the verified phasing runtime." >&2; exit 2;
}
if [[ "$resume" == true ]]; then
    if [[ ( "$staged" == true || "$guided" == true ) && -f "$experiment/reports/scoped-cache-cleanup.json" ]]; then
        echo "This staged comparison already passed both audits and its work caches were removed; do not rerun it with --resume." >&2
        exit 2
    fi
    history="$experiment/launch/.nextflow"
    [[ "$smoke_first" == true ]] && history="$experiment/smoke/launch/.nextflow"
    [[ "$staged" == true ]] && history="$experiment/planning/launch/.nextflow"
    [[ "$guided" == true ]] && history="$experiment/planning/launch/.nextflow"
    [[ -d "$history" ]] || {
        echo "No retained Nextflow history at $history" >&2; exit 2;
    }
    cmp --silent "$config" "$experiment/experiment.json" || {
        echo "Config changed; create a new run directory, not a misleading resume." >&2; exit 1;
    }
    previous_gate=false
    [[ -f "$experiment/smoke-first" ]] && previous_gate="$(cat "$experiment/smoke-first")"
    [[ "$previous_gate" == "$smoke_first" ]] || {
        echo "Resume must preserve the original smoke-first setting." >&2; exit 2;
    }
    previous_staged=false
    [[ -f "$experiment/staged-development" ]] && previous_staged="$(cat "$experiment/staged-development")"
    [[ "$previous_staged" == "$staged" ]] || {
        echo "Resume must preserve the frozen staged-development policy." >&2; exit 2;
    }
    previous_guided=false
    [[ -f "$experiment/guided-followup" ]] && previous_guided="$(cat "$experiment/guided-followup")"
    [[ "$previous_guided" == "$guided" ]] || {
        echo "Resume must preserve the frozen guided-followup policy." >&2; exit 2;
    }
else
    if [[ -e "$experiment" ]]; then
        echo "Refusing an existing run directory. Use --resume for this exact frozen run." >&2
        exit 1
    fi
    mkdir -p "$experiment"/{source,launch,work,results,reports,controller_home}
    cp "$source_dir/"*.py "$source_dir/"*.nf "$source_dir/nextflow.config" "$experiment/source/"
    cp "$config" "$experiment/experiment.json"
    printf '%s\n' "$smoke_first" > "$experiment/smoke-first"
    printf '%s\n' "$staged" > "$experiment/staged-development"
    printf '%s\n' "$guided" > "$experiment/guided-followup"
    (cd "$experiment/source"; sha256sum *.py *.nf nextflow.config) > "$experiment/source.sha256"
fi
[[ "$guided" == false || -z "${SNPSTER_REUSE_MANIFEST:-}" ]] || {
    echo "Guided follow-up generates its own phase-only reuse policy; do not supply a native model reuse manifest." >&2; exit 2;
}
if [[ "$staged" == true ]]; then
    [[ -z "${SNPSTER_REUSE_MANIFEST:-}" ]] || {
        echo "Staged development generates separate frozen subset-reuse manifests per stage." >&2; exit 2;
    }
    source_run="${SNPSTER_STAGE_SOURCE_RUN:-}"
    [[ -n "$source_run" ]] || { echo "Set SNPSTER_STAGE_SOURCE_RUN to the stopped superseded run." >&2; exit 2; }
    source_run="$(realpath "$source_run")"
    mkdir -p "$experiment/runtime"
    if [[ "$resume" == true ]]; then
        [[ "$(cat "$experiment/runtime/staged-source-run")" == "$source_run" ]] || {
            echo "A resume cannot change its frozen staged source run." >&2; exit 2;
        }
    else
        printf '%s\n' "$source_run" > "$experiment/runtime/staged-source-run"
    fi
fi
reuse_manifest="${SNPSTER_REUSE_MANIFEST:-}"
if [[ -n "$reuse_manifest" ]]; then
    reuse_manifest="$(realpath "$reuse_manifest")"
    if [[ "$resume" == true ]]; then
        cmp --silent "$reuse_manifest" "$experiment/runtime/reuse-manifest.json" || {
            echo "A resume cannot add or change its frozen reuse manifest." >&2; exit 1;
        }
    else
        mkdir -p "$experiment/runtime"
        cp "$reuse_manifest" "$experiment/runtime/reuse-manifest.json"
        sha256sum "$experiment/runtime/reuse-manifest.json" > "$experiment/runtime/reuse-manifest.sha256"
    fi
fi
if [[ -f "$experiment/runtime/reuse-manifest.json" ]]; then
    sha256sum --check "$experiment/runtime/reuse-manifest.sha256"
    python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); assert c["java_profile"]==sys.argv[2], "Reuse requires the original JVM profile"' \
        "$experiment/runtime/reuse-manifest.json" "$java_profile"
fi
(cd "$experiment/source"; sha256sum --check "$experiment/source.sha256")
if [[ "$java_profile" == c1 ]]; then
    mkdir -p "$experiment/runtime"
    if [[ -f "$experiment/runtime/java-c1-pinned.config" ]]; then
        cmp --silent "$source_dir/java_c1.config" "$experiment/runtime/java-c1-pinned.config" || {
            echo "The recorded JVM recovery profile changed; refusing an ambiguous resume." >&2; exit 1;
        }
    else
        cp "$source_dir/java_c1.config" "$experiment/runtime/java-c1-pinned.config"
    fi
fi
printf '%s\n' "$java_profile" > "$experiment/java-profile"
suffix="$(date +%Y%m%d_%H%M%S_%N)"
name="snpster-ancestry-standard-$suffix"
printf '%s\n' "$name" > "$experiment/reports/$suffix.controller.name"
cp "$source_dir/launch.sh" "$experiment/reports/$suffix.launcher.sh"
sha256sum "$source_dir/launch.sh" > "$experiment/reports/$suffix.launcher.sha256"
if [[ "$java_profile" == c1 ]]; then
    sha256sum "$experiment/runtime/java-c1-pinned.config" > "$experiment/reports/$suffix.java-profile.sha256"
fi
options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v /srv/dependencies:/srv/dependencies:ro"
# Add read-only absolute input roots outside this repository using a colon-free, newline-delimited file.
mounts=()
if [[ -n "${SNPSTER_INPUT_ROOTS_FILE:-}" ]]; then
    while IFS= read -r root; do
        root="$(realpath "$root")"
        [[ -d "$root" && "$root" != / && "$root" =~ ^/[A-Za-z0-9_./-]+$ ]] || {
            echo "Invalid read-only input directory: $root" >&2; exit 2;
        }
        [[ "$root" == "$repository" || "$root" == /srv/dependencies || "$root" == "$experiment/source" ]] && continue
        mounts+=(-v "$root:$root:ro")
        options="$options -v $root:$root:ro"
    done < "$SNPSTER_INPUT_ROOTS_FILE"
fi
image='sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Resource budget: $cpu task CPUs / $budget GiB, plus a 2-GiB controller"
echo "Hard task caps: one CPU / $task_memory GiB; PCA bases $pca_memory GiB."
echo "JVM runtime profile: $java_profile"
if [[ "$java_profile" == c1 ]]; then
    echo "C2 disabled in a pinned Java task image; PHASE/native-INFER caches change, non-Java fits remain reusable."
fi
echo "Frozen run: $experiment"
[[ "$smoke_first" == true ]] && echo "Stages: unchanged-model chr22 smoke -> independent audit -> full autosomes -> independent audit."
[[ "$staged" == true ]] && echo "Stages: 96-person development plan -> chr22 screen/audit -> five validation-selected FLARE finalists + inexpensive baselines genome-wide/audit. No test stage."
[[ "$guided" == true ]] && echo "Stages: founder-grounded eight-source plan -> bounded chr20/22 guided screen -> independent audit -> EAS/SEA and GME report. No test or genome expansion."
set +e
docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e SUFFIX="$suffix" \
    -e BUDGET_CPUS="$cpu" -e BUDGET_MEMORY="$budget" -e TASK_MEMORY="$task_memory" -e PCA_MEMORY="$pca_memory" \
    -e CONTAINER_OPTIONS="$options" -e RESUME="$resume" -e SMOKE_FIRST="$smoke_first" -e JAVA_PROFILE="$java_profile" \
    -e STAGED_DEVELOPMENT="$staged" \
    -e GUIDED_FOLLOWUP="$guided" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    -v /srv/dependencies:/srv/dependencies:ro "${mounts[@]}" \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment/launch" \
    --entrypoint /bin/bash "$image" -euo pipefail -c '
        stty cols 180 rows 55
        docker run --rm --network none --cpus=1 --memory=128m \
            --entrypoint /bin/true sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
        echo "Task-container engine probe passed"
        extra=()
        [[ "$RESUME" == true ]] && extra=(-resume)
        runtime_config="$EXPERIMENT/source/nextflow.config"
        [[ "$JAVA_PROFILE" == c1 ]] && runtime_config="$runtime_config,$EXPERIMENT/runtime/java-c1-pinned.config"
        pipeline() {
            local stage=$1 config=$2 root="$EXPERIMENT" prefix="$SUFFIX"
            local reuse_args=()
            local workflow="$EXPERIMENT/source/main.nf"
            if [[ "$GUIDED_FOLLOWUP" == true ]]; then
                workflow="$EXPERIMENT/source/guided_followup.nf"
                reuse_args=(--guided_stage run --phase_reuse "$EXPERIMENT/planning/design/phase-reuse.json")
            fi
            if [[ "$stage" != full && -f "$EXPERIMENT/runtime/reuse-manifest.json" ]]; then
                reuse_args=(--reuse_manifest "$EXPERIMENT/runtime/reuse-manifest.json")
            fi
            if [[ "$STAGED_DEVELOPMENT" == true ]]; then
                case "$stage" in
                    screen) reuse_args=(--reuse_manifest "$EXPERIMENT/planning/plan/screen.reuse.json") ;;
                    genome) reuse_args=(--reuse_manifest "$EXPERIMENT/advancement/selection/genome.reuse.json") ;;
                    *) echo "Unexpected staged comparison phase: $stage" >&2; exit 2 ;;
                esac
            fi
            if [[ "$stage" != single ]]; then
                root="$EXPERIMENT/$stage"
                prefix="$SUFFIX.$stage"
            fi
            mkdir -p "$root"/{launch,work,results}
            cd "$root/launch"
            echo "=== $stage: Nextflow model comparison ==="
            nextflow -log "$EXPERIMENT/reports/$prefix.nextflow.log" \
            -C "$runtime_config" run "$workflow" \
            "${extra[@]}" -pool-size "$BUDGET_CPUS" -ansi-log true -work-dir "$root/work" \
            -with-trace "$EXPERIMENT/reports/$prefix.trace.tsv" \
            -with-report "$EXPERIMENT/reports/$prefix.report.html" \
            -with-timeline "$EXPERIMENT/reports/$prefix.timeline.html" \
            -with-dag "$EXPERIMENT/reports/$prefix.pipeline.html" \
            --experiment "$config" --output "$root/results" \
            --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" \
            --task_memory "$TASK_MEMORY GB" --pca_memory "$PCA_MEMORY GB" --container_options "$CONTAINER_OPTIONS" \
            "${reuse_args[@]}"
        }
        audit() {
            local stage=$1
            echo "=== Independent durable-output audit: $stage ==="
            docker run --rm --network none --cpus=1 --memory="${TASK_MEMORY}g" \
                $CONTAINER_OPTIONS -v "$EXPERIMENT:$EXPERIMENT:rw" \
                -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -w "$EXPERIMENT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                audit.py "$EXPERIMENT/$stage/results" 2>&1 | tee "$EXPERIMENT/reports/$SUFFIX.$stage.independent-audit.log"
        }
        metadata() {
            local stage=$1 selection=$2
            local root="$EXPERIMENT/$stage" prefix="$SUFFIX.$stage"
            mkdir -p "$root"/{launch,work}
            cd "$root/launch"
            nextflow -log "$EXPERIMENT/reports/$prefix.nextflow.log" \
                -C "$runtime_config" run "$EXPERIMENT/source/staged.nf" --metadata_stage "$selection" \
                "${extra[@]}" -pool-size "$BUDGET_CPUS" -ansi-log true -work-dir "$root/work" \
                -with-trace "$EXPERIMENT/reports/$prefix.trace.tsv" \
                -with-report "$EXPERIMENT/reports/$prefix.report.html" \
                -with-timeline "$EXPERIMENT/reports/$prefix.timeline.html" \
                -with-dag "$EXPERIMENT/reports/$prefix.pipeline.html" \
                --experiment "$EXPERIMENT/experiment.json" --output "$root" \
                --source_run "$(cat "$EXPERIMENT/runtime/staged-source-run")" \
                --plan_dir "$EXPERIMENT/planning/plan" --screen_results "$EXPERIMENT/screen/results" \
                --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" \
                --task_memory "$TASK_MEMORY GB" --container_options "$CONTAINER_OPTIONS"
        }
        guided_metadata() {
            local stage=$1 selection=$2
            local root="$EXPERIMENT/$stage" prefix="$SUFFIX.$stage"
            mkdir -p "$root"/{launch,work}
            cd "$root/launch"
            nextflow -log "$EXPERIMENT/reports/$prefix.nextflow.log" \
                -C "$runtime_config" run "$EXPERIMENT/source/guided_followup.nf" --guided_stage "$selection" \
                "${extra[@]}" -pool-size "$BUDGET_CPUS" -ansi-log true -work-dir "$root/work" \
                -with-trace "$EXPERIMENT/reports/$prefix.trace.tsv" \
                -with-report "$EXPERIMENT/reports/$prefix.report.html" \
                -with-timeline "$EXPERIMENT/reports/$prefix.timeline.html" \
                -with-dag "$EXPERIMENT/reports/$prefix.pipeline.html" \
                --experiment "$EXPERIMENT/experiment.json" --output "$root" \
                --design "$EXPERIMENT/planning/design" \
                --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" \
                --task_memory "$TASK_MEMORY GB" --container_options "$CONTAINER_OPTIONS"
        }
        if [[ "$GUIDED_FOLLOWUP" == true ]]; then
            guided_metadata planning plan
            cat "$EXPERIMENT/planning/design/provenance.json"
            pipeline screen "$EXPERIMENT/planning/design/config.json"
            audit screen
            guided_metadata reporting report
            cat "$EXPERIMENT/reporting/discrimination/source-discrimination.json"
            docker run --rm --network none --cpus=1 --memory="${TASK_MEMORY}g" \
                $CONTAINER_OPTIONS -v "$EXPERIMENT:$EXPERIMENT:rw" \
                -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -w "$EXPERIMENT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                guided_followup.py cleanup --results "$EXPERIMENT" \
                --output "$EXPERIMENT/reports/scoped-cache-cleanup.json"
            echo "GUIDED SCREEN, EAS/SEA REPORT AND INDEPENDENT AUDIT PASSED. CACHES CLEANED. TEST REMAINS WITHHELD."
        elif [[ "$STAGED_DEVELOPMENT" == true ]]; then
            metadata planning plan
            pipeline screen "$EXPERIMENT/planning/plan/screen.json"
            audit screen
            metadata advancement advance
            echo "SCREEN AUDITED; ADVANCING EXACTLY FIVE VALIDATION-SELECTED FLARE CONFIGURATIONS."
            cat "$EXPERIMENT/advancement/selection/finalists.json"
            pipeline genome "$EXPERIMENT/advancement/selection/genome.json"
            audit genome
            echo "GENOME-WIDE DEVELOPMENT COMPARISON AND INDEPENDENT AUDIT PASSED. TEST REMAINS WITHHELD."
            docker run --rm --network none --cpus=1 --memory="${TASK_MEMORY}g" \
                $CONTAINER_OPTIONS -v "$EXPERIMENT:$EXPERIMENT:rw" \
                -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -w "$EXPERIMENT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                staged_comparison.py cleanup --source-run "$EXPERIMENT" \
                --output "$EXPERIMENT/reports/scoped-cache-cleanup.json"
            echo "EXACT SUCCESSFUL STAGE CACHES REMOVED; DURABLE SCIENTIFIC FILES REVERIFIED."
        elif [[ "$SMOKE_FIRST" == true ]]; then
            docker run --rm --network none --cpus=1 --memory=2g \
                $CONTAINER_OPTIONS -v "$EXPERIMENT:$EXPERIMENT:rw" \
                -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                -w "$EXPERIMENT/source" --entrypoint python3 \
                sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                prepare_comparison.py smoke --config "$EXPERIMENT/experiment.json" --output "$EXPERIMENT/smoke.json"
            pipeline smoke "$EXPERIMENT/smoke.json"
            audit smoke
            echo "REAL CHR22 MODEL SMOKE AND INDEPENDENT AUDIT PASSED; STARTING FULL AUTOSOMES."
            pipeline full "$EXPERIMENT/experiment.json"
            audit full
            echo "FULL WORKFLOW AND INDEPENDENT AUDIT PASSED."
        else
            pipeline single "$EXPERIMENT/experiment.json"
        fi
    '
status=$?
set -e
docker logs "$name" > "$experiment/reports/$suffix.controller.log" 2>&1
docker inspect "$name" > "$experiment/reports/$suffix.controller.json"
printf '%s\n' "$status" > "$experiment/reports/$suffix.controller.exit-status"
if [[ $status -eq 0 && "$smoke_first" == false && "$staged" == false && "$guided" == false ]]; then
    echo "Independently verifying durable truth, posterior integration, aggregation and metrics..."
    if docker run --rm --network none --cpus=1 --memory="${task_memory}g" \
        --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
        -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
        -v /srv/dependencies:/srv/dependencies:ro "${mounts[@]}" \
        -w "$experiment/source" --entrypoint python3 \
        sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
        audit.py "$experiment/results" 2>&1 | tee "$experiment/reports/$suffix.independent-audit.log"; then
        echo "Independent durable-output audit passed."
    else
        echo "Independent durable-output audit failed; retaining this run's work/history." >&2
        status=1
    fi
fi
printf '%s\n' "$status" > "$experiment/reports/$suffix.exit-status"
if [[ $status -eq 0 ]]; then
    echo "Workflow and durable-output audit passed. Retain history for resume, or clean only this exact successful run."
else
    echo "Workflow failed ($status); work/history retained for exact-run resume." >&2
fi
exit "$status"
