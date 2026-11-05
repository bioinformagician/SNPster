#!/usr/bin/env bash
set -euo pipefail

main() {
    [[ $# -eq 2 ]] || { echo "Usage: bash run_master_resources.sh ABSOLUTE_SPECIFICATION ABSOLUTE_NEW_RUN" >&2; return 2; }
    local source_dir repository specification design_name design prepared experiment cpu budget phase name status options
    source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repository="$(cd "$source_dir/../.." && pwd)"
    specification="$(realpath "$1")"
    design_name="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["founder_design"])' "$specification")"
    [[ "$design_name" =~ ^[A-Za-z0-9_.-]+\.json$ ]] || {
        echo "founder_design must name a JSON file beside the specification, without directory components." >&2; return 2;
    }
    design="$(realpath "$(dirname "$specification")/$design_name")"
    experiment="$(realpath -m "$2")"
    prepared="${SNPSTER_PREPARED_RESOURCES:-}"
    if [[ -n "$prepared" ]]; then
        prepared="$(realpath "$prepared")"
        [[ "$prepared" =~ ^/[A-Za-z0-9_./-]+$ && -f "$prepared/kinship/screen.json" ]] || {
            echo "Require an existing absolute all-resource preparation with genome-wide evidence." >&2; return 2;
        }
    fi
    cpu="${SNPSTER_CPUS:-12}"
    budget="${SNPSTER_MEMORY_GB:-96}"
    [[ "$cpu" =~ ^[1-9][0-9]*$ && "$cpu" -ge 4 && "$cpu" -le 22 &&
       "$budget" =~ ^[1-9][0-9]*$ && "$budget" -ge 32 && "$budget" -le 178 ]] || {
        echo "Require 4-22 task CPUs and 32-178 GiB memory; KING uses four CPUs/32 GiB." >&2; return 2;
    }
    for path in "$specification" "$design" "$experiment" "$repository"; do
        [[ "$path" =~ ^/[A-Za-z0-9_./-]+$ && "$path" != / && "$path" != "$HOME" ]] || {
            echo "Invalid absolute path: $path" >&2; return 2;
        }
    done
    [[ ! -e "$experiment" && "$experiment" != "$repository" ]] || { echo "Require a new dedicated run directory." >&2; return 2; }
    mkdir -p "$experiment"/{source,inputs,smoke,full,controller_home}
    cp "$source_dir/"*.py "$source_dir/master_resources.nf" "$source_dir/master_panels.nf" \
        "$source_dir/master_resources.config" "$experiment/source/"
    cp "${BASH_SOURCE[0]}" "$experiment/source/run_master_resources.sh"
    cp "$repository/data_analysis/FLARE/SNPster_v3_100model_pedigree/pedigree_benchmark/scripts/"{convert_plink_map,realized_ancestry,verify_founder_assignments}.py "$experiment/source/"
    cp "$specification" "$experiment/inputs/specification.json"
    cp "$design" "$experiment/inputs/$design_name"
    (cd "$experiment/source"; sha256sum *.py *.nf *.sh *.config) > "$experiment/source.sha256"
    sha256sum "$experiment/inputs/"*.json > "$experiment/inputs.sha256"
    options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v /srv/dependencies:/srv/dependencies:ro"
    for phase in smoke full; do
        mkdir -p "$experiment/$phase"/{launch,work,release,reports}
        name="snpster-resources-$phase-$(basename "$experiment")"
        echo "$name" > "$experiment/$phase/reports/controller.name"
        echo "Monitor: docker logs -f $name"
        echo "$phase: $cpu task CPUs/$budget GiB; one-core preparation, four-core KING; no ancestry-model fitting."
        set +e
        docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
            --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
            -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
            -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e PHASE="$phase" \
            -e PREPARED_RESOURCES="$prepared" \
            -e BUDGET_CPUS="$cpu" -e BUDGET_MEMORY="$budget" -e CONTAINER_OPTIONS="$options" \
            -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" -v /srv/dependencies:/srv/dependencies:ro \
            -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
            -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment/$phase/launch" \
            --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
            -euo pipefail -c '
                stty cols 180 rows 55
                (cd "$EXPERIMENT/source"; sha256sum --check "$EXPERIMENT/source.sha256")
                sha256sum --check "$EXPERIMENT/inputs.sha256"
                docker run --rm --network none --cpus=1 --memory=128m --entrypoint /bin/true \
                    sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
                smoke=false
                [[ "$PHASE" == smoke ]] && smoke=true
                reuse=()
                if [[ -n "$PREPARED_RESOURCES" ]]; then
                    echo "Reusing audited genotype preparation only; source hashes and all-pair coverage will be checked. Roles and descendants are new."
                    reuse=(--prepared_resources "$PREPARED_RESOURCES")
                fi
                nextflow -log "$EXPERIMENT/$PHASE/reports/nextflow.log" \
                    -C "$EXPERIMENT/source/master_resources.config" run "$EXPERIMENT/source/master_resources.nf" \
                    -ansi-log true -work-dir "$EXPERIMENT/$PHASE/work" \
                    -with-trace "$EXPERIMENT/$PHASE/reports/trace.tsv" \
                    -with-report "$EXPERIMENT/$PHASE/reports/report.html" \
                    -with-timeline "$EXPERIMENT/$PHASE/reports/timeline.html" \
                    -with-dag "$EXPERIMENT/$PHASE/reports/pipeline.html" \
                    --specification "$EXPERIMENT/inputs/specification.json" \
                    --output "$EXPERIMENT/$PHASE/release" --smoke_only "$smoke" \
                    --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" \
                    --container_options "$CONTAINER_OPTIONS" "${reuse[@]}"
                if [[ "$PHASE" == full ]]; then
                    echo "Independently checking durable source hashes, all people, biological splits and actual model-input interfaces..."
                    docker run --rm --network none --cpus=1 --memory=8g \
                        --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                        -v "'"$repository:$repository"':ro" -v "$EXPERIMENT:$EXPERIMENT:rw" \
                        -v /srv/dependencies:/srv/dependencies:ro -w "$EXPERIMENT/source" --entrypoint python3 \
                        sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                        audit_resources.py verify --dataset "$EXPERIMENT/full/release" \
                        --output "$EXPERIMENT/full/release/verification/independent.json"
                    echo "VERIFIED: full-resource dataset, persistent pedigree VCFs and platform-specific model configurations."
                else
                    if [[ -n "$PREPARED_RESOURCES" ]]; then
                        echo "PINNED-TASK IMPORT SMOKE PASSED: original 22-autosome GT audits and all-pair evidence are unchanged."
                    else
                        echo "REAL CHR22/CHR4 ALL-RESOURCE MEMORY-SCALE SMOKE AND ORIGINAL-GENOTYPE AUDITS PASSED."
                    fi
                fi
            '
        status=$?
        set -e
        docker logs "$name" > "$experiment/$phase/reports/controller.log" 2>&1
        docker inspect "$name" > "$experiment/$phase/reports/controller.json"
        echo "$status" > "$experiment/$phase/reports/exit-status"
        [[ "$status" -eq 0 ]] || { echo "$phase FAILED; frozen source and resumable work preserved." >&2; return "$status"; }
    done
    echo "COMPLETE AND VERIFIED: $experiment/full/release"
    echo "Only this verified run is eligible for scoped disposable-work cleanup."
}

main "$@"
