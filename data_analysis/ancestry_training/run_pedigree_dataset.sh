#!/usr/bin/env bash
set -euo pipefail

main() {
    if [[ $# -ne 3 ]]; then
        echo "Usage: bash run_pedigree_dataset.sh ABSOLUTE_VERIFIED_RELEASE ABSOLUTE_DESIGN_JSON ABSOLUTE_NEW_RUN" >&2
        return 2
    fi
    local source_dir repository parent design experiment cpu budget vendor path name phase status options
    source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repository="$(cd "$source_dir/../.." && pwd)"
    parent="$(realpath "$1")"
    design="$(realpath "$2")"
    experiment="$(realpath -m "$3")"
    cpu="${SNPSTER_CPUS:-12}"
    budget="${SNPSTER_MEMORY_GB:-96}"
    for path in "$parent" "$design" "$experiment" "$repository"; do
        [[ "$path" =~ ^/[A-Za-z0-9_./-]+$ && "$path" != / && "$path" != "$HOME" ]] || {
            echo "Invalid absolute path: $path" >&2; return 2;
        }
    done
    [[ "$cpu" =~ ^[1-9][0-9]*$ && "$cpu" -le 22 && "$budget" =~ ^[1-9][0-9]*$ && "$budget" -ge 8 && "$budget" -le 178 ]] || {
        echo "Require 1-22 task cores and 8-178 GiB task memory." >&2; return 2;
    }
    [[ -f "$parent/verification/independent.json" && ! -e "$experiment" && "$experiment" != "$repository" ]] || {
        echo "Require a verified source release and a new dedicated output directory." >&2; return 2;
    }
    mkdir -p "$experiment"/{source,inputs,smoke,full,controller_home}
    cp "$source_dir/"*.py "$source_dir/prepare_pedigree.nf" "$source_dir/nextflow.config" "$experiment/source/"
    cp "${BASH_SOURCE[0]}" "$experiment/source/run_pedigree_dataset.sh"
    vendor="$repository/data_analysis/FLARE/SNPster_v3_100model_pedigree/pedigree_benchmark/scripts"
    cp "$vendor/convert_plink_map.py" "$vendor/realized_ancestry.py" "$vendor/verify_founder_assignments.py" "$experiment/source/"
    cp "$design" "$experiment/inputs/design.json"
    (cd "$experiment/source"; sha256sum *.py *.nf *.sh nextflow.config) > "$experiment/source.sha256"
    sha256sum "$experiment/inputs/design.json" > "$experiment/inputs.sha256"
    options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v /srv/dependencies:/srv/dependencies:ro"
    for phase in smoke full; do
        mkdir -p "$experiment/$phase"/{launch,work,release,reports}
        name="snpster-pedigree-$phase-$(basename "$experiment")"
        printf '%s\n' "$name" > "$experiment/$phase/reports/controller.name"
        echo "Controller: $name"
        echo "Monitor: docker logs -f $name"
        echo "Stage: $phase; $cpu one-core slots, $budget GiB tasks plus 2-GiB controller; no model fitting."
        set +e
        docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
            --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
            -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
            -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e PARENT="$parent" -e PHASE="$phase" \
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
                nextflow -log "$EXPERIMENT/$PHASE/reports/nextflow.log" \
                    -C "$EXPERIMENT/source/nextflow.config" run "$EXPERIMENT/source/prepare_pedigree.nf" \
                    -ansi-log true -work-dir "$EXPERIMENT/$PHASE/work" \
                    -with-trace "$EXPERIMENT/$PHASE/reports/trace.tsv" \
                    -with-report "$EXPERIMENT/$PHASE/reports/report.html" \
                    -with-timeline "$EXPERIMENT/$PHASE/reports/timeline.html" \
                    -with-dag "$EXPERIMENT/$PHASE/reports/pipeline.html" \
                    --parent "$PARENT" --design "$EXPERIMENT/inputs/design.json" \
                    --output "$EXPERIMENT/$PHASE/release" --smoke_only "$smoke" \
                    --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" --task_memory "8 GB" \
                    --container_options "$CONTAINER_OPTIONS"
                if [[ "$PHASE" == full ]]; then
                    echo "Verifying durable hashes, biological families, three-generation ancestry and actual pipeline interfaces..."
                    docker run --rm --network none --cpus=1 --memory=8g \
                        --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
                        -v "'"$repository:$repository"':ro" -v "$EXPERIMENT:$EXPERIMENT:rw" \
                        -v /srv/dependencies:/srv/dependencies:ro -w "$EXPERIMENT/source" --entrypoint python3 \
                        sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
                        audit_pedigree.py verify --dataset "$EXPERIMENT/full/release" \
                        --output "$EXPERIMENT/full/release/verification/independent.json"
                    echo "VERIFIED: persistent phased founder/descendant VCFs, exact inheritance and model-ready inputs. No models fitted."
                else
                    echo "REAL CHR22 SMOKE AND INDEPENDENT EXACT-INHERITANCE AUDIT PASSED."
                fi
            '
        status=$?
        set -e
        docker logs "$name" > "$experiment/$phase/reports/controller.log" 2>&1
        docker inspect "$name" > "$experiment/$phase/reports/controller.json"
        printf '%s\n' "$status" > "$experiment/$phase/reports/exit-status"
        [[ "$status" -eq 0 ]] || { echo "$phase FAILED; frozen source, logs and work retained." >&2; return "$status"; }
    done
    echo "COMPLETED AND VERIFIED: $experiment/full/release"
    echo "Successful work/cache cleanup is permitted only after independent published-output verification."
}

main "$@"
