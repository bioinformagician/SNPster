#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 || ( $# -eq 3 && "$3" != --resume ) ]]; then
    echo "Usage: bash run_kinship_dataset.sh ABSOLUTE_COMPLETED_PARENT ABSOLUTE_NEW_DATASET_RUN [--resume]" >&2
    exit 2
fi
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
self="$(realpath "${BASH_SOURCE[0]}")"
if [[ "${SNPSTER_COHORT_FROZEN_WRAPPER:-}" != "$self" || "$self" != "$source_dir"/.cohort-launch.*.sh ]]; then
    snapshot="$(mktemp "$source_dir/.cohort-launch.XXXXXXXX.sh")"
    cp "$self" "$snapshot"
    export SNPSTER_COHORT_FROZEN_WRAPPER="$snapshot"
    exec bash "$snapshot" "$@"
fi
readonly snapshot_file="$self"
trap 'rm -- "$snapshot_file"' EXIT
repository="$(cd "$source_dir/../.." && pwd)"
parent="$(realpath "$1")"
experiment="$(realpath -m "$2")"
dense=/srv/dependencies/imputation_runner/imputer/beagle_references
cpu="${SNPSTER_CPUS:-20}"
budget="${SNPSTER_MEMORY_GB:-176}"
[[ "$cpu" =~ ^[1-9][0-9]*$ && "$cpu" -le 22 && "$budget" =~ ^[1-9][0-9]*$ && "$budget" -ge 8 && "$budget" -le 178 ]] || {
    echo "Require 1-22 task CPUs and 8-178 GiB task memory." >&2; exit 2;
}
for path in "$parent" "$experiment" "$repository"; do
    [[ "$path" =~ ^/[A-Za-z0-9_./-]+$ && "$path" != / && "$path" != "$HOME" ]] || {
        echo "Invalid absolute mount path: $path" >&2; exit 2;
    }
done
[[ -f "$parent/full/results/independent_verification.json" && "$experiment" != "$repository" && "$experiment" != "$source_dir" ]] || {
    echo "Require an independently verified completed parent and a dedicated output directory." >&2; exit 2;
}
resume=false
if [[ "${3:-}" == --resume ]]; then
    resume=true
    [[ -f "$experiment/source.sha256" && "$(cat "$experiment/parent.path")" == "$parent" ]] || {
        echo "Missing or mismatched exact-run history." >&2; exit 2;
    }
else
    [[ ! -e "$experiment" ]] || { echo "Output exists; choose a new run or exact --resume." >&2; exit 2; }
    mkdir -p "$experiment"/{source,inputs,launch,work,release,reports,controller_home}
    cp "$source_dir/"*.py "$source_dir/prepare_kinship.nf" "$source_dir/nextflow.config" "$experiment/source/"
    curl --fail --location --silent --show-error --max-time 120 \
        https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped \
        -o "$experiment/inputs/1000genomes.ped"
    printf '%s  %s\n' 38687335d053056d799891a37c42b1ac24a5a2a1eb5fccaf34db58188bccfff7 \
        "$experiment/inputs/1000genomes.ped" | sha256sum --check
    cp /srv/dependencies/ancestry/flare2_test/source_model/ancestry_role_audit.tsv "$experiment/inputs/ancestry_role_audit.tsv"
    printf '%s\n' "$parent" > "$experiment/parent.path"
    (cd "$experiment/source"; sha256sum *.py *.nf nextflow.config) > "$experiment/source.sha256"
    (cd "$experiment/inputs"; sha256sum *.ped *.tsv) > "$experiment/inputs.sha256"
fi
(cd "$experiment/source"; sha256sum --check "$experiment/source.sha256")
suffix="$(date +%Y%m%d_%H%M%S_%N)"
name="${SNPSTER_COHORT_NAME:-snpster-cohort-kinship-$suffix}"
printf '%s\n' "$name" > "$experiment/reports/$suffix.controller.name"
cp "$self" "$experiment/reports/$suffix.launcher.sh"
options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v /srv/dependencies:/srv/dependencies:ro"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Dataset: $experiment/release"
echo "Budget: $cpu one-core task slots, $budget GiB task budget plus 2-GiB controller."
echo "Stages: lineage catalog -> real chr22 smoke -> 22-autosome KING -> whole-family split -> reusable genotypes/metadata -> independent audit."
echo "No model training or sweep will be launched."
set +e
docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e PARENT="$parent" -e DENSE="$dense" \
    -e BUDGET_CPUS="$cpu" -e BUDGET_MEMORY="$budget" -e CONTAINER_OPTIONS="$options" -e RESUME="$resume" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment/launch" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 55
        (cd "$EXPERIMENT/inputs"; sha256sum --check "$EXPERIMENT/inputs.sha256")
        docker run --rm --network none --cpus=1 --memory=128m --entrypoint /bin/true \
            sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
        extra=()
        [[ "$RESUME" == true ]] && extra=(-resume)
        nextflow -log "$EXPERIMENT/reports/nextflow.log" \
            -C "$EXPERIMENT/source/nextflow.config" run "$EXPERIMENT/source/prepare_kinship.nf" \
            "${extra[@]}" -ansi-log true -work-dir "$EXPERIMENT/work" \
            -with-trace "$EXPERIMENT/reports/trace.tsv" -with-report "$EXPERIMENT/reports/report.html" \
            -with-timeline "$EXPERIMENT/reports/timeline.html" -with-dag "$EXPERIMENT/reports/pipeline.html" \
            --parent "$PARENT" --dense "$DENSE" --pedigree "$EXPERIMENT/inputs/1000genomes.ped" \
            --authority "$EXPERIMENT/inputs/ancestry_role_audit.tsv" --output "$EXPERIMENT/release" \
            --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" --task_memory "8 GB" \
            --container_options "$CONTAINER_OPTIONS"
        echo "Independently auditing published metadata, relatives, ancestry, hashes and both model-pipeline input routes..."
        docker run --rm --network none --cpus=1 --memory=8g \
            --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
            -v "'"$repository:$repository"':ro" -v "$EXPERIMENT:$EXPERIMENT:rw" \
            -v /srv/dependencies:/srv/dependencies:ro -w "$EXPERIMENT/source" --entrypoint python3 \
            sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d \
            audit_kinship.py verify --dataset "$EXPERIMENT/release" --interfaces \
            --output "$EXPERIMENT/release/verification/independent.json"
        echo "REUSABLE DATASET AND PIPELINE INPUT AUDIT PASSED. No model fitting was run."
    '
status=$?
set -e
docker logs "$name" > "$experiment/reports/$suffix.controller.log" 2>&1
docker inspect "$name" > "$experiment/reports/$suffix.controller.json"
printf '%s\n' "$status" > "$experiment/reports/$suffix.exit-status"
echo "Dataset controller exit: $status; work/history retained until independent verification and scoped cleanup."
exit "$status"
