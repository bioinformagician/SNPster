#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo "Usage: bash run_flare_diagnostic.sh ABSOLUTE_PARENT_RUN ABSOLUTE_NEW_RUN" >&2
    exit 2
fi
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
self="$(realpath "${BASH_SOURCE[0]}")"
if [[ "${SNPSTER_DIAGNOSTIC_FROZEN_WRAPPER:-}" != "$self" || "$self" != "$source_dir"/.flare-launch.*.sh ]]; then
    snapshot="$(mktemp "$source_dir/.flare-launch.XXXXXXXX.sh")"
    cp "$self" "$snapshot"
    export SNPSTER_DIAGNOSTIC_FROZEN_WRAPPER="$snapshot"
    exec bash "$snapshot" "$@"
fi
readonly snapshot_file="$self"
trap 'rm -- "$snapshot_file"' EXIT
repository="$(cd "$source_dir/../.." && pwd)"
parent="$(realpath "$1")"
experiment="$(realpath -m "$2")"
for path in "$parent" "$experiment" "$repository"; do
    [[ "$path" =~ ^/[A-Za-z0-9_./-]+$ && "$path" != / && "$path" != "$HOME" ]] || {
        echo "Invalid mount path: $path" >&2; exit 2;
    }
done
[[ -f "$parent/full/results/summary/summary/summary.json" && ! -e "$experiment" ]] || {
    echo "Require a completed parent workflow and a new diagnostic directory." >&2; exit 2;
}
mkdir -p "$experiment"/{source,launch,work,results,reports,controller_home}
cp "$source_dir/"*.py "$source_dir/diagnose_flare.nf" "$source_dir/nextflow.config" "$experiment/source/"
cp "${SNPSTER_DIAGNOSTIC_PROTOCOL:-$source_dir/flare_diagnostic_v1.json}" "$experiment/protocol.json"
(cd "$experiment/source"; sha256sum *.py *.nf nextflow.config) > "$experiment/source.sha256"
suffix="$(date +%Y%m%d_%H%M%S)"
name="${SNPSTER_DIAGNOSTIC_NAME:-snpster-flare-diagnosis-$suffix}"
cpu="${SNPSTER_DIAGNOSTIC_CPUS:-14}"
[[ "$cpu" =~ ^[1-9][0-9]*$ && "$cpu" -le 22 ]] || { echo "Require 1-22 task cores." >&2; exit 2; }
memory=$((cpu * 4))
printf '%s\n' "$name" > "$experiment/reports/controller.name"
cp "$self" "$experiment/reports/launcher.sh"
(cd "$experiment/reports"; sha256sum launcher.sh) > "$experiment/reports/launcher.sha256"
options="--user $(id -u):$(id -g) -v $repository:$repository:ro -v $experiment/source:$experiment/source:ro -v /srv/dependencies:/srv/dependencies:ro"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Diagnostic: $experiment"
echo "Budget: $cpu one-core/4-GiB tasks, $memory-GiB task budget plus 2-GiB controller; parent audit left unchanged."
set +e
docker run -t --name "$name" --init --network none --cpus="$cpu" --memory=2g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e EXPERIMENT="$experiment" -e PARENT="$parent" \
    -e CONTAINER_OPTIONS="$options" \
    -e BUDGET_CPUS="$cpu" -e BUDGET_MEMORY="$memory" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment/launch" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 55
        (cd "$EXPERIMENT/source"; sha256sum --check "$EXPERIMENT/source.sha256")
        docker run --rm --network none --cpus=1 --memory=128m --entrypoint /bin/true \
            sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d
        nextflow -log "$EXPERIMENT/reports/nextflow.log" \
            -C "$EXPERIMENT/source/nextflow.config" run "$EXPERIMENT/source/diagnose_flare.nf" \
            -ansi-log true -work-dir "$EXPERIMENT/work" \
            -with-trace "$EXPERIMENT/reports/trace.tsv" \
            -with-report "$EXPERIMENT/reports/report.html" \
            -with-timeline "$EXPERIMENT/reports/timeline.html" \
            -with-dag "$EXPERIMENT/reports/pipeline.html" \
            --experiment "$EXPERIMENT/protocol.json" --parent "$PARENT" --output "$EXPERIMENT/results" \
            --budget_cpus "$BUDGET_CPUS" --budget_memory "$BUDGET_MEMORY GB" --task_memory "4 GB" \
            --container_options "$CONTAINER_OPTIONS"
        echo "Diagnostic workflow and exact native-output verification passed."
    '
status=$?
set -e
docker logs "$name" > "$experiment/reports/controller.log" 2>&1
docker inspect "$name" > "$experiment/reports/controller.json"
printf '%s\n' "$status" > "$experiment/reports/exit-status"
echo "Diagnostic controller exit: $status; published evidence and exact failed work retained."
exit "$status"
