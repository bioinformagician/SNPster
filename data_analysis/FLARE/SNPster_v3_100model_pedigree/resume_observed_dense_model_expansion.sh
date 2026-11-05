#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
if [[ $# -ne 2 ]]; then
    echo "Usage: bash resume_observed_dense_model_expansion.sh ABSOLUTE_EXPERIMENT ABSOLUTE_SMOKE_RECOVERY" >&2
    exit 1
fi
experiment="$(realpath "$1")"
recovery="$(realpath "$2")"
[[ "$experiment" == "$project"/results_v3_observed_dense_model_expansion_* ]]
[[ "$recovery" == "$project"/results_v3_observed_dense_model_expansion_*/smoke/work/*/result ]]
test -d "$experiment/smoke/launch"
test -d "$experiment/full/launch"
test -s "$recovery/flare_full_array_fixed_g007/HGDP01374_inherited_phase.model"
source="$project/results_v3_observed_dense_01"
suffix="$(date +%Y%m%d_%H%M%S)"
name="snpster-observed-model96-resume-$suffix"
printf '%s\n' "$name" > "$experiment/resume-$suffix.controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=22 --memory=180g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PROJECT="$project" -e REPOSITORY="$repository" \
    -e EXPERIMENT="$experiment" -e SOURCE="$source" -e RECOVERY="$recovery" -e SUFFIX="$suffix" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    -v "$source:$source:ro" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 60
        options="--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v $SOURCE:$SOURCE:ro -v /srv/dependencies:/srv/dependencies:ro -v /srv/dependencies/ancestry:/data:ro"
        for mode in smoke full; do
            scope="$EXPERIMENT/$mode"
            [[ "$mode" == smoke ]] && smoke=true || smoke=false
            recovery_args=()
            if [[ "$mode" == smoke ]]; then
                recovery_args=(--recovery "$RECOVERY")
            fi
            cd "$scope/launch"
            nextflow -log "$scope/nextflow-resume-$SUFFIX.log" \
                -C "$PROJECT/observed_dense_model_expansion.config" \
                run "$PROJECT/observed_dense_model_expansion.nf" -resume -ansi-log true \
                -work-dir "$scope/work" -with-trace "$scope/trace-resume-$SUFFIX.tsv" \
                -with-report "$scope/report-resume-$SUFFIX.html" \
                -with-timeline "$scope/timeline-resume-$SUFFIX.html" \
                -with-dag "$scope/pipeline-resume-$SUFFIX.html" \
                --project "$PROJECT" --run_dir "$scope" --source "$SOURCE" \
                --plan "$PROJECT/config/dense_model_expansion.json" --smoke "$smoke" \
                "${recovery_args[@]}" --container_options "$options"
        done
    ' > "$experiment/resume-$suffix.controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $experiment"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
