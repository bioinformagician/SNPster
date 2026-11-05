#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
source="$project/results_v3_observed_dense_01"
recovery="${OBSERVED_EXPANSION_RECOVERY:-}"
test -s "$source/completion.json"
test -s "$source/recovery_state.json"
if [[ -n "$recovery" ]]; then
    recovery="$(realpath "$recovery")"
    [[ "$recovery" == "$project"/results_v3_observed_dense_model_expansion_*/smoke/work/*/result ]]
    test -d "$recovery"
fi
suffix="$(date +%Y%m%d_%H%M%S)"
experiment="$project/results_v3_observed_dense_model_expansion_$suffix"
name="snpster-observed-model96-$suffix"
mkdir -p "$experiment/smoke/launch" "$experiment/full/launch" "$experiment/controller_home"
sha256sum \
    "$project/observed_dense_model_expansion.nf" \
    "$project/observed_dense_model_expansion.config" \
    "$project/run_observed_dense_model_expansion.sh" \
    "$project/config/dense_model_expansion.json" \
    "$project/config/observed_dense_protocol.json" \
    "$project/pedigree_benchmark/scripts/observed_dense_model_expansion.py" \
    "$project/pedigree_benchmark/scripts/benchmark_dense_observed.py" \
    "$source/completion.json" "$source/recovery_state.json" \
    "$HOME/.nextflow/framework/26.04.4/nextflow-26.04.4-one.jar" > "$experiment/source_sha256.txt"
printf '%s\n' "$name" > "$experiment/controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=22 --memory=180g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PROJECT="$project" -e REPOSITORY="$repository" \
    -e EXPERIMENT="$experiment" -e SOURCE="$source" -e RECOVERY="$recovery" \
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
            if [[ "$mode" == smoke && -n "$RECOVERY" ]]; then
                recovery_args=(--recovery "$RECOVERY")
            fi
            cd "$scope/launch"
            nextflow -log "$scope/nextflow.log" -C "$PROJECT/observed_dense_model_expansion.config" \
                run "$PROJECT/observed_dense_model_expansion.nf" -resume -ansi-log true \
                -work-dir "$scope/work" -with-trace "$scope/trace.tsv" \
                -with-report "$scope/report.html" -with-timeline "$scope/timeline.html" \
                -with-dag "$scope/pipeline.html" --project "$PROJECT" --run_dir "$scope" \
                --source "$SOURCE" --plan "$PROJECT/config/dense_model_expansion.json" \
                --smoke "$smoke" "${recovery_args[@]}" --container_options "$options"
        done
    ' > "$experiment/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $experiment"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
