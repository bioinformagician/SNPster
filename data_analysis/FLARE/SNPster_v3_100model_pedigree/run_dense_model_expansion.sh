#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
if [[ $# -ne 1 ]]; then
    echo "Usage: bash run_dense_model_expansion.sh ABSOLUTE_SOURCE_EXPERIMENT" >&2; exit 1
fi
source_experiment="$(realpath "$1")"
[[ "$source_experiment" == "$project"/results_v3_dense_split_* ]]
test -s "$source_experiment/freeze/cohort/frozen/cohort.json"
upstream="$(cat "$source_experiment/controller.name")"
[[ "$upstream" =~ ^snpster-ancestry-v2-dense-split-[0-9_]+$ ]]
docker inspect "$upstream" >/dev/null
for active in $(docker ps --filter name=snpster-ancestry-v2- --format '{{.Names}}'); do
    if [[ "$active" != "$upstream" ]]; then
        echo "Refusing concurrent ancestry controller: $active" >&2; exit 1
    fi
done
sha256sum --check --status "$source_experiment/source_sha256.txt"
suffix="$(date +%Y%m%d_%H%M%S)"
experiment="$project/results_v3_dense_model_expansion_$suffix"
name="snpster-ancestry-v2-model98-$suffix"
mkdir "$experiment"
docker run --rm --network none --user "$(id -u):$(id -g)" \
    -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    --entrypoint python3 sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b \
    "$project/pedigree_benchmark/scripts/dense_model_expansion.py" freeze \
    --source-experiment "$source_experiment" --plan "$project/config/dense_model_expansion.json" \
    --output "$experiment/frozen"
for mode in smoke full; do
    mkdir -p "$experiment/$mode/launch"
    for folder in preparation founders truth raw standardized combined harmonized merged \
        imputed split qc indexed_qc preprocessing preprocessing_gate; do
        ln -s "$source_experiment/$mode/$folder" "$experiment/$mode/$folder"
    done
done
mkdir "$experiment/controller_home"
sha256sum "$project/dense_model_expansion.nf" "$project/run_dense_model_expansion.sh" \
    "$project/config/dense_model_expansion.json" "$project/production_matched_benchmark.nf" \
    "$project/production_matched.config" "$project/pedigree_benchmark/scripts/dense_model_expansion.py" \
    "$project/pedigree_benchmark/scripts/production_matched_models.py" \
    "$project/pedigree_benchmark/scripts/production_matched_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/verify_production_matched.py" \
    "$project/pedigree_benchmark/scripts/finalize_production_matched.py" \
    "$project/pedigree_benchmark/scripts/model_family_gnomix.py" \
    "$project/pedigree_benchmark/scripts/dense_split_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/pedigree_benchmark/scripts/screen_model_families.py" \
    "$project/pedigree_benchmark/scripts/realized_ancestry.py" \
    "$project/results_v3_global_precision_01/precision.jar" \
    "$experiment/frozen/cohort.json" "$experiment/frozen/candidate_plan.json" \
    "$HOME/.nextflow/framework/26.04.4/nextflow-26.04.4-one.jar" > "$experiment/source_sha256.txt"
cat "$source_experiment/instrumentation/assets.sha256" >> "$experiment/source_sha256.txt"
printf '%s\n' "$name" > "$experiment/controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e EXPERIMENT="$experiment" \
    -e SOURCE_EXPERIMENT="$source_experiment" -e UPSTREAM="$upstream" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" \
    -v "$source_experiment:$source_experiment:ro" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 60
        monitor="-v $SOURCE_EXPERIMENT/instrumentation:/opt/snpster-nextflow-metrics:ro -v $SOURCE_EXPERIMENT/instrumentation/ps:/usr/local/bin/ps:ro"
        options="--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v $SOURCE_EXPERIMENT:$SOURCE_EXPERIMENT:ro -v /srv/dependencies:/srv/dependencies:ro $monitor"
        echo "Frozen98 distinct settings; original ten outputs give eight unique baselines and two aliases;67 external SGDP people remain reserved."
        echo "Upstream progress: docker logs -f $UPSTREAM"
        echo "The Nextflow wait task uses docker wait; fitting cannot overlap the upstream analysis or bypass its gates."
        for mode in smoke full; do
            scope="$EXPERIMENT/$mode"
            if [[ "$mode" == smoke ]]; then smoke=true; else smoke=false; fi
            cd "$scope/launch"
            nextflow -log "$scope/nextflow.log" -C "$PROJECT/production_matched.config" \
                run "$PROJECT/dense_model_expansion.nf" -ansi-log true -qs 4 -work-dir "$scope/work" \
                -with-trace "$scope/trace.tsv" -with-report "$scope/report.html" \
                -with-timeline "$scope/timeline.html" -with-dag "$scope/pipeline.html" \
                --project "$PROJECT" --run_dir "$scope" --cohort "$EXPERIMENT/frozen/cohort.json" \
                --smoke "$smoke" --source_scope "$SOURCE_EXPERIMENT/$mode" \
                --source_experiment "$SOURCE_EXPERIMENT" --upstream_controller "$UPSTREAM" \
                --socket_gid "$(stat -c "%g" /var/run/docker.sock)" --container_options "$options"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$scope"
            session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
            [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
            nextflow clean -n "$session" > "$scope/cache_cleanup.plan"
            if ! grep -q "^Would remove " "$scope/cache_cleanup.plan" ||
                grep "^Would remove " "$scope/cache_cleanup.plan" | grep -v -F "Would remove $scope/work/" >/dev/null; then
                echo "Missing or unsafe exact-session cleanup plan" >&2; exit 1
            fi
            nextflow clean -f "$session" > "$scope/cache_cleanup.log"
            find "$scope/work" -depth -mindepth 1 -type d -empty -delete
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$scope" --post-cleanup
        done
        echo "Complete:98 verified configurations, all22 autosomes and both masks; native8 and EAS+SEA-merged7 TV; no production promotion."
    ' > "$experiment/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $experiment"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
