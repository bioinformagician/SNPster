#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
if [[ $# -ne 1 ]]; then echo "Usage: bash resume_production_matched_benchmark.sh ABSOLUTE_EXPERIMENT" >&2; exit 1; fi
run="$(realpath "$1")"
[[ "$run" == "$project"/results_v3_production_matched_* ||
   "$run" == "$project"/results_v3_observed_truth_support_*/benchmark ]]
test -s "$run/source_sha256.txt"
fresh=1
if [[ -f "$run/controller.name" ]]; then
    fresh=0
    previous_controller="$(cat "$run/controller.name")"
    if [[ "$(docker inspect --format '{{.State.Running}}' "$previous_controller")" == true ]]; then
        echo "Refusing to resume an experiment with a running controller: $previous_controller" >&2; exit 1
    fi
elif [[ -s "$run/smoke/launch/.nextflow/history" ]]; then
    echo "Existing experiment has no controller tracking; refusing an untracked resume." >&2; exit 1
fi
cohort="${SNPSTER_MATCHED_COHORT:-$project/results_v3_production_cohort_20261022_110441/cohort/frozen/cohort.json}"
if ! cmp -s "$cohort" "$run/frozen_cohort.json"; then
    echo "Requested cohort differs from this experiment; provide its original SNPSTER_MATCHED_COHORT." >&2; exit 1
fi
suffix="$(date +%Y%m%d_%H%M%S)"
name="snpster-ancestry-v2-production-recovery-$suffix"
if [[ "$fresh" == 1 ]]; then name="snpster-ancestry-v2-production-$suffix"; fi
attempt="$run/attempts/$suffix"
mkdir -p "$attempt"
cp "$run/source_sha256.txt" "$attempt/previous_source_sha256.txt"
bash "$project/prepare_nextflow_metrics.sh" "$run/instrumentation"
bash "$project/prepare_production_executables.sh" "$run/production_executables"
while read -r previous filename; do sha256sum "$filename"; done < "$attempt/previous_source_sha256.txt" > "$run/source_sha256.txt"
sha256sum "$project/resume_production_matched_benchmark.sh" "$project/prepare_nextflow_metrics.sh" \
    "$project/pedigree_benchmark/scripts/nextflow_ps.sh" "$project/production_runtime_smoke.nf" \
    "$project/prepare_production_executables.sh" "$run/production_executables/assets.sha256" \
    "$run/instrumentation/assets.sha256" >> "$run/source_sha256.txt"
cat "$run/instrumentation/assets.sha256" >> "$run/source_sha256.txt"
cat "$run/production_executables/assets.sha256" >> "$run/source_sha256.txt"
if [[ -d "$run/failures" ]]; then
    find "$run/failures" -type f -print0 | sort -z | xargs -0 -r sha256sum >> "$run/source_sha256.txt"
fi
printf '%s\n' "$name" > "$attempt/controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e EXPERIMENT="$run" -e COHORT="$cohort" \
    -e ATTEMPT="$attempt" \
    -e COLUMNS=180 -e LINES=60 \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 60
        monitor="-v $EXPERIMENT/instrumentation:/opt/snpster-nextflow-metrics:ro -v $EXPERIMENT/instrumentation/ps:/usr/local/bin/ps:ro"
        options="--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro $monitor"
        harmonizer="-v $EXPERIMENT/production_executables/plink:/app/dependencies/plink:ro -v $EXPERIMENT/production_executables/plink2:/app/dependencies/plink2:ro"
        cleanup_session() {
            local run="$1" prefix="$2" session
            session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
            [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
            nextflow clean -n "$session" > "$prefix.cleanup.plan"
            if ! grep -q "^Would remove " "$prefix.cleanup.plan" ||
                grep "^Would remove " "$prefix.cleanup.plan" | grep -v -F "Would remove $run/work/" >/dev/null; then
                echo "Missing or unsafe exact-session cleanup plan" >&2; exit 1
            fi
            nextflow clean -f "$session" > "$prefix.cleanup.log"
            find "$run/work" -depth -mindepth 1 -type d -empty -delete
        }
        runtime="$ATTEMPT/runtime"
        mkdir -p "$runtime/launch"
        cd "$runtime/launch"
        echo "Testing required subprocess tools and real resource monitoring in all10 original pinned task images."
        nextflow -log "$runtime/nextflow.log" -C "$PROJECT/production_matched.config" \
            run "$PROJECT/production_runtime_smoke.nf" -ansi-log true -work-dir "$runtime/work" \
            -with-trace "$runtime/trace.tsv" -with-report "$runtime/report.html" \
            -with-timeline "$runtime/timeline.html" -with-dag "$runtime/pipeline.html" \
            --project "$PROJECT" --run_dir "$runtime" --container_options "$options" --harmonizer_options "$harmonizer"
        awk -F "\t" "NR>1 && (\$5!=\"COMPLETED\" || \$6!=0) {bad=1} END {if(NR!=11 || bad) exit 1}" "$runtime/trace.tsv"
        [[ "$(find "$runtime/probes" -maxdepth 1 -type f -name "*.txt" | wc -l)" == 10 ]]
        sha256sum "$runtime/probes/"*.txt "$runtime/trace.tsv" "$runtime/report.html" \
            "$runtime/timeline.html" "$runtime/pipeline.html" > "$runtime/published.sha256"
        sha256sum --check "$runtime/published.sha256"
        cleanup_session "$runtime" "$runtime/cache"
        sha256sum --check "$runtime/published.sha256"
        echo "All10 original task images passed subprocess and metrics checks; probe cache cleaned and evidence reverified."
        for mode in smoke full; do
            run="$EXPERIMENT/$mode"
            if [[ -s "$run/post_cleanup_verification.json" ]]; then
                python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$run" --post-cleanup
                continue
            fi
            if [[ "$mode" == smoke ]]; then smoke=true; else smoke=false; fi
            cd "$run/launch"
            resume=()
            if [[ -s .nextflow/history ]]; then
                session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
                [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
                resume=(-resume "$session")
            fi
            echo "Starting/resuming $mode; original images, scientific parameters and input genotypes unchanged."
            nextflow -log "$ATTEMPT/$mode.nextflow.log" -C "$PROJECT/production_matched.config" \
                run "$PROJECT/production_matched_benchmark.nf" "${resume[@]}" -ansi-log true -work-dir "$run/work" \
                -with-trace "$ATTEMPT/$mode.trace.tsv" -with-report "$ATTEMPT/$mode.report.html" \
                -with-timeline "$ATTEMPT/$mode.timeline.html" -with-dag "$ATTEMPT/$mode.pipeline.html" \
                --project "$PROJECT" --run_dir "$run" --cohort "$COHORT" --smoke "$smoke" \
                --container_options "$options" --harmonizer_options "$harmonizer"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$run"
            cleanup_session "$run" "$ATTEMPT/$mode"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$run" --post-cleanup
            du -sh "$run/work"
        done
        echo "Full frozen ten-model production-matched development comparison independently verified."
    ' > "$attempt/controller.id"
printf '%s\n' "$name" > "$run/controller.name"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
