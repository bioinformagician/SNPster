#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
if [[ $# -ne 1 ]]; then
    echo "Usage: bash run_dense_split_sweep.sh ABSOLUTE_FROZEN_EXPERIMENT" >&2; exit 1
fi
experiment="$(realpath "$1")"
[[ "$experiment" == "$project"/results_v3_dense_split_* ]]
cohort="$experiment/freeze/cohort/frozen/cohort.json"
test -s "$cohort"
if [[ -e "$experiment/controller.name" ]]; then
    echo "Existing sweep controller recorded; inspect its evidence before an explicit resume." >&2; exit 1
fi
active="$(docker ps --filter name=snpster-ancestry-v2- --format '{{.Names}}')"
if [[ -n "$active" ]]; then
    printf 'Refusing a concurrent ancestry controller:\n%s\n' "$active" >&2; exit 1
fi
suffix="$(date +%Y%m%d_%H%M%S)"
name="snpster-ancestry-v2-dense-split-$suffix"
mkdir -p "$experiment/runtime/launch" "$experiment/references_smoke/launch" \
    "$experiment/references_full/launch" "$experiment/references" \
    "$experiment/smoke/launch" "$experiment/full/launch"
ln -s /srv/dependencies/imputation_runner/imputer/plink.GRCh38.map "$experiment/references/plink.GRCh38.map"
bash "$project/prepare_nextflow_metrics.sh" "$experiment/instrumentation"
bash "$project/prepare_production_executables.sh" "$experiment/production_executables"
sha256sum "$project/production_matched_benchmark.nf" "$project/production_cohort.nf" \
    "$project/production_matched.config" "$project/production_runtime_smoke.nf" \
    "$project/run_dense_split_sweep.sh" "$project/config/dense_split_protocol.json" \
    "$project/config/production_matched_protocol.json" "$cohort" \
    "$project/pedigree_benchmark/scripts/freeze_dense_split.py" \
    "$project/pedigree_benchmark/scripts/dense_split_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/production_matched_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/production_matched_models.py" \
    "$project/pedigree_benchmark/scripts/verify_production_matched.py" \
    "$project/pedigree_benchmark/scripts/finalize_production_matched.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/pedigree_benchmark/scripts/model_family_gnomix.py" \
    "$project/pedigree_benchmark/scripts/realized_ancestry.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/pedigree_benchmark/scripts/screen_model_families.py" \
    "$project/prepare_nextflow_metrics.sh" "$project/prepare_production_executables.sh" \
    "$project/pedigree_benchmark/scripts/nextflow_ps.sh" \
    "$HOME/.nextflow/framework/26.04.4/nextflow-26.04.4-one.jar" > "$experiment/source_sha256.txt"
cat "$experiment/instrumentation/assets.sha256" "$experiment/production_executables/assets.sha256" \
    >> "$experiment/source_sha256.txt"
printf '%s\n' "$name" > "$experiment/controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$experiment/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e EXPERIMENT="$experiment" -e COHORT="$cohort" \
    -v "$repository:$repository:ro" -v "$experiment:$experiment:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$experiment/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$experiment" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 60
        monitor="-v $EXPERIMENT/instrumentation:/opt/snpster-nextflow-metrics:ro -v $EXPERIMENT/instrumentation/ps:/usr/local/bin/ps:ro"
        options="--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro $monitor"
        harmonizer="-v $EXPERIMENT/production_executables/plink:/app/dependencies/plink:ro -v $EXPERIMENT/production_executables/plink2:/app/dependencies/plink2:ro"
        cleanup_session() {
            local scope="$1" session
            [[ "$scope" == "$EXPERIMENT/"* ]]
            session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
            [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
            nextflow clean -n "$session" > "$scope/cache_cleanup.plan"
            if ! grep -q "^Would remove " "$scope/cache_cleanup.plan" ||
                grep "^Would remove " "$scope/cache_cleanup.plan" | grep -v -F "Would remove $scope/work/" >/dev/null; then
                echo "Missing or unsafe exact-session cleanup plan" >&2; exit 1
            fi
            nextflow clean -f "$session" > "$scope/cache_cleanup.log"
            find "$scope/work" -depth -mindepth 1 -type d -empty -delete
        }
        run_nf() {
            local scope="$1" workflow="$2"
            shift 2
            cd "$scope/launch"
            nextflow -log "$scope/nextflow.log" -C "$PROJECT/production_matched.config" \
                run "$PROJECT/$workflow" -ansi-log true -qs 4 -work-dir "$scope/work" \
                -with-trace "$scope/trace.tsv" -with-report "$scope/report.html" \
                -with-timeline "$scope/timeline.html" -with-dag "$scope/pipeline.html" \
                --project "$PROJECT" --run_dir "$scope" --container_options "$options" "$@"
        }
        echo "Approved synthetic development split; production resources and application methods remain unchanged."
        echo "Reference/training counts after identity-family screening:"
        python3 -c "import json,os; c=json.load(open(os.environ[\"COHORT\"])); print(c[\"retained_counts\"]); print(\"Reserved external observed-data people:\",len(c[\"reserved_external_ids\"]))"
        run_nf "$EXPERIMENT/runtime" production_runtime_smoke.nf --harmonizer_options "$harmonizer"
        awk -F "\t" "NR>1 && (\$5!=\"COMPLETED\" || \$6!=0) {bad=1} END {if(NR!=11 || bad) exit 1}" "$EXPERIMENT/runtime/trace.tsv"
        sha256sum "$EXPERIMENT/runtime/probes/"*.txt > "$EXPERIMENT/runtime/published.sha256"
        sha256sum --check "$EXPERIMENT/runtime/published.sha256"
        cleanup_session "$EXPERIMENT/runtime"
        sha256sum --check "$EXPERIMENT/runtime/published.sha256"
        for mode in smoke full; do
            references="$EXPERIMENT/references_$mode"
            if [[ "$mode" == smoke ]]; then
                smoke=true
                chromosomes=(22)
            else
                smoke=false
                chromosomes=({1..21})
                echo "chr22 smoke passed every genotype, inheritance, model and scoring check; extending to all22 autosomes and both masks."
                available="$(df -PB1 "$EXPERIMENT" | awk "NR==2 {print \$4}")"
                required="$(python3 -c "import json,os; from pathlib import Path; p=Path(os.environ[\"EXPERIMENT\"]); c=json.load(open(os.environ[\"COHORT\"])); n=sum(x.stat().st_size for x in (p/\"smoke\").rglob(\"*\") if x.is_file()); spans=c[\"chromosome_spans\"]; ratio=sum(b-a for a,b in spans.values())/(spans[\"22\"][1]-spans[\"22\"][0]); refs=sum(Path(x).stat().st_size for x in c[\"template_sha256\"] if Path(x).suffix==\".gz\"); print(int(n*ratio*4+refs+20*1024**3))")"
                if (( available < required )); then
                    echo "Full-scope capacity gate failed: need $required bytes, available $available. Verified smoke and resumable evidence retained." >&2; exit 1
                fi
            fi
            echo "Preparing donor-purged benchmark reference copies ($mode); production references are read-only."
            chromosome_csv="$(IFS=,; echo "${chromosomes[*]}")"
            run_nf "$references" production_cohort.nf --prepare_references true \
                --cohort "$COHORT" --chromosomes "$chromosome_csv" --reference_dir "$EXPERIMENT/references"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" \
                --run "$references" --references "$EXPERIMENT/references" --cohort "$COHORT" \
                --chromosomes "${chromosomes[@]}"
            cleanup_session "$references"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$references" --post-cleanup
            echo "Starting $mode production preprocessing and ten-model sweep; no model runs before all scoped genotype gates pass."
            run_nf "$EXPERIMENT/$mode" production_matched_benchmark.nf --cohort "$COHORT" --smoke "$smoke" \
                --preprocessing_script "$PROJECT/pedigree_benchmark/scripts/dense_split_preprocessing.py" \
                --imputation_root "$EXPERIMENT/references" --harmonizer_options "$harmonizer"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$EXPERIMENT/$mode"
            cleanup_session "$EXPERIMENT/$mode"
            python3 "$PROJECT/pedigree_benchmark/scripts/finalize_production_matched.py" --run "$EXPERIMENT/$mode" --post-cleanup
        done
        cd "$EXPERIMENT/freeze/launch"
        find "$EXPERIMENT/freeze/cohort" -type f -print0 | sort -z | xargs -0 sha256sum > "$EXPERIMENT/freeze/published.sha256"
        sha256sum --check "$EXPERIMENT/freeze/published.sha256"
        cleanup_session "$EXPERIMENT/freeze"
        sha256sum --check "$EXPERIMENT/freeze/published.sha256"
        echo "Complete: verified synthetic development comparison,15 evaluation controls, both masks, all22 autosomes, native8 and merged7 TV error. No production promotion;67 external observed-data people remain reserved."
    ' > "$experiment/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $experiment"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
