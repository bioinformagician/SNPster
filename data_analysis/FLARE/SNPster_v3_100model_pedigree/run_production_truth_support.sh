#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
cohort="$project/results_v3_production_cohort_20261022_110441/cohort/frozen/cohort.json"
suffix="$(date +%Y%m%d_%H%M%S)"
if [[ $# -gt 1 ]]; then echo "Usage: bash run_production_truth_support.sh [EXISTING_FAILED_RUN]" >&2; exit 1; fi
resume=0
if [[ $# -eq 1 ]]; then
    run="$(realpath "$1")"
    [[ "$run" == "$project"/results_v3_observed_truth_support_* ]]
    test -s "$run/launch/.nextflow/history"
    previous="$(cat "$run/controller.name")"
    if [[ "$(docker inspect --format '{{.State.Running}}' "$previous")" == true ]]; then
        echo "Refusing to resume running controller: $previous" >&2; exit 1
    fi
    resume=1
    name="snpster-ancestry-v2-observed-support-recovery-$suffix"
else
    run="$project/results_v3_observed_truth_support_$suffix"
    name="snpster-ancestry-v2-observed-support-$suffix"
fi
attempt="$run/attempts/$suffix"
mkdir -p "$run/launch" "$run/controller_home"
mkdir -p "$attempt"
if [[ "$resume" == 1 ]]; then
    cp "$run/source_sha256.txt" "$attempt/previous_source_sha256.txt"
fi
sha256sum "$project/production_truth_support.nf" "$project/production_matched.config" \
    "$project/pedigree_benchmark/scripts/prepare_production_truth_support.py" \
    "$project/pedigree_benchmark/scripts/production_matched_preprocessing.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/run_production_truth_support.sh" "$cohort" > "$run/source_sha256.txt"
if [[ -d "$run/failures" ]]; then
    find "$run/failures" -type f -print0 | sort -z | xargs -0 -r sha256sum >> "$run/source_sha256.txt"
fi
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 -e HOME="$HOME" \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e RUN="$run" -e COHORT="$cohort" \
    -e RESUME="$resume" -e ATTEMPT="$attempt" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        resume=()
        if [[ "$RESUME" == 1 ]]; then
            session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
            [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
            resume=(-resume "$session")
        fi
        echo "Recovering original observed assay support; score>=0.8, genotype completeness>=90% after assay QC, typed markers>=100/mask."
        echo "Original missing and low-quality calls remain excluded; no imputed evaluation truth."
        nextflow -log "$ATTEMPT/nextflow.log" -C "$PROJECT/production_matched.config" \
            run "$PROJECT/production_truth_support.nf" "${resume[@]}" -ansi-log false -work-dir "$RUN/work" \
            -with-trace "$ATTEMPT/trace.tsv" -with-report "$ATTEMPT/report.html" \
            --project "$PROJECT" --run_dir "$RUN" --cohort "$COHORT" \
            --container_options "--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro"
        python3 -c "import json,sys;from pathlib import Path;sys.path.insert(0,str(Path(\"$PROJECT\")/\"pedigree_benchmark/scripts\"));from verify_sgdp_full_manifest import require_hashes,source_hashes; p=Path(\"$RUN\");r=json.loads((p/\"frozen/frozen/complete.json\").read_text());assert r[\"status\"]==\"all22_observed_truth_support_verified\";require_hashes(r[\"published_sha256\"]);require_hashes({str(p/\"frozen/frozen/cohort.json\"):r[\"cohort_sha256\"]});require_hashes(source_hashes(p/\"source_sha256.txt\"));print(\"All22 durable observed truth files and source hashes verified.\",flush=True)"
        session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
        [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
        nextflow clean -n "$session" > "$RUN/cache_cleanup.plan"
        if ! grep -q "^Would remove " "$RUN/cache_cleanup.plan" ||
            grep "^Would remove " "$RUN/cache_cleanup.plan" | grep -v -F "Would remove $RUN/work/" >/dev/null; then
            echo "Unsafe observed-support cleanup plan" >&2; exit 1
        fi
        nextflow clean -f "$session" > "$RUN/cache_cleanup.log"
        find "$RUN/work" -depth -mindepth 1 -type d -empty -delete
        python3 -c "import json,sys;from pathlib import Path;sys.path.insert(0,str(Path(\"$PROJECT\")/\"pedigree_benchmark/scripts\"));from verify_sgdp_full_manifest import require_hashes,source_hashes;from audit_production_ancestry import write_json;p=Path(\"$RUN\");r=json.loads((p/\"frozen/frozen/complete.json\").read_text());require_hashes(r[\"published_sha256\"]);require_hashes({str(p/\"frozen/frozen/cohort.json\"):r[\"cohort_sha256\"]});require_hashes(source_hashes(p/\"source_sha256.txt\"));write_json(p/\"post_cleanup_verification.json\",dict(status=\"all_durable_observed_truth_verified_after_cleanup\",files_verified=len(r[\"published_sha256\"])));print(\"Observed support reverified after exact successful-cache cleanup.\",flush=True)"
        SNPSTER_MATCHED_COHORT="$RUN/frozen/frozen/cohort.json" SNPSTER_MATCHED_RUN="$RUN/benchmark" \
            bash "$PROJECT/run_production_matched_benchmark.sh"
        child="$(cat "$RUN/benchmark/controller.name")"
        echo "Following production preprocessing and ten-model controller: $child"
        docker logs -f "$child"
        status="$(docker wait "$child")"
        docker inspect "$child" > "$RUN/benchmark/controller.inspect.json"
        docker logs "$child" > "$RUN/benchmark/controller.log" 2>&1
        [[ "$status" == 0 ]]
    ' > "$attempt/controller.id"
printf '%s\n' "$name" > "$run/controller.name"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
