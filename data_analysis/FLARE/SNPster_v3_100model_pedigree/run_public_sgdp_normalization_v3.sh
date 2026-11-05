#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
run="$project/results_v3_sgdp_normalization_full_v2_20261022_095828"
smoke="$project/results_v3_sgdp_normalization_smoke_20261022_095448/smoke_passed.json"
if [[ $# -ne 1 || "$1" != --resume ]]; then
    echo "Usage: bash run_public_sgdp_normalization_v3.sh --resume" >&2
    exit 1
fi
test -s "$smoke"
test -s "$run/launch/.nextflow/history"
test "$(docker inspect --format '{{.State.Status}}:{{.State.ExitCode}}' snpster-sgdp-normalization-v2-full-20261022_095828)" = exited:1
resume="$(awk -F '\t' 'END {print $6}' "$run/launch/.nextflow/history")"
[[ "$resume" =~ ^[a-f0-9-]{36}$ ]]
suffix="$(date +%Y%m%d_%H%M%S)"
name="snpster-sgdp-normalization-v3-full-$suffix"
attempt="$run/attempts/$suffix"
mkdir "$attempt"
cp "$run/source_sha256.txt" "$attempt/previous_source_sha256.txt"
sha256sum "$project/public_sgdp_normalization_v3.nf" \
    "$project/run_public_sgdp_normalization_v3.sh" > "$attempt/recovery_source_sha256.txt"
cat "$attempt/previous_source_sha256.txt" "$attempt/recovery_source_sha256.txt" > "$run/source_sha256.txt"
framework="$HOME/.nextflow/framework/26.04.4"
test -r "$framework/nextflow-26.04.4-one.jar"
printf '%s\n' "$name" > "$attempt/controller.name"
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e RUN_DIR="$run" \
    -e SMOKE_MANIFEST="$smoke" -e ATTEMPT="$attempt" -e RESUME="$resume" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$framework:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /bin/bash "$image" -euo pipefail -c '
        helper="$PROJECT/pedigree_benchmark/scripts/finalize_sgdp_normalization.py"
        python3 "$helper" require-smoke --manifest "$SMOKE_MANIFEST"
        nextflow -log "$ATTEMPT/nextflow.log" -C "$PROJECT/public_sgdp_normalization.config" \
            run "$PROJECT/public_sgdp_normalization_v3.nf" -resume "$RESUME" \
            -ansi-log false -work-dir "$RUN_DIR/work" --project "$PROJECT" --run_dir "$RUN_DIR" \
            -with-trace "$ATTEMPT/trace.tsv" -with-report "$ATTEMPT/report.html" -with-timeline "$ATTEMPT/timeline.html" \
            --input_manifest "$PROJECT/results_v3_public_sgdp_recovery_20261022_082600/summary/published_verification.json" \
            --chain "$REPOSITORY/data pipeline/standardizer_module/dependencies/hg19ToHg38.over.chain.gz" \
            --smoke false \
            --container_options "--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro"
        python3 "$helper" full --run "$RUN_DIR" --output "$RUN_DIR/summary/pre_cleanup_verification.json"
        nextflow clean -n "$RESUME" > "$ATTEMPT/cache_cleanup.plan"
        if ! grep -q "^Would remove " "$ATTEMPT/cache_cleanup.plan" ||
            grep "^Would remove " "$ATTEMPT/cache_cleanup.plan" | grep -v -F "Would remove $RUN_DIR/work/" >/dev/null; then
            echo "Missing or unsafe exact-session cleanup plan" >&2; exit 1
        fi
        echo "Verified all22 normalization and genetic QC. Cleaning only session $RESUME."
        nextflow clean -f "$RESUME" > "$ATTEMPT/cache_cleanup.log"
        if find "$RUN_DIR/work" -mindepth 1 \( -type f -o -type l \) -print -quit | grep -q .; then
            echo "Retaining unindexed failed-task evidence in work." >&2
        fi
        find "$RUN_DIR/work" -depth -mindepth 1 -type d -empty -delete
        python3 "$helper" full --run "$RUN_DIR" --output "$RUN_DIR/summary/post_cleanup_verification.json"
        du -sh "$RUN_DIR/work"
    ' > "$attempt/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
