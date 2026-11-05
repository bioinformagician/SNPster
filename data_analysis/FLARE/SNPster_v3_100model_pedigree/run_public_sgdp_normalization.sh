#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
controller_image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
task_image=sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b
if [[ $# -eq 1 && "$1" == --smoke ]]; then
    mode=smoke
    smoke=true
    smoke_manifest=
elif [[ $# -eq 2 && "$1" == --run && -s "$2/smoke_passed.json" ]]; then
    mode=full
    smoke=false
    smoke_manifest="$(realpath "$2/smoke_passed.json")"
else
    echo "Usage: bash run_public_sgdp_normalization.sh --smoke | --run PASSED_SMOKE_DIRECTORY" >&2
    exit 1
fi
framework="$HOME/.nextflow/framework/26.04.4"
test -r "$framework/nextflow-26.04.4-one.jar"
test -S /var/run/docker.sock
docker image inspect "$controller_image" >/dev/null
docker image inspect "$task_image" >/dev/null
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_sgdp_normalization_${mode}_$suffix"
name="snpster-sgdp-normalization-${mode}-$suffix"
mkdir "$run"
mkdir -p "$run/launch" "$run/work" "$run/controller_home/framework/26.04.4"
chain="$repository/data pipeline/standardizer_module/dependencies/hg19ToHg38.over.chain.gz"
sha256sum "$framework/nextflow-26.04.4-one.jar" "$project/public_sgdp_normalization.nf" \
    "$project/public_sgdp_normalization.config" "$project/production_preflight.config" \
    "$project/run_public_sgdp_normalization.sh" "$project/config/sgdp_normalization.json" \
    "$project/requirements-sgdp-normalization.txt" "$project/Dockerfile.sgdp-normalization" \
    "$project/pedigree_benchmark/scripts/normalize_sgdp.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/finalize_sgdp_normalization.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$chain" > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
printf '%s\n' "$task_image" > "$run/task_image.txt"
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e RUN_DIR="$run" -e SMOKE="$smoke" \
    -e SMOKE_MANIFEST="$smoke_manifest" -e CHAIN="$chain" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$framework:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /bin/bash "$controller_image" -euo pipefail -c '
        helper="$PROJECT/pedigree_benchmark/scripts/finalize_sgdp_normalization.py"
        if [[ "$SMOKE" == false ]]; then
            python3 "$helper" require-smoke --manifest "$SMOKE_MANIFEST"
        fi
        nextflow -C "$PROJECT/public_sgdp_normalization.config" run "$PROJECT/public_sgdp_normalization.nf" \
            -ansi-log false -work-dir "$RUN_DIR/work" --project "$PROJECT" --run_dir "$RUN_DIR" \
            --input_manifest "$PROJECT/results_v3_public_sgdp_recovery_20261022_082600/summary/published_verification.json" \
            --chain "$CHAIN" --smoke "$SMOKE" \
            --container_options "--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro"
        if [[ "$SMOKE" == true ]]; then
            python3 "$helper" smoke --run "$RUN_DIR" --output "$RUN_DIR/smoke_passed.json"
        else
            python3 "$helper" full --run "$RUN_DIR" --output "$RUN_DIR/summary/pre_cleanup_verification.json"
        fi
        session="$(awk -F "\t" "END {print \$6}" "$RUN_DIR/launch/.nextflow/history")"
        if [[ ! "$session" =~ ^[a-f0-9-]{36}$ ]]; then
            echo "Unable to resolve exact successful Nextflow session" >&2
            exit 1
        fi
        nextflow clean -n "$session" > "$RUN_DIR/cache_cleanup.plan"
        if ! grep -q "^Would remove " "$RUN_DIR/cache_cleanup.plan" ||
            grep "^Would remove " "$RUN_DIR/cache_cleanup.plan" | grep -v -F "Would remove $RUN_DIR/work/" >/dev/null; then
            echo "Missing or unsafe exact-session cleanup plan" >&2
            exit 1
        fi
        echo "Verified published data. Cleaning only normalization session $session."
        nextflow clean -f "$session" > "$RUN_DIR/cache_cleanup.log"
        if find "$RUN_DIR/work" -mindepth 1 \( -type f -o -type l \) -print -quit | grep -q .; then
            echo "Normalization work files remain" >&2
            exit 1
        fi
        find "$RUN_DIR/work" -depth -mindepth 1 -type d -empty -delete
        if [[ "$SMOKE" == true ]]; then
            python3 "$helper" require-smoke --manifest "$RUN_DIR/smoke_passed.json"
        else
            python3 "$helper" full --run "$RUN_DIR" --output "$RUN_DIR/summary/post_cleanup_verification.json"
        fi
        du -sh "$RUN_DIR/work"
    ' > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
