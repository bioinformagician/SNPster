#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
smoke="$project/results_v3_public_sgdp_smoke_20261021_203540"
import="$project/results_v3_public_sgdp_full_v2_20261021_204615"
if [[ $# -ne 1 || "$1" != --run ]]; then
    echo "Usage: bash run_public_sgdp_verify.sh --run" >&2
    exit 1
fi
test -s "$import/summary/full_import.json"
test "$(docker inspect --format '{{.State.Status}}:{{.State.ExitCode}}' snpster-ancestry-v2-sgdp-full-v2-20261021_204615)" = exited:0
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_public_sgdp_verify_$suffix"
name="snpster-ancestry-v2-sgdp-verify-$suffix"
framework="$HOME/.nextflow/framework/26.04.4"
test -r "$framework/nextflow-26.04.4-one.jar"
test -S /var/run/docker.sock
docker image inspect "$image" >/dev/null
docker image inspect sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d >/dev/null
mkdir "$run"
mkdir -p "$run/launch" "$run/work" "$run/controller_home/framework/26.04.4"
sha256sum "$framework/nextflow-26.04.4-one.jar" "$project/public_sgdp_verify.nf" \
    "$project/production_preflight.config" "$project/run_public_sgdp_verify.sh" \
    "$project/config/sgdp_public_import.json" \
    "$project/pedigree_benchmark/scripts/verify_public_sgdp.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$import/summary/full_import.json" > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e IMPORT_RUN="$import" -e SMOKE_RUN="$smoke" -e VERIFY_RUN="$run" \
    -e REPOSITORY="$repository" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v "$import:$import:rw" \
    -v "$framework:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /bin/bash "$image" -euo pipefail -c '
        nextflow -C "$PROJECT/production_preflight.config" run "$PROJECT/public_sgdp_verify.nf" \
            -ansi-log false -work-dir "$VERIFY_RUN/work" --project "$PROJECT" --run_dir "$VERIFY_RUN" \
            --import_dir "$IMPORT_RUN" --smoke_dir "$SMOKE_RUN" \
            --container_options "--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro"
        manifest="$VERIFY_RUN/summary/published_verification.json"
        test -s "$manifest"
        cd "$IMPORT_RUN/launch"
        session=2d2f3b41-7812-4261-aa63-91fe306fe82f
        nextflow clean -n "$session" > "$VERIFY_RUN/import_cleanup.plan"
        if grep "^Would remove " "$VERIFY_RUN/import_cleanup.plan" | grep -v -F "Would remove $IMPORT_RUN/work/" >/dev/null; then
            echo "Cleanup escaped the exact completed import work directory" >&2
            exit 1
        fi
        echo "Independent verification passed. Cleaning only import session $session."
        nextflow clean -f "$session" > "$VERIFY_RUN/import_cleanup.log"
        if find "$IMPORT_RUN/work" -mindepth 1 \( -type f -o -type l \) -print -quit | grep -q .; then
            echo "Unexpected import work files remain; cleanup incomplete" >&2
            exit 1
        fi
        find "$IMPORT_RUN/work" -depth -mindepth 1 -type d -empty -delete
        python3 "$PROJECT/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" recheck \
            --manifest "$manifest" --output "$VERIFY_RUN/summary/post_cleanup_verification.json"
        du -sh "$IMPORT_RUN/work"
    ' > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
