#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_production_cohort_$suffix"
name="snpster-ancestry-v2-cohort-$suffix"
mkdir -p "$run/launch" "$run/controller_home"
sha256sum "$project/production_cohort.nf" "$project/production_matched.config" \
    "$project/config/production_matched_protocol.json" "$project/run_production_cohort.sh" \
    "$project/pedigree_benchmark/scripts/freeze_production_cohort.py" > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
docker run --detach --name "$name" --init --network none --cpus=2 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e RUN_DIR="$run" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        nextflow -log "$RUN_DIR/nextflow.log" -C "$PROJECT/production_matched.config" \
            run "$PROJECT/production_cohort.nf" -ansi-log false -work-dir "$RUN_DIR/work" \
            -with-trace "$RUN_DIR/trace.tsv" -with-report "$RUN_DIR/report.html" \
            --project "$PROJECT" --run_dir "$RUN_DIR" \
            --container_options "--user $(id -u):$(id -g) -v $REPOSITORY:$REPOSITORY:ro -v /srv/dependencies:/srv/dependencies:ro"
    ' > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
