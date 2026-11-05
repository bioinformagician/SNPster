#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
if [[ $# -eq 2 && "$1" == --controller ]]; then
    run="$(realpath "$2")"
    if [[ "$run" != "$project"/results_v3_production_preflight_* || ! -d "$run/launch" ]]; then
        echo "Controller requires its dedicated preflight run directory." >&2
        exit 1
    fi
    cd "$run/launch"
    echo "Containerized Nextflow preflight: runtime identity and 22 reference audits."
    echo "This is a readiness audit, not an ancestry accuracy sweep."
    set +e
    nextflow -C "$project/production_preflight.config" run "$project/production_preflight.nf" \
        -ansi-log false -work-dir "$run/work" --project "$project" --run_dir "$run" \
        --container_options "--network none --user $(id -u):$(id -g) -v $repository:$repository:ro -v /srv/dependencies:/srv/dependencies:ro"
    status=$?
    set -e
    printf '%s\n' "$status" > "$run/controller.exitcode"
    echo "Nextflow preflight exited with status $status; inspect published summary before any scientific launch."
    exit "$status"
fi
if [[ $# -ne 1 || "$1" != --run ]]; then
    echo "Usage: bash run_production_preflight.sh --run" >&2
    exit 1
fi
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_production_preflight_$suffix"
name="snpster-ancestry-v2-preflight-$suffix"
framework="$HOME/.nextflow/framework/26.04.4"
if [[ ! -r "$framework/nextflow-26.04.4-one.jar" || ! -S /var/run/docker.sock ]]; then
    echo "Pinned Nextflow framework and accessible Docker daemon are required." >&2
    exit 1
fi
docker image inspect "$image" >/dev/null
docker image inspect sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d >/dev/null
mkdir "$run"
mkdir -p "$run/launch" "$run/work" "$run/controller_home/framework/26.04.4"
sha256sum "$framework/nextflow-26.04.4-one.jar" "$project/production_preflight.nf" \
    "$project/production_preflight.config" "$project/run_production_preflight.sh" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/config/ancestry_benchmark_contract_v2.json" > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -v "$repository:$repository:ro" -v "$run:$run:rw" \
    -v "$framework:$run/controller_home/framework/26.04.4:ro" \
    -v /srv/dependencies:/srv/dependencies:ro \
    -v /var/run/docker.sock:/var/run/docker.sock \
    --entrypoint /bin/bash "$image" "$project/run_production_preflight.sh" \
    --controller "$run" > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
