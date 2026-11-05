#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
image=sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1
if [[ $# -ne 1 || "$1" != --smoke ]]; then
    echo "Usage: bash run_public_sgdp_import.sh --smoke" >&2
    exit 1
fi
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_public_sgdp_smoke_$suffix"
name="snpster-ancestry-v2-sgdp-$suffix"
framework="$HOME/.nextflow/framework/26.04.4"
if [[ ! -r "$framework/nextflow-26.04.4-one.jar" || ! -S /var/run/docker.sock ]]; then
    echo "Pinned Nextflow framework and accessible Docker daemon required." >&2
    exit 1
fi
docker image inspect "$image" >/dev/null
docker image inspect sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d >/dev/null
mkdir "$run"
mkdir -p "$run/launch" "$run/work" "$run/controller_home/framework/26.04.4"
sha256sum "$framework/nextflow-26.04.4-one.jar" "$project/public_sgdp_import.nf" \
    "$project/public_sgdp_import.config" "$project/production_preflight.config" \
    "$project/run_public_sgdp_import.sh" "$project/config/sgdp_public_import.json" \
    "$project/pedigree_benchmark/scripts/import_sgdp_public.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
docker run --detach --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -v "$repository:$repository:ro" -v "$run:$run:rw" \
    -v "$framework:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run/launch" \
    --entrypoint /usr/local/bin/nextflow "$image" \
    -C "$project/public_sgdp_import.config" run "$project/public_sgdp_import.nf" \
    -ansi-log false -work-dir "$run/work" --project "$project" --run_dir "$run" \
    --container_options "--user $(id -u):$(id -g) -v $repository:$repository:ro" \
    > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}} {{.Id}}' "$name"
