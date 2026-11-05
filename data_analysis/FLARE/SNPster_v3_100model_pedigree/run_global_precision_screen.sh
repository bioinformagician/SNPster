#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${GLOBAL_PRECISION_OUTDIR:-results_v3_global_precision_01}"
if [[ ! "$output_name" =~ ^results_v3_global_precision_[a-zA-Z0-9_]+$ ]]; then
    echo "GLOBAL_PRECISION_OUTDIR must name an isolated results_v3_global_precision_* directory" >&2
    exit 1
fi
if [[ -e "$project_dir/$output_name" ]]; then
    echo "Output already exists; choose a new precision experiment directory" >&2
    exit 1
fi
jdk="$(dirname "$(dirname "$(readlink -f "$(command -v javac)")")")"
image="snpster-ancestry-benchmark:pedigree-100"
docker image inspect "$image" --format '{{.Id}}' >/dev/null
mkdir -p "$project_dir/$output_name/input"
docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/container_image.txt"
cp "$0" "$project_dir/$output_name/input/launcher.sh"
curl --fail --location --silent --show-error \
    https://raw.githubusercontent.com/browning-lab/flare/d6a6c8fd4479ebdb68529824acd11d2b6bbac4a1/src/admix/EstimatedGlobalAncProportions.java \
    -o "$project_dir/$output_name/input/EstimatedGlobalAncProportions.java"
for stage in prepare stock recompiled precision summarize; do
    docker run --rm --network none --name "snpster-${output_name}-${stage}" \
        --cpus=2 --memory=14g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -v "$jdk:/jdk:ro" -w /tmp "$image" \
        python3 /workspace/pedigree_benchmark/scripts/screen_global_precision.py \
            "$stage" --outdir /results \
        2>&1 | tee "$project_dir/$output_name/$stage.log"
done
