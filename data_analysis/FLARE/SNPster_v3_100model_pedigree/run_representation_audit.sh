#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${REPRESENTATION_AUDIT_OUTDIR:-results_v3_representation_audit_01}"
if [[ ! "$output_name" =~ ^results_v3_representation_audit_[a-zA-Z0-9_]+$ ]]; then
    echo "REPRESENTATION_AUDIT_OUTDIR must name an isolated results_v3_representation_audit_* directory" >&2
    exit 1
fi
if [[ -e "$project_dir/$output_name" ]]; then
    echo "Output already exists; choose a new audit directory" >&2
    exit 1
fi
mkdir "$project_dir/$output_name"
image="snpster-ancestry-benchmark:pedigree-100"
docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/container_image.txt"
for index in 0 1 2 3 4 5 6 7 8 9 summary; do
    args=(--summarize --outdir /results)
    if [[ "$index" != summary ]]; then
        printf -v folder 'rep_%02d' "$index"
        args=(--index "$index" --outdir "/results/$folder")
    fi
    docker run --rm --network none --name "snpster-${output_name}-${index}" \
        --cpus=2 --memory=8g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image" \
        python3 /workspace/pedigree_benchmark/scripts/audit_representations.py "${args[@]}" \
        2>&1 | tee "$project_dir/$output_name/$index.log"
done

docker run --rm --network none --name "snpster-${output_name}-collisions" \
    --cpus=2 --memory=8g --user "$(id -u):$(id -g)" \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/$output_name:/results:rw" \
    -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
    -w /tmp "$image" \
    python3 /workspace/pedigree_benchmark/scripts/audit_feature_collisions.py \
        --audit-dir "/workspace/$output_name" --outdir /results/collisions \
    2>&1 | tee "$project_dir/$output_name/collisions.log"
