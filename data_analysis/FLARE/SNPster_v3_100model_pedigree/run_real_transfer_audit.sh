#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${REAL_TRANSFER_AUDIT_OUTDIR:-results_v3_real_transfer_audit_01}"
if [[ ! "$output_name" =~ ^results_v3_real_transfer_audit_[a-zA-Z0-9_]+$ ||
      -e "$project_dir/$output_name" ]]; then
    echo "Select a new isolated REAL_TRANSFER_AUDIT_OUTDIR" >&2
    exit 1
fi
image="snpster-ancestry-benchmark:gnomix-screen"
docker image inspect "$image" --format '{{.Id}}' >/dev/null
mkdir "$project_dir/$output_name"
docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/image.txt"
docker run --rm --network none --name "snpster-${output_name}" --cpus=6 --memory=40g \
    --user "$(id -u):$(id -g)" \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
    -e NUMBA_CACHE_DIR=/tmp/numba-cache -e MPLCONFIGDIR=/tmp/matplotlib \
    -e XDG_CACHE_HOME=/tmp/cache -e TMPDIR=/tmp \
    -v "$project_dir:/workspace:ro" \
    -v "$project_dir/results_v3_model_family_screen_05:/results:ro" \
    -v "$project_dir/$output_name:/audit:rw" \
    -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
    -v /srv/imputed/Lonni/325:/upstream:ro \
    -w /tmp "$image" \
    python3 /workspace/pedigree_benchmark/scripts/audit_real_transfer.py \
        --screen /results --outdir /audit --upstream /upstream \
    2>&1 | tee "$project_dir/$output_name/audit.log"
