#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${DECODER_SCREEN_OUTDIR:-results_v3_decoder_screen_01}"
if [[ ! "$output_name" =~ ^results_v3_decoder_screen_[a-zA-Z0-9_]+$ ]]; then
    echo "DECODER_SCREEN_OUTDIR must be an isolated results_v3_decoder_screen_* directory name" >&2
    exit 1
fi
if [[ -e "$project_dir/$output_name" ]]; then
    echo "Output already exists; choose a new DECODER_SCREEN_OUTDIR rather than overwriting experiments" >&2
    exit 1
fi
mkdir "$project_dir/$output_name"
image="snpster-ancestry-benchmark:pedigree-100"
for mode in audit ridge_pure ridge_mixed stochastic_pure stochastic_mixed extratrees_mixed summarize; do
    destination="/results/$mode"
    if [[ "$mode" == summarize ]]; then destination=/results; fi
    docker run --rm --network none \
        --name "snpster-${output_name}-${mode}" \
        --cpus=2 --memory=8g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image" \
        python3 /workspace/pedigree_benchmark/scripts/screen_v1_decoders.py \
            --mode "$mode" --outdir "$destination" \
        2>&1 | tee "$project_dir/$output_name/$mode.log"
done
