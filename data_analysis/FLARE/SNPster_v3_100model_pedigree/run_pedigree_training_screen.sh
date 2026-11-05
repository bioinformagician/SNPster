#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${TRAINING_SCREEN_OUTDIR:-results_v3_pedigree_training_screen_01}"
if [[ ! "$output_name" =~ ^results_v3_pedigree_training_screen_[a-zA-Z0-9_]+$ ]]; then
    echo "TRAINING_SCREEN_OUTDIR must be an isolated results_v3_pedigree_training_screen_* name" >&2
    exit 1
fi
if [[ -e "$project_dir/$output_name" || -e "$project_dir/work_$output_name" ]]; then
    echo "Output/work already exists; choose a new TRAINING_SCREEN_OUTDIR" >&2
    exit 1
fi
mkdir "$project_dir/$output_name"
image="snpster-ancestry-benchmark:pedigree-100"
docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/container_image.txt"
fingerprint="$(find "$project_dir/pedigree_benchmark/scripts" -name '*.py' -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')"
printf '%s\n' "$fingerprint" > "$project_dir/$output_name/python_sources_fingerprint.txt"
sha256sum "$project_dir/pedigree_training.nf" "$project_dir/frozen_inputs.nf" \
    "$project_dir/main.nf" "$project_dir/config/params.json" "$project_dir/nextflow.config" \
    "$project_dir/config/frozen_training_v1.tsv" "$project_dir/run_pedigree_training_screen.sh" \
    > "$project_dir/$output_name/launch_sources.sha256"

docker run --rm --network none --name "snpster-${output_name}-cohort" \
    --cpus=20 --memory=180g --user "$(id -u):$(id -g)" \
    -e NXF_HOME=/workspace/.nextflow -e NXF_OFFLINE=true \
    -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    -v "$project_dir:/workspace:rw" \
    -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
    -w /workspace "$image" nextflow run pedigree_training.nf \
    -params-file config/params.json \
    -work-dir "/workspace/work_$output_name" \
    -with-trace "/workspace/$output_name/trace.txt" \
    -with-report "/workspace/$output_name/report.html" \
    -with-timeline "/workspace/$output_name/timeline.html" \
    --outdir "/workspace/$output_name/cohort" \
    --frozen_candidates_manifest /workspace/config/frozen_training_v1.tsv \
    --pedigree_reuse_dir /workspace/results_v3_sweep_100/pedigree/reusable_inputs \
    --pedigree_design full --pedigree_families 108 --pedigree_batch_size 10 \
    --pedigree_seed 20261014 --pedigree_keep_high_density false \
    --pedigree_fingerprint "$fingerprint" \
    2>&1 | tee "$project_dir/$output_name/cohort.log"

sha256sum --check "$project_dir/$output_name/launch_sources.sha256"
current_fingerprint="$(find "$project_dir/pedigree_benchmark/scripts" -name '*.py' -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')"
if [[ "$current_fingerprint" != "$fingerprint" ]]; then
    echo "Python experiment sources changed during cohort generation; refusing mixed-version scoring" >&2
    exit 1
fi

for mode in baseline ridge_pure ridge_pedigree stochastic_pure stochastic_pedigree extratrees_pure extratrees_pedigree summarize; do
    destination="/results/$mode"
    if [[ "$mode" == summarize ]]; then destination=/results; fi
    docker run --rm --network none --name "snpster-${output_name}-${mode}" \
        --cpus=2 --memory=8g --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e SKLEARN_N_JOBS=2 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image" \
        python3 /workspace/pedigree_benchmark/scripts/screen_training_decoders.py \
            --mode "$mode" --outdir "$destination" \
            --training-dir "/workspace/$output_name/cohort" \
        2>&1 | tee "$project_dir/$output_name/$mode.log"
done
