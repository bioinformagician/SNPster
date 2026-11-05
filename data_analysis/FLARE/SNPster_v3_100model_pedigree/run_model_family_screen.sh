#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${MODEL_FAMILY_SCREEN_OUTDIR:-results_v3_model_family_screen_01}"
if [[ ! "$output_name" =~ ^results_v3_model_family_screen_[a-zA-Z0-9_]+$ ]]; then
    echo "MODEL_FAMILY_SCREEN_OUTDIR must name an isolated results_v3_model_family_screen_* directory" >&2
    exit 1
fi
image="snpster-ancestry-benchmark:pedigree-100"
name="snpster-${output_name}-${1:-run}-${2:-stage}"
run_container() {
    local stage="$1" candidate="${2:-}" smoke="${3:-}"
    local args=(--stage "$stage" --outdir /results)
    if [[ "$stage" == recover ]]; then
        args+=(--previous "/workspace/$candidate")
    elif [[ -n "$candidate" ]]; then
        args+=(--candidate "$candidate")
    fi
    if [[ "$smoke" == smoke ]]; then args+=(--smoke); fi
    local command=(python3 /workspace/pedigree_benchmark/scripts/screen_model_families.py "${args[@]}")
    if [[ "$stage" == verify-inputs ]]; then
        command=(python3 -c 'from pathlib import Path; import sys; sys.path.insert(0, "/workspace/pedigree_benchmark/scripts"); from screen_model_families import verify; verify(Path("/results")); print("Verified all model-family input hashes")')
    fi
    docker run --rm --network none --name "$name" --cpus=6 --memory=40g \
        --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -v "$project_dir:/workspace:ro" \
        -v "$project_dir/$output_name:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$image" \
        "${command[@]}"
}
case "${1:-}" in
    --prepare)
        if [[ -e "$project_dir/$output_name" ]]; then
            echo "Output already exists; select a new model-family screen directory" >&2
            exit 1
        fi
        docker image inspect "$image" --format '{{.Id}}' >/dev/null
        mkdir "$project_dir/$output_name"
        docker image inspect "$image" --format '{{.Id}}' > "$project_dir/$output_name/flare_image.txt"
        run_container prepare 2>&1 | tee "$project_dir/$output_name/prepare.log"
        ;;
    --initialize)
        run_container initialize 2>&1 | tee "$project_dir/$output_name/initialize.log"
        ;;
    --recover)
        if [[ -e "$project_dir/$output_name" || ! "${2:-}" =~ ^results_v3_model_family_screen_[a-zA-Z0-9_]+$ ]]; then
            echo "Recovery requires a new output and an explicit prior screen directory" >&2
            exit 1
        fi
        mkdir "$project_dir/$output_name"
        run_container recover "$2" 2>&1 | tee "$project_dir/$output_name/recovery.log"
        ;;
    --matched-baseline)
        run_container matched-baseline 2>&1 | tee "$project_dir/$output_name/matched_baseline.log"
        ;;
    --flare)
        run_container flare "$2" "${3:-}" 2>&1 |
            tee "$project_dir/$output_name/${2}${3:+_smoke}.log"
        ;;
    --score)
        run_container score "$2" 2>&1 | tee "$project_dir/$output_name/${2}_score.log"
        ;;
    --summarize)
        run_container summarize 2>&1 | tee "$project_dir/$output_name/summary.log"
        ;;
    --verify-inputs)
        run_container verify-inputs
        ;;
    *)
        echo "Usage: bash run_model_family_screen.sh --prepare|--initialize|--matched-baseline|--flare CANDIDATE [smoke]|--score CANDIDATE|--summarize|--verify-inputs" >&2
        exit 1
        ;;
esac
