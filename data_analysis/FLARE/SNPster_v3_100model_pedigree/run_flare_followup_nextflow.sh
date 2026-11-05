#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_name="${FLARE_FOLLOWUP_OUTDIR:-results_v3_flare_nextflow_01}"
if [[ ! "$output_name" =~ ^results_v3_flare_nextflow_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run && "$1" != --resume) ]]; then
    echo "Usage: bash run_flare_followup_nextflow.sh --prepare|--smoke|--run|--resume" >&2
    exit 1
fi
output="$project/$output_name"
image="sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d"
mounts=(-v "$project:/workspace:ro"
    -v "$project/results_v3_observed_dense_01:/benchmark:ro"
    -v "$project/results_v3_flare_reference_density_01:/screen:ro"
    -v "$project/results_v3_model_family_screen_05:/results:ro"
    -v /srv/dependencies/ancestry:/data:ro)
if [[ "$1" == --prepare ]]; then
    mkdir "$output"
    printf '%s\n' "$image" > "$output/image.txt"
    docker run --rm --network none --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 \
        "${mounts[@]}" -v "$output:/followup:rw" "$image" python3 -c '
import sys
from pathlib import Path
sys.path.insert(0,"/workspace/pedigree_benchmark/scripts")
from pilot_matched_imputation import fingerprint_files
from model_family_gnomix import write_json
project=Path("/workspace")
files=[project / name for name in ("flare_followup.nf","flare_followup.config","run_flare_followup_nextflow.sh",
    "config/flare_nextflow_followup.json","pedigree_benchmark/scripts/flare_nextflow_tasks.py",
    "pedigree_benchmark/scripts/finalize_flare_nextflow.py","tests/test_flare_nextflow.py")]
write_json("/followup/source_sha256.json",fingerprint_files(files))
'
    NXF_OFFLINE=true nextflow -version > "$output/nextflow_version.txt"
    exit
fi
if [[ ! -f "$output/source_sha256.json" ]]; then
    echo "Freeze sources using --prepare first." >&2
    exit 1
fi
mode=full
[[ "$1" == --smoke ]] && mode=smoke
if [[ "$mode" == full && ! -f "$output/smoke/published_verification.json" ]]; then
    echo "Verified published smoke results are required." >&2
    exit 1
fi
run="$output/$mode"
mkdir -p "$run/launch"
if [[ -f "$run/cache_cleanup.complete" ]]; then
    echo "This completed run was verified and its cache removed; choose a new output version." >&2
    exit 1
fi
cd "$run/launch"
if [[ -f "$run/run_name.txt" && "$1" != --resume && "$mode" == full ]]; then
    echo "Existing run requires explicit --resume." >&2
    exit 1
fi
name="flare_followup_${mode}_$(date +%Y%m%d_%H%M%S)"
printf '%s\n' "$name" > "$run/run_name.txt"
options="--network none --cpus=6 --memory=12g --user $(id -u):$(id -g)"
for mount in "${mounts[@]}"; do
    options+=" $mount"
done
resume=()
[[ "$1" == --resume || -f .nextflow/history ]] && resume=(-resume)
export NXF_OFFLINE=true NXF_OPTS="-Xms256m -Xmx2g"
nextflow -C "$project/flare_followup.config" run "$project/flare_followup.nf" \
    -name "$name" "${resume[@]}" -w "$run/work" \
    --mode "$mode" --outdir "$output" \
    --source_manifest "$output/source_sha256.json" --container_options "$options"
docker run --rm --network none --cpus=2 --memory=6g --user "$(id -u):$(id -g)" \
    -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
    "${mounts[@]}" -v "$output:$output:rw" "$image" \
    python3 /workspace/pedigree_benchmark/scripts/finalize_flare_nextflow.py --run-dir "$run"
if [[ "${FLARE_KEEP_CACHE_FOR_RESUME_TEST:-0}" == 1 ]]; then
    echo "Keeping this verified smoke cache temporarily for a bounded resume test." >&2
    exit
fi
if [[ ! -f "$run/published_verification.json" || "$PWD" != "$run/launch" ]]; then
    echo "Refusing cache cleanup without independent published-output verification and exact launch directory." >&2
    exit 1
fi
nextflow clean -n "$name" > "$run/cache_cleanup.plan"
nextflow clean -f "$name" > "$run/cache_cleanup.log"
printf '%s\n' "Verified results retained; cleaned only Nextflow run $name" > "$run/cache_cleanup.complete"
