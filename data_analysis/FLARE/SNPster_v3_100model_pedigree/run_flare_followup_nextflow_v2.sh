#!/usr/bin/env bash
set -euo pipefail
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export FLARE_FOLLOWUP_OUTDIR="${FLARE_FOLLOWUP_OUTDIR:-results_v3_flare_nextflow_02}"
if [[ ! "$FLARE_FOLLOWUP_OUTDIR" =~ ^results_v3_flare_nextflow_[a-zA-Z0-9_]+$ ||
      $# -ne 1 || ("$1" != --prepare && "$1" != --smoke && "$1" != --run && "$1" != --resume) ]]; then
    echo "Usage: bash run_flare_followup_nextflow_v2.sh --prepare|--smoke|--run|--resume" >&2
    exit 1
fi
keep="${FLARE_KEEP_CACHE_FOR_RESUME_TEST:-0}"
FLARE_KEEP_CACHE_FOR_RESUME_TEST=1 bash "$project/run_flare_followup_nextflow.sh" "$1"
output="$project/$FLARE_FOLLOWUP_OUTDIR"
if [[ "$1" == --prepare ]]; then
    docker run --rm --network none --user "$(id -u):$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 \
        -v "$project:/workspace:ro" -v "$output:/followup:rw" \
        sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d python3 -c '
import json,sys
from pathlib import Path
sys.path.insert(0,"/workspace/pedigree_benchmark/scripts")
from pilot_matched_imputation import fingerprint_files
from model_family_gnomix import write_json
path=Path("/followup/source_sha256.json")
values=json.loads(path.read_text())
values.update(fingerprint_files([Path("/workspace/run_flare_followup_nextflow_v2.sh")]))
write_json(path,values)
'
    exit
fi
[[ "$keep" == 1 ]] && exit
mode=full
[[ "$1" == --smoke ]] && mode=smoke
run="$output/$mode"
cd "$run/launch"
if [[ ! -f "$run/published_verification.json" || "$PWD" != "$run/launch" ]]; then
    echo "Refusing cleanup: verified durable output and exact launch directory required." >&2
    exit 1
fi
session_id="$(awk -F '\t' 'END {print $6}' .nextflow/history)"
if [[ ! "$session_id" =~ ^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$ ]]; then
    echo "Cannot identify the exact completed Nextflow session." >&2
    exit 1
fi
NXF_OFFLINE=true nextflow clean -n "$session_id" > "$run/cache_cleanup.plan"
if grep '^Would remove ' "$run/cache_cleanup.plan" | grep -v -F "Would remove $run/work/" >/dev/null; then
    echo "Cleanup dry-run escaped this run's dedicated work directory." >&2
    exit 1
fi
NXF_OFFLINE=true nextflow clean -f "$session_id" > "$run/cache_cleanup.log"
if find "$run/work" -mindepth 3 -type f -print -quit | grep -q .; then
    echo "Cleanup left cached task files; refusing to mark cleanup complete." >&2
    exit 1
fi
printf '%s\n' "Verified published results retained; cleaned exact Nextflow session $session_id" \
    > "$run/cache_cleanup.complete"
