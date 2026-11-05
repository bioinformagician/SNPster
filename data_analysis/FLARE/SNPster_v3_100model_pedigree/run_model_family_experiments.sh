#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MODEL_FAMILY_SCREEN_OUTDIR="${MODEL_FAMILY_SCREEN_OUTDIR:-results_v3_model_family_screen_05}"
output_name="$MODEL_FAMILY_SCREEN_OUTDIR"
if [[ ! "$output_name" =~ ^results_v3_model_family_screen_[a-zA-Z0-9_]+$ ]]; then
    echo "Invalid isolated model-family output directory" >&2
    exit 1
fi
if [[ $# -gt 1 || ($# != 0 && "${1:-}" != --smoke-only &&
      "${1:-}" != --fast-smoke-verified && "${1:-}" != --resume-verified-smokes &&
      "${1:-}" != --resume && "${1:-}" != --verify-resume && "${1:-}" != --pause) ]]; then
    echo "Usage: bash run_model_family_experiments.sh [--smoke-only|--fast-smoke-verified|--resume-verified-smokes|--resume|--verify-resume|--pause]" >&2
    exit 1
fi
results="$project_dir/$output_name"
if [[ ! -f "$results/input_audit.json" || ! -f "$results/initialization/initialization_audit.json" ||
      ! -f "$results/baseline/metrics.json" ]]; then
    echo "Prepare inputs, initialize FLARE models, and run the matched baseline first" >&2
    exit 1
fi
mode="${1:-}"
if [[ "$mode" == --pause ]]; then
    python3 - "$results" <<'PY'
import datetime
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
request = root / "pause_requested.json"
staging = request.with_suffix(".pending")
staging.write_text(json.dumps({
    "status": "pause_requested",
    "requested_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}, indent=2) + "\n")
staging.replace(request)
PY
    echo "Pause requested. Gnomix fits stop after saving/evaluating their current chromosome."
    echo "The scheduler will not start another candidate; completed work remains on disk."
    exit 0
fi
gnomix_image="snpster-ancestry-benchmark:gnomix-screen"
image_id="$(docker image inspect "$gnomix_image" --format '{{.Id}}')"
if [[ "$mode" == --resume || "$mode" == --verify-resume ]]; then
    if [[ ! -f "$results/gnomix_image.txt" || "$image_id" != "$(cat "$results/gnomix_image.txt")" ]]; then
        echo "Gnomix image changed or original image provenance is missing; refusing resume" >&2
        exit 1
    fi
else
    printf '%s\n' "$image_id" > "$results/gnomix_image.txt"
fi
mkdir "$results/experiment_scheduler.lock"
active_names=()
cleanup() {
    local status="$?"
    trap - EXIT INT TERM
    if [[ "$status" != 0 ]]; then
        for container in "${active_names[@]}"; do
            if [[ "$(docker inspect --format '{{.State.Running}}' "$container" 2>/dev/null || true)" == true ]]; then
                docker stop "$container" >/dev/null
            fi
        done
        if [[ "$status" == 75 ]]; then
            echo "Screen paused at saved checkpoints; use --resume to continue." >&2
        else
            echo "Screen stopped with status $status; incomplete results are not a leaderboard." >&2
        fi
    fi
    rmdir "$results/experiment_scheduler.lock"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

verify_flare_result() {
    python3 - "$results" "$1" "$2" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root, candidate, phase = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
directory = root / ("smoke" if phase == "smoke" else "candidates") / candidate
audit = json.loads((directory / "persistence_audit.json").read_text())
if audit.get("model_roundtrip_verified") is not True or audit["smoke_only"] != (phase == "smoke"):
    raise ValueError("FLARE result is not a verified fit for the requested phase")
for name, key in (("predictions.tsv", "prediction_sha256"),
                  ("predictions_physical.tsv", "physical_prediction_sha256"),
                  ("model.final.model", "model_sha256")):
    if hashlib.sha256((directory / name).read_bytes()).hexdigest() != audit[key]:
        raise ValueError(f"FLARE result changed: {candidate}/{name}")
if phase != "smoke":
    json.loads((directory / "metrics.json").read_text())
print("Verified saved FLARE", phase, candidate)
PY
}

run_candidate() {
    local candidate="$1" phase="$2"
    local args=(--outdir /results --candidate "$candidate")
    local name="snpster-${output_name}-gnomix-${candidate}-${phase}"
    local score_required=1
    if [[ "$candidate" == flare_* ]]; then
        if [[ "$mode" == --verify-resume ||
              ("$mode" == --resume && -f "$results/candidates/$candidate/metrics.json") ]]; then
            verify_flare_result "$candidate" "$phase"
            return
        fi
        if [[ "$phase" == smoke ]]; then
            bash "$project_dir/run_model_family_screen.sh" --flare "$candidate" smoke
        else
            bash "$project_dir/run_model_family_screen.sh" --flare "$candidate"
        fi
        return
    fi
    if [[ "$phase" == smoke ]]; then args+=(--smoke); fi
    if [[ "$mode" == --resume || "$mode" == --verify-resume ]]; then args+=(--resume); fi
    if [[ "$mode" == --verify-resume ]]; then args+=(--verify-only); fi
    if [[ "$mode" == --resume && "$phase" == full &&
          -f "$results/candidates/$candidate/metrics.json" &&
          -f "$results/candidates/$candidate/audit.json" ]]; then
        score_required="$(python3 - "$results/candidates/$candidate/audit.json" <<'PY'
import json
import sys
with open(sys.argv[1]) as handle:
    audit = json.load(handle)
print(0 if audit["status"] == "complete" and audit["stage"] != "fit" else 1)
PY
)"
    fi
    docker run --rm --network none --name "$name" --cpus=6 --memory=40g \
        --user "$(id -u):$(id -g)" \
        -e PYTHONUNBUFFERED=1 -e PYTHONDONTWRITEBYTECODE=1 \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        -e NUMBA_CACHE_DIR=/tmp/numba-cache -e MPLCONFIGDIR=/tmp/matplotlib \
        -e XDG_CACHE_HOME=/tmp/cache -e TMPDIR=/results \
        -v "$project_dir:/workspace:ro" \
        -v "$results:/results:rw" \
        -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
        -w /tmp "$gnomix_image" \
        python3 /workspace/pedigree_benchmark/scripts/model_family_gnomix.py "${args[@]}" \
        2>&1 | tee -a "$results/${candidate}_${phase}.log"
    if [[ "$phase" == full && "$mode" != --verify-resume && "$score_required" == 1 ]]; then
        bash "$project_dir/run_model_family_screen.sh" --score "$candidate"
    fi
}

run_pair() {
    local first="$1" second="$2" phase="$3" first_pid second_pid first_status=0 second_status=0
    active_names=(
        "snpster-${output_name}---flare-${first}" "snpster-${output_name}---flare-${second}"
        "snpster-${output_name}-gnomix-${first}-${phase}"
        "snpster-${output_name}-gnomix-${second}-${phase}"
    )
    echo "Starting $phase: $first and $second (maximum 12 CPUs, 80 GiB limits)"
    run_candidate "$first" "$phase" &
    first_pid="$!"
    run_candidate "$second" "$phase" &
    second_pid="$!"
    wait "$first_pid" || first_status="$?"
    wait "$second_pid" || second_status="$?"
    if [[ ("$first_status" == 75 || "$second_status" == 75) &&
          ("$first_status" == 0 || "$first_status" == 75) &&
          ("$second_status" == 0 || "$second_status" == 75) ]]; then
        return 75
    fi
    if [[ "$first_status" != 0 || "$second_status" != 0 ]]; then
        echo "Failed $phase: $first status=$first_status; $second status=$second_status" >&2
        return 1
    fi
    active_names=()
}

flare=(flare_direct_fixed flare_direct_refined flare_fine_fixed flare_fine_refined
       flare_subclusters_fixed flare_subclusters_refined)
gnomix=(gnomix_fast_020 gnomix_fast_050 gnomix_default_020 gnomix_default_050
        gnomix_best_020 gnomix_best_050)
all_candidates=("${flare[@]}" "${gnomix[@]}")
if [[ "$mode" == --resume || "$mode" == --verify-resume ]]; then
    bash "$project_dir/run_model_family_screen.sh" --verify-inputs
    python3 - "$results" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
state = root / "pause_state.json"
if state.exists():
    for relative, expected in json.loads(state.read_text()).get("saved_result_sha256", {}).items():
        path = root / relative
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("Invalid checkpoint snapshot path")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Saved result changed since pause: {relative}")
PY
    if [[ ! -f "$results/smoke_complete.txt" ]]; then
        echo "Resume requires all twelve completed compatibility checks" >&2
        exit 1
    fi
    for candidate in "${all_candidates[@]}"; do
        if [[ "$candidate" == flare_* ]]; then
            verify_flare_result "$candidate" smoke
        else
            mode=--verify-resume run_candidate "$candidate" smoke
        fi
    done
    if [[ "$mode" == --verify-resume ]]; then
        for candidate in "${all_candidates[@]}"; do
            if [[ "$candidate" == flare_* ]]; then
                if [[ -f "$results/candidates/$candidate/metrics.json" ]]; then
                    verify_flare_result "$candidate" full
                else
                    echo "Pending full candidate: $candidate"
                fi
            elif [[ -f "$results/candidates/$candidate/audit.json" ]]; then
                run_candidate "$candidate" full
            elif [[ -d "$results/candidates/$candidate/models" ]]; then
                echo "Model directory lacks fit provenance: $candidate" >&2
                exit 1
            else
                echo "Pending full candidate: $candidate"
            fi
        done
        echo "Resume verification completed without fitting or changing saved predictions."
        exit 0
    fi
    python3 - "$results" <<'PY'
from pathlib import Path
import sys
(Path(sys.argv[1]) / "pause_requested.json").unlink(missing_ok=True)
PY
fi
if [[ "$mode" != --resume ]]; then
for ((i=0; i<12; i+=2)); do
    if [[ "${1:-}" == --resume-verified-smokes && "$i" -lt 4 ]]; then
        python3 - "$results" "${all_candidates[i]}" "${all_candidates[i+1]}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
recovery = json.loads((root / "recovery_audit.json").read_text())
for candidate in sys.argv[2:]:
    directory = root / "smoke" / candidate
    audit = json.loads((directory / "persistence_audit.json").read_text())
    if candidate not in recovery["reused_smokes"] or not audit["model_roundtrip_verified"]:
        raise ValueError("Requested smoke was not part of verified recovery")
    for name, key in (("predictions.tsv", "prediction_sha256"),
                      ("predictions_physical.tsv", "physical_prediction_sha256"),
                      ("model.final.model", "model_sha256")):
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != audit[key]:
            raise ValueError("Recovered smoke checksum changed")
print("Reusing verified smokes:", ", ".join(sys.argv[2:]))
PY
        continue
    fi
    if [[ ("${1:-}" == --fast-smoke-verified || "${1:-}" == --resume-verified-smokes) &&
          "${all_candidates[i]}" == gnomix_fast_020 ]]; then
        if [[ ! -f "$results/smoke/gnomix_fast_020/output_audit.json" ||
              ! -f "$results/smoke/gnomix_fast_020/predictions.tsv" ]]; then
            echo "The independently verified fast020 smoke artifacts are missing" >&2
            exit 1
        fi
        python3 - "$results" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
audit = json.loads((root / "smoke/gnomix_fast_020/output_audit.json").read_text())
implementation = root.parent / "pedigree_benchmark/scripts/model_family_gnomix.py"
if (audit.get("status") != "complete" or audit.get("smoke_only") is not True
        or audit.get("model_roundtrip_verified") is not True
        or audit.get("candidate") != "gnomix_fast_020"
        or audit.get("adapter_sha256") != hashlib.sha256(implementation.read_bytes()).hexdigest()):
    raise ValueError("Existing fast020 smoke is incomplete or its implementation changed")
PY
        active_names=("snpster-${output_name}-gnomix-gnomix_fast_050-smoke")
        run_candidate gnomix_fast_050 smoke
        active_names=()
        continue
    fi
    run_pair "${all_candidates[i]}" "${all_candidates[i+1]}" smoke
done
fi
printf 'All twelve chromosome-22 compatibility tests completed; not an accuracy ranking.\n' \
    > "$results/smoke_complete.txt"
if [[ "${1:-}" == --smoke-only ]]; then exit 0; fi
for ((i=0; i<12; i+=2)); do
    if [[ -f "$results/pause_requested.json" ]]; then exit 75; fi
    run_pair "${all_candidates[i]}" "${all_candidates[i+1]}" full
done
if [[ -f "$results/pause_requested.json" ]]; then exit 75; fi
bash "$project_dir/run_model_family_screen.sh" --summarize
