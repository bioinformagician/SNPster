#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${IMAGE:-snpster-ancestry-benchmark:pedigree-100}"
HOST_UID="${SUDO_UID:-$(id -u)}"
HOST_GID="${SUDO_GID:-$(id -g)}"
PEDIGREE_DESIGN="${PEDIGREE_DESIGN:-reduced}"
RUN_OUTDIR="${RUN_OUTDIR:-/workspace/results_v3_sweep_100_disjoint_${PEDIGREE_DESIGN}}"
RUN_WORKDIR="${RUN_WORKDIR:-/workspace/work}"
[[ "$RUN_OUTDIR" == /workspace/* ]] || { echo 'RUN_OUTDIR must be under /workspace' >&2; exit 2; }
[[ "$RUN_WORKDIR" == /workspace/* ]] || { echo 'RUN_WORKDIR must be under /workspace' >&2; exit 2; }
mkdir -p "${ROOT}/${RUN_OUTDIR#/workspace/}" "${ROOT}/${RUN_WORKDIR#/workspace/}" "${ROOT}/.nextflow" "${ROOT}/config/generated/model_configs"
PDB_DOCKER_ARGS=()
if [[ "${PEDIGREE_ENABLED:-true}" == true ]]; then
 [[ -d /srv/dependencies/ancestry/hgdp_1kg_references ]] || { echo 'Missing HGDP + 1kGP reference directory' >&2; exit 2; }
 [[ -f /srv/dependencies/ancestry/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map ]] || { echo 'Missing GRCh38 genetic map' >&2; exit 2; }
 PDB_DOCKER_ARGS=(-v /srv/dependencies/ancestry/hgdp_1kg_references:/data/hgdp_1kg_references:ro)
 if [[ -n "${PEDSIM_HOST_DIR:-}" ]]; then
  [[ -x "$PEDSIM_HOST_DIR/ped-sim" ]] || { echo "Missing Ped-sim: $PEDSIM_HOST_DIR/ped-sim" >&2; exit 2; }
  PDB_DOCKER_ARGS+=(-v "$PEDSIM_HOST_DIR:/opt/ped-sim:ro")
 fi
fi
TTY=(); [[ -t 0 && -t 1 ]] && TTY=(-it)
docker run --rm "${TTY[@]}" --user "${HOST_UID}:${HOST_GID}" --cpus=20 --memory=180g \
 -e NXF_HOME=/workspace/.nextflow -e MAX_MODELS="${MAX_MODELS:-0}" \
 -e MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-6000}" \
 -e SEARCH_SPACE="${SEARCH_SPACE:-config/search_space_100.json}" \
 -e RUN_OUTDIR="$RUN_OUTDIR" -e RESUME_SESSION="${RESUME_SESSION:-}" \
 -e RUN_WORKDIR="$RUN_WORKDIR" -e FROZEN_CANDIDATES_MANIFEST="${FROZEN_CANDIDATES_MANIFEST:-}" \
 -e PEDIGREE_REUSE_DIR="${PEDIGREE_REUSE_DIR:-}" \
 -e RECOVERY_STAGE2_MANIFEST="${RECOVERY_STAGE2_MANIFEST:-}" \
 -e RECOVERY_SOURCE_OUTDIR="${RECOVERY_SOURCE_OUTDIR:-$RUN_OUTDIR}" \
 -e PEDIGREE_DESIGN="$PEDIGREE_DESIGN" \
 -e PEDIGREE_ENABLED="${PEDIGREE_ENABLED:-true}" -e PEDIGREE_FAMILIES="${PEDIGREE_FAMILIES:-90}" \
 -e PEDIGREE_BATCH_SIZE="${PEDIGREE_BATCH_SIZE:-10}" -e PEDIGREE_SEED="${PEDIGREE_SEED:-20261008}" \
 -e PEDIGREE_ALLOW_REUSE="${PEDIGREE_ALLOW_REUSE:-true}" \
 -e PEDIGREE_KEEP_HIGH_DENSITY="${PEDIGREE_KEEP_HIGH_DENSITY:-false}" \
 -e PEDIGREE_CHROMOSOMES="${PEDIGREE_CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}" \
 -v "$ROOT:/workspace" -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
 -v /srv/imputed:/data/imputed:ro "${PDB_DOCKER_ARGS[@]}" \
 "$IMAGE" bash -lc 'bash ./run_benchmark_v3.sh'
