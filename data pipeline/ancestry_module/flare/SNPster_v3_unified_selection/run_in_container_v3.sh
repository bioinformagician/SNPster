#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${IMAGE:-snpster-ancestry-benchmark:v1}"
MAX_MODELS="${MAX_MODELS:-8}"
MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-80}"
SEARCH_SPACE="${SEARCH_SPACE:-config/search_space.json}"
HOST_UID="${SUDO_UID:-$(id -u)}"
HOST_GID="${SUDO_GID:-$(id -g)}"

# Ped-sim reference mount is opt-in. The original command is unchanged when
# PEDIGREE_ENABLED=false (the default).
PDB_DOCKER_ARGS=()
if [[ "${PEDIGREE_ENABLED:-true}" == "true" ]]; then
  PEDSIM_HOST_DIR="${PEDSIM_HOST_DIR:-/srv/dependencies/ancestry/ped-sim}"
  [[ -x "${PEDSIM_HOST_DIR}/ped-sim" ]] || { echo "Missing Ped-sim: ${PEDSIM_HOST_DIR}/ped-sim" >&2; exit 2; }
  [[ -d /srv/dependencies/ancestry/hgdp_1kg_references ]] || { echo 'Missing HGDP+1kGP' >&2; exit 2; }
  [[ -f /srv/dependencies/imputation_runner/imputer/plink.GRCh38.map ]] || { echo 'Missing GRCh38 PLINK map' >&2; exit 2; }
  PDB_DOCKER_ARGS=(
    -v /srv/dependencies/ancestry/hgdp_1kg_references:/data/hgdp_1kg_references:ro
    -v /srv/dependencies/imputation_runner/imputer:/data/imputer:ro
    -v "${PEDSIM_HOST_DIR}:/opt/ped-sim:ro"
  )
fi

mkdir -p "${ROOT}/results_v3" "${ROOT}/work" "${ROOT}/.nextflow" "${ROOT}/config/generated/model_configs"

echo "======================================================================"
echo "SNPster ancestry strategy laboratory v3"
echo "======================================================================"
echo "Image:                   ${IMAGE}"
echo "Benchmark dir:           ${ROOT}"
echo "Search space:            ${SEARCH_SPACE}"
echo "Trainable FLARE models:  ${MAX_MODELS} (0 = all)"
echo "Stage2 candidates cap:   ${MAX_STAGE2_CANDIDATES} (0 = all)"
echo

docker run --rm -it \
  --user "${HOST_UID}:${HOST_GID}" \
  --cpus=20 \
  --memory=180g \
  -e NXF_HOME=/workspace/.nextflow \
  -e MAX_MODELS="${MAX_MODELS}" \
  -e MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES}" \
  -e SEARCH_SPACE="${SEARCH_SPACE}" \
  -e PEDIGREE_ENABLED="${PEDIGREE_ENABLED:-true}" \
  -e PEDIGREE_FAMILIES="${PEDIGREE_FAMILIES:-90}" \
  -e PEDIGREE_BATCH_SIZE="${PEDIGREE_BATCH_SIZE:-10}" \
  -e PEDIGREE_SEED="${PEDIGREE_SEED:-20261008}" \
  -e PEDIGREE_ALLOW_REUSE="${PEDIGREE_ALLOW_REUSE:-true}" \
  -e PEDIGREE_KEEP_HIGH_DENSITY="${PEDIGREE_KEEP_HIGH_DENSITY:-false}" \
  -e PEDIGREE_CHROMOSOMES="${PEDIGREE_CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}" \
  -v "${ROOT}:/workspace" \
  -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
  -v /srv/imputed:/data/imputed:ro \
  "${PDB_DOCKER_ARGS[@]}" \
  "${IMAGE}" \
  bash -lc './run_benchmark_v3.sh'
