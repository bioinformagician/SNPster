#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${IMAGE:-snpster-ancestry-benchmark:v1}"
MAX_MODELS="${MAX_MODELS:-8}"
MAX_STAGE2_CANDIDATES="${MAX_STAGE2_CANDIDATES:-80}"
SEARCH_SPACE="${SEARCH_SPACE:-config/search_space.json}"
HOST_UID="${SUDO_UID:-$(id -u)}"
HOST_GID="${SUDO_GID:-$(id -g)}"

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
  -v "${ROOT}:/workspace" \
  -v /srv/dependencies/ancestry/flare2_test:/data/flare2_test:ro \
  -v /srv/imputed:/data/imputed:ro \
  "${IMAGE}" \
  bash -lc './run_benchmark_v3.sh'
