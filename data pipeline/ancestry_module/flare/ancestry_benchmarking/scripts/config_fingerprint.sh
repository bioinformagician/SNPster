#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
SEARCH_SPACE="${1:-config/search_space.json}"

{
  sha256sum "$SEARCH_SPACE"
  find config -maxdepth 1 -type f \
    \( -name '*.json' -o -name '*.tsv' \) ! -name 'search_space*.json' -print0 \
    | sort -z \
    | xargs -0 -r sha256sum
  find scripts -maxdepth 1 -type f -name '*.py' -print0 \
    | sort -z \
    | xargs -0 -r sha256sum
} | sha256sum | awk '{print $1}'
