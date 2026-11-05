#!/usr/bin/env bash
set -euo pipefail

SEARCH_SPACE="$1"
OUTDIR="$2"
MAX_MODELS="${3:-0}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$OUTDIR/model_configs"
OUTDIR="$(cd "$OUTDIR" && pwd)"
TMP="$(mktemp -d "$OUTDIR/.model-configs.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

python3 "$ROOT/scripts/generate_model_configs.py" \
  --search-space "$SEARCH_SPACE" --outdir "$TMP" --max-models "$MAX_MODELS"

# The manifest must point at persistent configs, not the temporary generation.
awk -F '\t' -v OFS='\t' -v configs="$OUTDIR/model_configs/" '
  NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == "config_path") column = i
    if (!column) exit 1
  }
  NR > 1 {
    count = split($column, parts, "/")
    $column = configs parts[count]
  }
  { print }
' "$TMP/model_manifest.tsv" > "$TMP/model_manifest.fixed.tsv"
mv "$TMP/model_manifest.fixed.tsv" "$TMP/model_manifest.tsv"

for source in "$TMP"/*.tsv "$TMP"/*.json "$TMP"/model_configs/*.json; do
  target="$OUTDIR/${source#"$TMP"/}"
  if [[ -f "$target" ]]; then
    if cmp -s "$source" "$target"; then
      continue
    else
      status=$?
      [[ "$status" -eq 1 ]] || exit "$status"
    fi
  fi
  cp "$source" "$target"
done

for target in "$OUTDIR"/model_configs/*.json; do
  [[ -e "$TMP/model_configs/${target##*/}" ]] || rm "$target"
done
