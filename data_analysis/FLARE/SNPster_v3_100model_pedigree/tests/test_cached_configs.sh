#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENERATOR="${1:-$ROOT/scripts/generate_model_configs_cached.sh}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

generate() {
  bash "$GENERATOR" "$ROOT/config/search_space_100.json" "$TMP/generated" "$1"
}

generate 0
test "$(find "$TMP/generated/model_configs" -name '*.json' | wc -l)" -eq 100
find "$TMP/generated" -type f -exec stat -c '%n %s %y' {} + | sort > "$TMP/before"
generate 0
find "$TMP/generated" -type f -exec stat -c '%n %s %y' {} + | sort > "$TMP/after"
diff -u "$TMP/before" "$TMP/after"
awk -F '\t' -v dir="$TMP/generated/model_configs/" '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "config_path") column = i }
  NR > 1 { if (index($column, dir) != 1) exit 1 }
  END { if (!column) exit 1 }
' "$TMP/generated/model_manifest.tsv"

CONFIG="$(find "$TMP/generated/model_configs" -name '*.json' -print -quit)"
cp "$CONFIG" "$TMP/expected.json"
printf 'changed\n' > "$CONFIG"
generate 0
cmp "$CONFIG" "$TMP/expected.json"

generate 8
test "$(find "$TMP/generated/model_configs" -name '*.json' | wc -l)" -eq 8
find "$TMP/generated" -type f -exec stat -c '%n %s %y' {} + | sort > "$TMP/before"
generate 8
find "$TMP/generated" -type f -exec stat -c '%n %s %y' {} + | sort > "$TMP/after"
diff -u "$TMP/before" "$TMP/after"
test "$(find "$TMP/generated" -maxdepth 1 -name '.model-configs.*' | wc -l)" -eq 0
echo 'PASS: unchanged configs retain timestamps; changed configs update, stale configs are removed, and manifest paths remain persistent.'
