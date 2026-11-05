#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import pandas as pd


def product_dict(grid):
    if not grid:
        yield {}
        return
    keys = list(grid)
    for vals in itertools.product(*(grid[k] for k in keys)):
        yield dict(zip(keys, vals))


def stable_hash(obj, n=10):
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(raw.encode()).hexdigest()[:n]


def sanitize(x):
    return str(x).replace(".", "p").replace("-", "m").replace("+", "p")


ap = argparse.ArgumentParser()
ap.add_argument("--search-space", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--max-models", type=int, default=0)
args = ap.parse_args()

cfg = json.loads(Path(args.search_space).read_text())
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=True)
model_dir = out / "model_configs"
if model_dir.exists():
    for p in model_dir.glob("*.json"):
        p.unlink()
else:
    model_dir.mkdir(parents=True)

# ---------------------------------------------------------------------------
# FLARE panel-copying configurations
# ---------------------------------------------------------------------------
p = cfg["panel_search"]
panel_rows = []
for panel_cm, min_maf, min_mac in itertools.product(
    p["panel_cm"], p["min_maf"], p["min_mac"]
):
    raw = {"panel_cm": panel_cm, "min_maf": min_maf, "min_mac": min_mac}
    panel_rows.append(
        {
            "panel_id": f"panel_{stable_hash(raw, 8)}",
            **raw,
        }
    )

# ---------------------------------------------------------------------------
# Strategy-specific clustering grids
# ---------------------------------------------------------------------------
models = []
for panel in panel_rows:
    for strategy in cfg["flare_model_strategies"]:
        name = strategy["name"]
        family = strategy["family"]
        fixed = strategy.get("fixed", {})
        grid = strategy.get("grid", {})

        for params in product_dict(grid):
            row = {
                "strategy_name": name,
                "family": family,
                **fixed,
                **params,
                **panel,
            }
            if "k" not in row:
                raise SystemExit(f"Strategy {name} does not define k")
            row.setdefault("model_seed", cfg.get("seed", 20261005))
            row.setdefault("flare_seed", 20261002)
            row.setdefault("pca_components", 0)
            row.setdefault("compact_inactive", False)
            row.setdefault("refinement_mode", "update_p")
            row.setdefault("refine_em", True)

            h = stable_hash(row, 10)
            rep_id = f"{name}_K{row['k']}_{h}"
            row["rep_id"] = rep_id
            models.append(row)

# Deterministic ordering.
models.sort(
    key=lambda r: (
        r["strategy_name"],
        int(r["k"]),
        int(r.get("pca_components", 0) or 0),
        int(r.get("model_seed", 0)),
        r["rep_id"],
    )
)

# Smoke caps should span strategy families first, then K values within each
# strategy, instead of consuming one lexicographic bucket.
if args.max_models and args.max_models > 0:
    by_strategy = {}
    for row in models:
        by_strategy.setdefault(row["strategy_name"], []).append(row)

    for name in by_strategy:
        by_strategy[name] = sorted(
            by_strategy[name],
            key=lambda r: (int(r["k"]), int(r.get("pca_components", 0) or 0), r["rep_id"]),
        )

    strategy_order = sorted(by_strategy)
    selected = []
    depth = 0
    while len(selected) < args.max_models:
        added = False
        for name in strategy_order:
            if depth < len(by_strategy[name]):
                selected.append(by_strategy[name][depth])
                added = True
                if len(selected) >= args.max_models:
                    break
        if not added:
            break
        depth += 1
    models = selected

used_panels = {r["panel_id"] for r in models}
panel_rows = [r for r in panel_rows if r["panel_id"] in used_panels]

pd.DataFrame(panel_rows).to_csv(out / "panel_configs.tsv", sep="\t", index=False)

manifest = []
for row in models:
    path = model_dir / f"{row['rep_id']}.json"
    path.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
    manifest.append(
        {
            "panel_id": row["panel_id"],
            "rep_id": row["rep_id"],
            "strategy_name": row["strategy_name"],
            "family": row["family"],
            "k": row["k"],
            "pca_components": row.get("pca_components", 0),
            "model_seed": row.get("model_seed", ""),
            "flare_seed": row.get("flare_seed", ""),
            "config_path": str(path.resolve()),
        }
    )

pd.DataFrame(manifest).to_csv(out / "model_manifest.tsv", sep="\t", index=False)
(out / "model_configs.json").write_text(json.dumps(models, indent=2, sort_keys=True) + "\n")

print(
    f"Generated {len(models)} trainable FLARE models across "
    f"{len(set(r['strategy_name'] for r in models))} strategies and "
    f"{len(panel_rows)} panel-probability configurations"
)
for name, grp in pd.DataFrame(models).groupby("strategy_name"):
    print(f"  {name:24s} {len(grp):6d}")
