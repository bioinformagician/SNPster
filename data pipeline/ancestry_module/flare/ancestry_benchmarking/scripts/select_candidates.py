#!/usr/bin/env python3
import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--root")
    inputs.add_argument("--manifest")
    parser.add_argument("--search-space", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--table", required=True)
    return parser.parse_args()


def metric_paths(root, manifest):
    if root:
        yield from Path(root).glob("candidate_*/metrics.json")
        return

    for line_number, raw_path in enumerate(
        Path(manifest).read_text().splitlines(), start=1
    ):
        if not raw_path.strip():
            continue
        metrics_path = Path(raw_path.strip()) / "metrics.json"
        if not metrics_path.is_file():
            raise FileNotFoundError(
                f"Manifest line {line_number} has no metrics.json: {metrics_path}"
            )
        yield metrics_path


a = parse_args()
s = json.loads(Path(a.search_space).read_text())
n = int(s["stage2"]["retain_top_candidates"])
fixed_ids = {
    x["candidate_id"] for x in s["stage2"].get("fixed_candidates", [])
}
rows = []
for p in metric_paths(a.root, a.manifest):
    m = json.loads(p.read_text())
    rows.append(
        {
            "candidate_id": m["candidate_id"],
            "pass": bool(m.get("pass", False)),
            "real_hard_pass": bool(m.get("real_hard_pass", False)),
            "score": float(m.get("selection_score", 999)),
            "dir": str(p.parent),
        }
    )
if not rows:
    raise SystemExit("No stage2 results found")
df = pd.DataFrame(rows).sort_values(
    ["pass", "real_hard_pass", "score"], ascending=[False, False, True]
)
keep = list(df.head(n).candidate_id)
for x in fixed_ids:
    if x in set(df.candidate_id) and x not in keep:
        keep.append(x)
outdir = Path(a.outdir)
outdir.mkdir(parents=True, exist_ok=True)
selected = df[df.candidate_id.isin(keep)].copy()
for _, r in selected.iterrows():
    shutil.copy2(
        Path(r.dir) / "candidate.json", outdir / f"{r.candidate_id}.json"
    )
selected.drop(columns=["dir"]).to_csv(a.table, sep="\t", index=False)
print(selected.drop(columns=["dir"]).to_string(index=False))
