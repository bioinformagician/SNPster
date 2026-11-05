#!/usr/bin/env python3
"""Validate durable published results before permitting exact-run cache cleanup."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from verify_flare_reference_density import digest


SOURCES = ("EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA")
AXES = ["panel", "chromosome", "candidate", "stage", "sample_id"]


def read(path):
    return json.loads(Path(path).read_text())


def verify(run):
    prepared, summary = run / "preflight/prepared", run / "summary/summary"
    characterization = run / "characterization/characterization/verification.json"
    plan = read(prepared / "jobs.json")
    experiment = read(prepared / "experiment.json")
    audit = read(summary / "verification.json")
    if audit["new_tasks"] != len(plan) or audit["mode"] != experiment["mode"]:
        raise ValueError("Published result task/mode counts differ")
    selection = read(prepared / "selection.json")
    old = pd.read_csv(prepared / "reused_predictions.tsv", sep="\t", float_precision="round_trip")
    metadata = {row["sample_id"]: row for row in selection["people"]}
    hashes, records = {}, []
    for job in plan:
        directory = run / "calls" / job["id"]
        record = read(directory / "prediction.json")
        row = record["row"]
        if tuple(row[k] for k in AXES) != tuple(job[k] for k in AXES):
            raise ValueError("Published task axes differ")
        for original, expected in record["output_sha256"].items():
            path = directory / Path(original).name
            if digest(path) != expected:
                raise ValueError(f"Published local ancestry/model differs: {path}")
            hashes[str(path)] = expected
        hashes[str(directory / "prediction.json")] = digest(directory / "prediction.json")
        p = np.array([row[f"pred_{s}"] for s in SOURCES])
        truth = np.array(selection["truth"][str(row["chromosome"])][row["sample_id"]])
        np.testing.assert_allclose(truth, [row[f"truth_{s}"] for s in SOURCES], rtol=0, atol=1e-12)
        if not np.isfinite(p).all() or (p < 0).any() or abs(p.sum() - 1) > 1e-10:
            raise ValueError("Invalid durable posterior fractions")
        error = np.abs(p - truth)
        for key, value in (("tv", error.sum() / 2), ("max_component_error", error.max()),
                           ("off_target_mass", p[truth == 0].sum())):
            if abs(row[key] - value) > 1e-12:
                raise ValueError(f"Durable metric differs: {key}")
        for key in ("track", "family_id", "source_group", "fine", "scenario"):
            if row[key] != metadata[row["sample_id"]][key]:
                raise ValueError("Durable metadata differs from training selection")
        records.append(row)
    new = pd.DataFrame(records)
    expected_recovery = pd.concat([old, new[new.kind.eq("recovery")]], ignore_index=True)
    native = old[old.chromosome.eq(12) & old.panel.isin(["lonni", "frederik"])
                 & old.candidate.isin(["full_em", "finecap16_em", "size_control_em"])]
    if experiment["mode"] == "smoke":
        native = native[native.sample_id.isin(new.sample_id)]
    expected_theta = pd.concat([native, new[new.kind.eq("theta")]], ignore_index=True)
    for prefix, expected in (("recovered_screen", expected_recovery), ("crossed_theta", expected_theta)):
        predictions = pd.read_csv(summary / f"{prefix}_predictions.tsv", sep="\t", float_precision="round_trip")
        if predictions.duplicated(AXES).any() or len(predictions) != len(expected):
            raise ValueError("Durable summary duplicate/missing predictions")
        a, b = predictions.sort_values(AXES), expected.sort_values(AXES)
        if list(a[AXES].itertuples(index=False, name=None)) != list(b[AXES].itertuples(index=False, name=None)):
            raise ValueError("Durable prediction keys differ")
        for key in [f"pred_{s}" for s in SOURCES] + ["tv", "max_component_error", "off_target_mass"]:
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-12)
        summaries = pd.read_csv(summary / f"{prefix}_comparison.tsv", sep="\t", float_precision="round_trip")
        group_axes = ["panel", "chromosome", "candidate", "stage", "track"]
        groups = dict(tuple(predictions.groupby(group_axes)))
        if len(groups) != len(summaries):
            raise ValueError("Durable summary group coverage differs")
        seen = set()
        for row in summaries.to_dict("records"):
            group_key = tuple(row[key] for key in group_axes)
            if group_key in seen:
                raise ValueError("Duplicate durable summary group")
            seen.add(group_key)
            group = groups[group_key]
            values = dict(mean_tv=group.tv.mean(), p95_tv=group.tv.quantile(.95), max_tv=group.tv.max(),
                          coverage10=(group.max_component_error <= .1 + 1e-12).mean(),
                          coverage5=(group.max_component_error <= .05 + 1e-12).mean(),
                          family_balanced_tv=group.groupby("family_id").tv.mean().mean(),
                          off_target_mass=group.off_target_mass.mean())
            if row["n_people"] != len(group) or any(abs(row[key] - value) > 1e-12 for key, value in values.items()):
                raise ValueError("Durable summary metrics differ")
    if experiment["mode"] == "full":
        wanted = set(itertools.product(
            ("lonni", "frederik", "frederik_thin"), (22, 12, 1),
            ("full_fixed", "full_em", "finecap16_fixed", "finecap16_em", "size_control_fixed", "size_control_em"),
            ("inherited_phase", "array_phase_only"), metadata))
        actual = set(expected_recovery[AXES].itertuples(index=False, name=None))
        if actual != wanted or len(expected_theta) != 7584:
            raise ValueError("Full follow-up Cartesian coverage incomplete")
    for path in itertools.chain(prepared.rglob("*"), summary.rglob("*"),
                                (characterization,), (run / "trace.tsv",)):
        if path.is_file():
            hashes[str(path)] = digest(path)
    for original, expected in audit["output_sha256"].items():
        if digest(summary / Path(original).name) != expected:
            raise ValueError("Published summary hash differs")
    result = dict(status="verified_durable_outputs_safe_for_exact_run_cleanup",
                  mode=experiment["mode"], new_calls=len(plan),
                  recovered_screen_predictions=len(expected_recovery),
                  theta_with_native_baseline_predictions=len(expected_theta),
                  original_failure_preserved=True, whole_genome_score=False,
                  source_sha256=experiment["source_sha256"], verifier_sha256=digest(__file__),
                  published_output_sha256=hashes)
    (run / "published_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if not k.endswith("sha256")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    verify(args.run_dir)
