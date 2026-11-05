#!/usr/bin/env python3
"""Report user-requested mean-TV and per-person percentage-point accuracy bands."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import model_family_gnomix as adapter


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)


def accuracy_summary(truth, prediction, targets):
    if (truth.shape != prediction.shape or truth.ndim != 2
            or truth.shape[0] == 0 or truth.shape[1] != len(SOURCES)):
        raise ValueError("Accuracy reporting requires matching nonempty seven-source arrays")
    tolerance = targets["comparison_tolerance"]
    for values in (truth, prediction):
        if (not np.isfinite(values).all() or np.any(values < -tolerance)
                or np.any(values > 1 + tolerance)
                or not np.allclose(values.sum(axis=1), 1, atol=tolerance, rtol=0)):
            raise ValueError("Accuracy reporting requires finite normalized ancestry fractions")
    absolute = np.abs(prediction - truth)
    tv, worst = absolute.sum(axis=1) / 2, absolute.max(axis=1)
    mean_tv = float(tv.mean())
    acceptable = worst <= targets["per_person_all_components"]["acceptable_absolute_error_max"] + tolerance
    desirable = worst <= targets["per_person_all_components"]["desirable_absolute_error_max"] + tolerance
    summary = dict(
        n_people=len(truth), mean_tv=mean_tv, p95_tv=float(np.quantile(tv, .95)),
        mean_tv_acceptable=mean_tv <= targets["mean_tv"]["acceptable_max"] + tolerance,
        mean_tv_desirable=mean_tv <= targets["mean_tv"]["desirable_max"] + tolerance,
        n_all_components_within_10pp=int(acceptable.sum()),
        fraction_all_components_within_10pp=float(acceptable.mean()),
        n_all_components_within_5pp=int(desirable.sum()),
        fraction_all_components_within_5pp=float(desirable.mean()),
        p95_worst_component_error=float(np.quantile(worst, .95)),
        coverage_minimum_specified=False)
    return summary, dict(tv=tv, worst_component_error=worst,
                         all_components_within_10pp=acceptable,
                         all_components_within_5pp=desirable)


def report(stack, arrays, output):
    if output.exists() and any(output.iterdir()):
        raise ValueError("Accuracy report requires a new or empty output")
    output.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def record(path, expected=None):
        actual = adapter.digest(path)
        if expected is not None and actual != expected:
            raise ValueError(f"Accuracy-report provenance changed: {path}")
        inputs[str(path)] = actual
        return path

    target_path = record(PROJECT / "config/phasefree_accuracy_targets.json")
    targets = json.loads(target_path.read_text())
    if (targets["sources"] != SOURCES
            or targets["per_person_all_components"]["minimum_coverage_fraction"] is not None):
        raise ValueError("This report expects seven-source bands without a cohort coverage gate")
    selected = json.loads(record(stack / "selected_candidate.json").read_text())
    verification = json.loads(record(
        stack / "verification/final_verification.json", selected["verification_sha256"]).read_text())
    completed = json.loads(record(
        stack / "completion.json", verification["completion_sha256"]).read_text())
    candidate = stack / "candidates" / selected["candidate"]
    record(candidate / "model.pkl", selected["model_sha256"])
    record(Path(__file__))
    summaries, people = [], []

    def add(frame, cohort, prediction_columns):
        if frame.sample_id.duplicated().any() or frame.sample_id.isna().any():
            raise ValueError(f"Duplicate or missing reporting individuals in {cohort}")
        truth = frame[[f"truth_{source}" for source in SOURCES]].to_numpy(float)
        prediction = frame[prediction_columns].to_numpy(float)
        summary, errors = accuracy_summary(truth, prediction, targets)
        summaries.append(dict(cohort=cohort, **summary))
        detail = pd.DataFrame(dict(cohort=cohort, sample_id=frame.sample_id.to_numpy(), **errors))
        if "family_id" in frame:
            detail["family_id"] = frame.family_id.to_numpy()
        people.append(detail)

    for track in ("pure", "pedigree"):
        path = candidate / f"{track}_predictions_scored.tsv"
        logical = str(Path("/stack") / path.relative_to(stack))
        record(path, completed["output_sha256"][logical])
        frame = pd.read_csv(path, sep="\t", float_precision="round_trip")
        add(frame, f"original_clean_{track}", [f"pred_{source}" for source in SOURCES])
    array_status = dict(status="pending_independent_verification",
                        note="No incomplete or unverified array predictions are scored.")
    if (arrays / "verification.json").exists():
        checked = json.loads(record(arrays / "verification.json").read_text())
        if (checked["status"] != "complete_verified_seen_donor_stress_test"
                or checked["model_sha256"] != selected["model_sha256"]):
            raise ValueError("Array accuracy report requires the matching verified frozen model")
        array_completion = json.loads(record(
            arrays / "completion.json", checked["completion_sha256"]).read_text())
        path = arrays / "predictions.tsv"
        record(path, array_completion["output_sha256"]["/validation/predictions.tsv"])
        frame = pd.read_csv(path, sep="\t", float_precision="round_trip")
        for (mask, track, stage), group in frame.groupby(["mask", "track", "stage"]):
            add(group, f"array_control_{mask}_{track}_{stage}", SOURCES)
        array_status = dict(
            status="included_independently_verified_seen_donor_stress_test",
            n_unique_people=int(frame.sample_id.nunique()),
            note="Masks and stages are separate paired conditions, not independent extra people.")
    summary_path, people_path = output / "summary.tsv", output / "per_person.tsv"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
    pd.concat(people, ignore_index=True).to_csv(people_path, sep="\t", index=False)
    for path, expected in inputs.items():
        if adapter.digest(path) != expected:
            raise ValueError(f"Accuracy-report input changed during reporting: {path}")
    result = dict(
        status="complete_supplemental_accuracy_report", candidate=selected["candidate"],
        model_sha256=selected["model_sha256"], targets=targets, cohorts=summaries,
        arrays=array_status, mean_tv_and_worst_component_error_are_distinct=True,
        reporting_targets_not_frozen_selection_gates=True,
        coverage_gate_evaluable=False, coverage_gate_reason="No minimum coverage fraction specified.",
        minor_detection_threshold_changed=False, original_selection_gates_changed=False,
        independent_final_test=False, production_promotion=False,
        input_sha256=inputs,
        output_sha256={str(path): adapter.digest(path) for path in (summary_path, people_path)})
    adapter.write_json(output / "accuracy_report.json", result)
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)
    print(f"Array status: {array_status['status']}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=Path("/stack"))
    parser.add_argument("--arrays", type=Path, default=Path("/validation"))
    parser.add_argument("--output", type=Path, default=Path("/accuracy"))
    args = parser.parse_args()
    report(args.stack, args.arrays, args.output)
