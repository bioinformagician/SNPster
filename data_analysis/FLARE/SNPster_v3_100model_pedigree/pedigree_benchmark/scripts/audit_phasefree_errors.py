#!/usr/bin/env python3
"""Diagnose frozen development errors without fitting models or changing gates."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from screen_v1_decoders import predict_decoder


SOURCES = list(adapter.SOURCES)
MODEL_SHA = "95f3e3c9abaf76fb357d4ef51b2e29bf212a033002c02a04b825ec01b0231686"
FEATURE_COLUMNS = [f"{view}_{source}" for view in ("classifier", "frequency") for source in SOURCES]


def read_table(path):
    # Histogram-tree split boundaries require the exact serialized feature values.
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def validate_condition_axis(frame):
    if (frame[["group", "sample_id"]].isna().any().any()
            or frame.duplicated(["group", "sample_id"]).any()):
        raise ValueError("Duplicate or incomplete frozen condition/sample axis")


def ancestry_cells(frame, cfg):
    required = ["sample_id", "family_id",
                *[f"{prefix}_{source}" for prefix in ("truth", "pred", "classifier", "frequency")
                  for source in SOURCES]]
    if frame.sample_id.duplicated().any() or frame[required].isna().any().any():
        raise ValueError("Duplicate or incomplete diagnostic sample/feature axis")
    records = []
    for row in frame.to_dict("records"):
        for prefix in ("truth", "pred", "classifier", "frequency"):
            if not np.isclose(sum(float(row[f"{prefix}_{source}"]) for source in SOURCES),
                              1, atol=1e-8, rtol=0):
                raise ValueError("Invalid diagnostic ancestry fraction sum")
        for source in SOURCES:
            truth = float(row[f"truth_{source}"])
            values = {view: float(row[f"{prefix}_{source}"]) for view, prefix in (
                ("decoder", "pred"), ("classifier", "classifier"), ("frequency", "frequency"))}
            if not 0 <= truth <= 1 or any(
                    not np.isfinite(value) or not -1e-12 <= value <= 1 + 1e-12
                    for value in values.values()):
                raise ValueError("Invalid diagnostic ancestry fraction")
            records.append(dict(
                sample_id=row["sample_id"], family_id=row["family_id"], source=source,
                truth=truth, minor=cfg["minor_min"] <= truth <= cfg["minor_max"],
                absent=truth < cfg["absent_ancestry_max"], **values))
    return pd.DataFrame(records)


def cell_summary(cells, cfg):
    minor, absent = cells[cells.minor], cells[cells.absent]
    views = {}
    for view in ("classifier", "frequency", "decoder"):
        detected = int((minor[view] >= cfg["minor_detection_threshold"]).sum())
        false_positive = int((absent[view] >= cfg["minor_detection_threshold"]).sum())
        views[view] = dict(
            minor_detected=detected,
            minor_recall=detected / len(minor) if len(minor) else None,
            minor_mean_absolute_error=float(np.abs(minor[view] - minor.truth).mean())
            if len(minor) else None,
            absent_flagged=false_positive,
            absent_flagged_rate=false_positive / len(absent) if len(absent) else None,
        )
    return dict(
        minor_contributions=len(minor), minor_people=minor.sample_id.nunique(),
        minor_families=minor.family_id.nunique(), absent_contributions=len(absent),
        minor_gate_evaluable=len(minor) >= cfg["minor_gate_min_support"],
        views=views)


def audit(stack, output):
    if output.exists() and any(output.iterdir()):
        raise ValueError("Diagnostic output must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def record(path):
        inputs[str(path)] = adapter.digest(path)
        return path

    selected = json.loads(record(stack / "selected_candidate.json").read_text())
    verified_path = record(stack / "verification/final_verification.json")
    if inputs[str(verified_path)] != selected["verification_sha256"]:
        raise ValueError("Independent candidate verification changed")
    verified = json.loads(verified_path.read_text())
    completion_path = record(stack / "completion.json")
    if inputs[str(completion_path)] != verified["completion_sha256"]:
        raise ValueError("Verified stacking completion changed")
    completed = json.loads(completion_path.read_text())
    directory = stack / "candidates" / selected["candidate"]
    model_path = record(directory / "model.pkl")
    if inputs[str(model_path)] != MODEL_SHA or selected["model_sha256"] != MODEL_SHA:
        raise ValueError("Frozen candidate digest changed")
    paths = [stack / "features.tsv", stack / "training.tsv", directory / "metrics.json",
             directory / "pure_predictions_scored.tsv",
             directory / "pedigree_predictions_scored.tsv"]
    for path in paths:
        logical = str(Path("/stack") / path.relative_to(stack))
        record(path)
        if inputs[str(path)] != completed["output_sha256"][logical]:
            raise ValueError(f"Verified diagnostic input changed: {path}")
    cfg_path = Path(__file__).resolve().parents[2] / "config/unified_selection.json"
    cfg = json.loads(record(cfg_path).read_text())
    experiment = json.loads(record(stack / "experiment.json").read_text())
    if inputs[str(cfg_path)] != experiment["implementation_sha256"]["/workspace/config/unified_selection.json"]:
        raise ValueError("Unchanged gate configuration differs from frozen experiment")
    record(Path(__file__))
    record(Path(adapter.__file__))
    record(Path(__file__).with_name("screen_v1_decoders.py"))
    feature_frame = read_table(stack / "features.tsv")
    training = read_table(stack / "training.tsv")
    validate_condition_axis(feature_frame)
    validate_condition_axis(training)
    with model_path.open("rb") as stream:
        model = pickle.load(stream)
    metrics = json.loads((directory / "metrics.json").read_text())
    cells, group_metrics = [], []
    max_discrepancy = 0.
    with threadpool_limits(limits=2):
        for track, expected_n in (("pure", 159), ("pedigree", 123)):
            scored = read_table(directory / f"{track}_predictions_scored.tsv")
            values = feature_frame[feature_frame.group.eq(track)]
            merged = scored.merge(values, on="sample_id", validate="one_to_one",
                                  how="outer", indicator=True)
            if len(merged) != expected_n or not merged._merge.eq("both").all():
                raise ValueError(f"Diagnostic features do not cover exact {track} evaluation")
            reproduced = predict_decoder(model["decoder"], merged[FEATURE_COLUMNS].to_numpy(float))
            expected = merged[[f"pred_{source}" for source in SOURCES]].to_numpy(float)
            discrepancy = float(np.max(np.abs(expected - reproduced)))
            if discrepancy > 1e-12:
                raise ValueError("Frozen decoder predictions do not reproduce")
            max_discrepancy = max(max_discrepancy, discrepancy)
            truth = merged[[f"truth_{source}" for source in SOURCES]].to_numpy(float)
            tv = np.abs(reproduced - truth).sum(axis=1) / 2
            if abs(float(tv.mean()) - metrics[track]["mean_tv"]) > 1e-12:
                raise ValueError("Diagnostic TV differs from independently verified metrics")
            cells.append(ancestry_cells(merged, cfg).assign(track=track))
            for column in (("source_group", "fine") if track == "pure" else ("scenario",)):
                for value, group in merged.groupby(column):
                    group_metrics.append(dict(
                        track=track, grouping=column, value=value, n=len(group),
                        mean_tv=float(group.TV_distance.mean()),
                        p95_tv=float(group.TV_distance.quantile(.95)),
                        mean_off_target=float(group.off_target_mass.mean())))
    cells = pd.concat(cells, ignore_index=True)
    summary = cell_summary(cells, cfg)
    if (summary["minor_contributions"] != metrics["unified"]["minor_2_5_n"]
            or abs(summary["views"]["decoder"]["minor_recall"]
                   - metrics["unified"]["minor_2_5_recall"]) > 1e-12):
        raise ValueError("Minor-ancestry diagnostic differs from unchanged gate")
    clean = training[training.weight_clean > 0]
    if len(clean) != 799 or not np.isclose(clean.weight_clean.sum(), 1):
        raise ValueError("Unexpected frozen clean training policy")
    supports = []
    for source in SOURCES:
        minor = clean[f"truth_{source}"].between(cfg["minor_min"], cfg["minor_max"])
        supports.append(dict(
            source=source, minor_training_people=int(minor.sum()),
            minor_clean_training_weight=float(clean.loc[minor, "weight_clean"].sum()),
            evaluation=cell_summary(cells[cells.source.eq(source)], cfg)))
    cells.to_csv(output / "ancestry_cells.tsv", sep="\t", index=False)
    cells[cells.minor].to_csv(output / "minor_contributions.tsv", sep="\t", index=False)
    pd.DataFrame(group_metrics).to_csv(output / "group_errors.tsv", sep="\t", index=False)
    for path, expected in inputs.items():
        if adapter.digest(path) != expected:
            raise ValueError(f"Diagnostic input changed during audit: {path}")
    report = dict(
        status="complete_frozen_development_error_audit", candidate=selected["candidate"],
        model_sha256=MODEL_SHA, model_fitted=False, thresholds_changed=False,
        new_validation_or_accuracy_claim=False, independent_final_test=False,
        classifier_frequency_views_are_uncalibrated_not_alternative_validated_models=True,
        detection_threshold=cfg["minor_detection_threshold"],
        minor_range=[cfg["minor_min"], cfg["minor_max"]],
        n_clean_training=799, n_development_people=282,
        max_saved_feature_decoder_reproduction_error=max_discrepancy,
        minor_and_absent_summary=summary, by_source=supports,
        input_sha256=inputs,
        output_sha256={str(path): adapter.digest(path) for path in output.iterdir() if path.is_file()})
    adapter.write_json(output / "error_audit.json", report)
    print(json.dumps({key: report[key] for key in (
        "minor_and_absent_summary", "by_source")}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=Path("/stack"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.stack, args.output)
