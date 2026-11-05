#!/usr/bin/env python3
"""Matched pure-only versus pure-plus-pedigree decoder development screen."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from screen_v1_decoders import (
    SOURCES, REPRESENTATION, SEED, digest, exact_join, family_donors,
    feature_frame, fit_stochastic, fractions, frozen_prediction, load_inputs,
    predict_decoder, reference_metrics, report, scored, training_weights, write_json,
)
from fit_evaluate_candidate import evaluate_constraints, make_estimator
from select_unified_ancestry import calculate

EXPERIMENTS = tuple(f"{kind}_{training}" for kind in ("ridge", "stochastic", "extratrees")
                    for training in ("pure", "pedigree"))
PEDIGREE_WEIGHT = .25


def check_disjoint(train_donors, pure_ids, eval_donors, train_ids, eval_ids):
    if not train_donors or not train_donors.issubset(pure_ids):
        raise ValueError("Training pedigrees must use the audited pure-training donor pool")
    if train_donors & eval_donors:
        raise ValueError("Training and evaluation founder donors overlap")
    if train_ids & eval_ids:
        raise ValueError("Training and evaluation node IDs collide")


def load_training(directory, cols, pure, evaluation, eval_donors, root):
    plan_dir = directory / "pedigree/plans/plans"
    audit_path = plan_dir / "training_cohort_audit.json"
    audit = json.loads(audit_path.read_text())
    if (audit.get("cohort_role") != "decoder_training"
            or audit.get("reference_disjoint") is not True
            or audit.get("labels_validated") is not True
            or audit.get("uses_old_mosaic_training_examples") is not False
            or audit.get("development_donor_overlap") != 0
            or audit.get("pure_validation_donor_overlap") != 0):
        raise ValueError("Explicit reference/development-disjoint training audit is required")
    for path, expected in audit["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Training cohort input changed: {path}")
    plan_path = plan_dir / "expected_ancestry.tsv"
    plan = pd.read_csv(plan_path, sep="\t")
    donors = set().union(*family_donors(plan).values())
    check_disjoint(donors, set(pure.sample_id), eval_donors,
                   set(plan.sample_id), set(evaluation.sample_id))
    if donors != set(audit["used_donor_ids"]):
        raise ValueError("Training audit and plan donor IDs differ")
    truth_path = directory / "pedigree/truth/all_samples_ancestry.tsv"
    truth = pd.read_csv(truth_path, sep="\t")
    truth = truth[truth.is_founder.astype(str).str.lower().isin(["false", "0"])]
    features = sorted((directory / "pedigree/features").glob(f"{REPRESENTATION}__*.features.tsv"))
    manifest_path = plan_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("design") != "full" or len(features) != manifest["n_batches"]:
        raise ValueError("Incomplete full-design training cohort")
    training = exact_join(feature_frame(features, cols), truth)
    plan_desc = plan[plan.is_founder.astype(str).str.lower().isin(["false", "0"])]
    if set(plan_desc.sample_id) != set(training.sample_id):
        raise ValueError("Training plan and feature IDs differ")
    matching = training[["sample_id", "family_id", "class"]].merge(
        plan_desc[["sample_id", "family_id", "class"]], on="sample_id", validate="one_to_one")
    if not (matching.family_id_x.eq(matching.family_id_y)
            & matching.class_x.eq(matching.class_y)).all():
        raise ValueError("Training truth and plan family/class assignments differ")
    fractions(training[[f"realized_{s}" for s in SOURCES]], "Training realized truth")
    if (training.family_id.nunique() != audit["n_families"]
            or len(training) != audit["n_descendants"]
            or set(training.family_id) != set(plan_desc.family_id)):
        raise ValueError("Training cohort coverage differs from its audit")
    audit_file = directory / "pilot_inputs/pilot_input_audit.json"
    frozen = json.loads(audit_file.read_text())
    if frozen["candidate_count"] != 1 or frozen["representation_count"] != 1:
        raise ValueError("Training inference must use only frozen v1")
    pilot = root / "results_v3_pedigree_pilot_10/pilot_inputs/pilot_input_audit.json"
    original = json.loads(pilot.read_text())
    if (frozen["frozen_artifact_sha256"][REPRESENTATION]
            != original["frozen_artifact_sha256"][REPRESENTATION]):
        raise ValueError("Training and development FLARE artifacts differ")
    return training, audit, [audit_path, plan_path, manifest_path, truth_path, audit_file, *features]


def fit_matched(kind, x, y, weights):
    if kind == "stochastic":
        return fit_stochastic(x, y, weights)
    config = {
        "ridge": {"type": "ridge", "params": {"alpha": .001, "fit_intercept": False}},
        "extratrees": {"type": "extratrees", "params": {
            "n_estimators": 200, "min_samples_leaf": 10, "max_features": 1.,
        }},
    }
    if kind not in config:
        raise ValueError(f"Unknown decoder: {kind}")
    model = make_estimator(config[kind], SEED)
    return model.fit(x, y, sample_weight=weights)


def assemble_training(pure, pedigree, cols, include_pedigree):
    x, y = pure[cols].to_numpy(float), pure[list(SOURCES)].to_numpy(float)
    if not include_pedigree:
        return x, y, training_weights(pure)
    weights = training_weights(pure, pedigree)
    weights[:len(pure)] *= 2 * (1 - PEDIGREE_WEIGHT)
    weights[len(pure):] *= 2 * PEDIGREE_WEIGHT
    return (np.vstack([x, pedigree[cols].to_numpy(float)]),
            np.vstack([y, pedigree[[f"realized_{s}" for s in SOURCES]].to_numpy(float)]),
            weights)


def unified_metrics(ped, reference, cfg):
    def convert(frame, reference_set):
        out = frame[["sample_id", "family_id"]].copy()
        out["scenario"] = "reference_source" if reference_set else frame["class"]
        for source in SOURCES:
            out[f"pred_{source}"] = frame[f"pred_{source}"]
            out[f"truth_{source}"] = frame[f"realized_{source}"]
        return out
    return calculate(pd.concat([convert(ped, False), convert(reference, True)], ignore_index=True), cfg)


def experiment(root, source, training_dir, mode, out):
    saved, cols, evaluation, pure, reference, real, donors, folds, baseline, provenance = load_inputs(root, source)
    eval_donors = set().union(*donors.values())
    training, cohort_audit, paths = load_training(training_dir, cols, pure, evaluation, eval_donors, root)
    cfg_path = root / "config/unified_selection.json"
    paths.extend([cfg_path, Path(__file__)])
    provenance["input_sha256"].update({str(p): digest(p) for p in paths})
    provenance.update({
        "evaluation": "fixed training donor pool versus previously inspected donor-disjoint pilot",
        "uses_old_mosaic_training_examples": False,
        "training_n_families": cohort_audit["n_families"],
        "training_n_descendants": len(training), "training_design": "full",
        "reference_evaluation": "all validation donors disjoint from all training donors",
        "pedigree_training_weight": PEDIGREE_WEIGHT,
    })
    write_json(out / "input_audit.json", provenance)
    if mode == "baseline":
        ped_pred = frozen_prediction(saved, evaluation[cols].to_numpy(float))
        ref_pred = frozen_prediction(saved, reference[cols].to_numpy(float))
        real_pred = frozen_prediction(saved, real[cols].to_numpy(float))
        parameters = {"decoder": "historical frozen v1 with its original calibrator"}
    else:
        kind, arm = mode.split("_")
        x, y, weights = assemble_training(pure, training, cols, arm == "pedigree")
        model = fit_matched(kind, x, y, weights)
        parameters = ({"kind": "row_stochastic", "alpha": model["alpha"]}
                      if isinstance(model, dict) else model.get_params())
        joblib.dump({
            "decoder": model, "calibrator": None, "sources": SOURCES, "feature_columns": cols,
            "experiment": mode, "development_screen_only": True,
            "uses_old_mosaic_training_examples": False, "input_sha256": provenance["input_sha256"],
        }, out / "development_model.joblib")
        ped_pred = predict_decoder(model, evaluation[cols].to_numpy(float))
        ref_pred = predict_decoder(model, reference[cols].to_numpy(float))
        real_pred = predict_decoder(model, real[cols].to_numpy(float))
    ped = scored(evaluation, ped_pred)
    ped["pedigree_design"] = "reduced"
    ped["reference_disjoint"] = True
    ped["founder_labels_validated"] = True
    ped["evaluation"] = "fixed_donor_disjoint_development_cohort"
    ped.to_csv(out / "development_predictions.tsv", sep="\t", index=False)
    ref = reference_metrics(reference, ref_pred)
    ref.to_csv(out / "reference_predictions.tsv", sep="\t", index=False)
    rules = pd.read_csv(root / "config/real_benchmarks.tsv", sep="\t")
    rules = rules[rules.dataset_id.eq("lonni_real")]
    if len(rules) != 1 or len(real) != 1:
        raise ValueError("Expected one known real regression case")
    passed, details = evaluate_constraints(dict(zip(SOURCES, real_pred[0])),
                                          json.loads(rules.iloc[0].constraints_json))
    delta = ped.groupby("family_id").TV_distance.mean() - baseline.groupby("family_id").TV_distance.mean()
    delta.rename("mean_tv_delta_vs_v1").to_csv(out / "paired_family_deltas.tsv", sep="\t")
    write_json(out / "experiment.json", {
        "experiment": mode, "parameters": parameters, "seed": SEED,
        "calibrator": "original_frozen" if mode == "baseline" else None,
        "pure_training_weight": 1 - PEDIGREE_WEIGHT if mode.endswith("_pedigree") else 1.,
        "pedigree_training_weight": PEDIGREE_WEIGHT if mode.endswith("_pedigree") else 0.,
        "development_screen_only": True, "promoted": False,
    })
    write_json(out / "metrics.json", {
        "experiment": mode, "development_screen_only": True, "promoted": False,
        "pedigree": report(ped), "reference": report(ref),
        "unified": unified_metrics(ped, ref, json.loads(cfg_path.read_text())),
        "paired_family_mean_tv_delta": float(delta.mean()),
        "families_improved": int((delta < -1e-8).sum()),
        "families_worse": int((delta > 1e-8).sum()),
        "real_regression": {"pass": bool(passed), "details": details,
                            "prediction": dict(zip(SOURCES, real_pred[0]))},
    })
    print(f"{mode}: COMPLETE", flush=True)


def summarize(directory):
    modes = ("baseline", *EXPERIMENTS)
    baseline_audit = json.loads((directory / "baseline/input_audit.json").read_text())
    rows = []
    for mode in modes:
        audit = json.loads((directory / mode / "input_audit.json").read_text())
        if audit != baseline_audit:
            raise ValueError(f"{mode}: development inputs changed between experiments")
        m = json.loads((directory / mode / "metrics.json").read_text())
        rows.append({
            "experiment": mode, "mean_tv": m["pedigree"]["mean_tv"],
            "p95_tv": m["pedigree"]["p95_tv"],
            "family_macro_mean_tv": m["pedigree"]["family_macro_mean_tv"],
            "reference_mean_tv": m["reference"]["mean_tv"],
            "balanced_mean_tv": m["unified"]["balanced_mean_tv"],
            "balanced_p95_tv": m["unified"]["balanced_p95_tv"],
            "off_target_mass": m["pedigree"]["off_target_mass"],
            "minor_2_5_recall": m["pedigree"]["minor_2_5_recall"],
            "real_regression_pass": m["real_regression"]["pass"],
            "families_improved": m["families_improved"], "families_worse": m["families_worse"],
            "development_screen_only": True, "promoted": False,
        })
    for kind in ("ridge", "stochastic", "extratrees"):
        left = json.loads((directory / f"{kind}_pure/experiment.json").read_text())
        right = json.loads((directory / f"{kind}_pedigree/experiment.json").read_text())
        if left["parameters"] != right["parameters"] or left["seed"] != right["seed"]:
            raise ValueError(f"{kind}: matched decoder parameters differ between training arms")
    pd.DataFrame(rows).sort_values("family_macro_mean_tv").to_csv(
        directory / "matched_decoder_screen.tsv", sep="\t", index=False)
    print("Matched screen COMPLETE; no model promoted", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--training-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=("baseline", "summarize", *EXPERIMENTS))
    args = parser.parse_args()
    if args.mode == "summarize":
        summarize(args.outdir)
        return
    args.outdir.mkdir(parents=True, exist_ok=False)
    experiment(args.root, args.source, args.training_dir, args.mode, args.outdir)


if __name__ == "__main__":
    main()
