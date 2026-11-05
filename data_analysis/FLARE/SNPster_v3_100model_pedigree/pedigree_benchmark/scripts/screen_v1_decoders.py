#!/usr/bin/env python3
"""Audit corrected pedigrees and screen decoders with founder-purged family folds."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge

from evaluate_finalist import SOURCES, decode, metrics, simplex_rows
from prepare_founders import mosaic_donor_ids, read_reference_panel

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from fit_evaluate_candidate import evaluate_constraints, load_truth

EXPERIMENTS = (
    "ridge_pure", "ridge_mixed", "stochastic_pure", "stochastic_mixed",
    "extratrees_mixed",
)
SEED = 20261013
REPRESENTATION = "v1_frozen_k18"


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fractions(values, label):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 2 or not len(values) or not np.isfinite(values).all()
            or np.any(values < -1e-8)
            or not np.allclose(values.sum(axis=1), 1, atol=1e-6)):
        raise ValueError(f"{label}: expected finite, nonnegative unit-sum vectors")
    return values


def family_donors(plan):
    founders = plan[plan.is_founder.astype(str).str.lower().isin(["true", "1"])]
    if founders.empty or founders[["founder_donor", "founder_source"]].isna().any().any():
        raise ValueError("Missing founder assignments")
    if founders.groupby("founder_donor").founder_source.nunique().max() != 1:
        raise ValueError("Donor has inconsistent ancestry labels")
    if founders.duplicated(["family_id", "founder_donor"]).any():
        raise ValueError("Donor reused within a family")
    donors = {str(k): set(g.founder_donor) for k, g in founders.groupby("family_id")}
    if set(donors) != set(plan.family_id):
        raise ValueError("Missing founders for a family")
    return donors


def purged_folds(frame, donors):
    if set(frame.family_id) != set(donors):
        raise ValueError("Prediction and founder families differ")
    result = []
    for family in sorted(donors):
        train_families = [f for f in sorted(donors) if not donors[f] & donors[family]]
        if not train_families:
            raise ValueError(f"{family}: no founder-disjoint training families")
        train = np.flatnonzero(frame.family_id.isin(train_families))
        test = np.flatnonzero(frame.family_id.eq(family))
        train_donors = set().union(*(donors[f] for f in train_families))
        if train_donors & donors[family] or not len(train) or not len(test):
            raise ValueError(f"{family}: invalid founder-purged split")
        result.append((family, train, test, {
            "heldout_family": family, "test_donors": sorted(donors[family]),
            "training_donors": sorted(train_donors),
            "training_families": train_families,
            "n_train": len(train), "n_test": len(test), "donor_overlap": 0,
        }))
    return result


def reference_indices(reference, training_donors, mixed):
    mask = (~reference.sample_id.isin(training_donors)
            if mixed else np.ones(len(reference), dtype=bool))
    return np.flatnonzero(mask)


def feature_frame(paths, cols):
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t", dtype={"sample_id": str})
        if "sample_id" not in frame or not set(cols).issubset(frame):
            raise ValueError(f"{path}: feature schema differs from frozen v1")
        frames.append(frame[["sample_id", *cols]])
    if not frames:
        raise ValueError("No feature files")
    frame = pd.concat(frames, ignore_index=True)
    if frame.sample_id.isna().any() or frame.sample_id.duplicated().any():
        raise ValueError("Missing/duplicate feature IDs")
    fractions(frame[cols], "FLARE features")
    return frame


def exact_join(features, truth):
    frame = features.merge(truth, on="sample_id", how="outer",
                           validate="one_to_one", indicator=True)
    if not frame._merge.eq("both").all():
        raise ValueError("Feature and truth sample IDs differ")
    return frame.drop(columns="_merge")


def frozen_prediction(saved, x, calibrate=True):
    pred = decode(saved["decoder"], x)
    if calibrate and saved["calibrator"] is not None:
        pred = simplex_rows(saved["calibrator"].predict(pred))
    return fractions(pred, "Frozen predictions")


def scored(frame, pred):
    out = frame.copy()
    pred = fractions(pred, "Decoder predictions")
    y = fractions(out[[f"realized_{s}" for s in SOURCES]], "Realized truth")
    for j, source in enumerate(SOURCES):
        out[f"pred_{source}"] = pred[:, j]
    out["TV_distance"] = np.abs(pred - y).sum(axis=1) / 2
    out["dominant_correct"] = pred.argmax(axis=1) == y.argmax(axis=1)
    out["off_target_mass"] = (pred * (y < .01)).sum(axis=1)
    return out


def report(frame):
    result = metrics(frame)
    result["family_macro_mean_tv"] = float(frame.groupby("family_id").TV_distance.mean().mean())
    result["class_macro_mean_tv"] = float(frame.groupby("class").TV_distance.mean().mean())
    result["by_class"] = {str(k): metrics(g) for k, g in frame.groupby("class")}
    y = frame[[f"realized_{s}" for s in SOURCES]].to_numpy(float)
    p = frame[[f"pred_{s}" for s in SOURCES]].to_numpy(float)
    minor = (y >= .02) & (y <= .05)
    result["minor_2_5_n"] = int(minor.sum())
    result["minor_2_5_recall"] = float((p[minor] >= .01).mean()) if minor.any() else None
    result["absent_component_calls_5pct"] = int((p[y < .01] >= .05).sum())
    return result


def load_inputs(root, source):
    pilot = root / "results_v3_pedigree_pilot_10"
    historical = root / "results_v3_sweep_100"
    candidate = historical / "stage3_full/candidate_v1_frozen_k18"
    paths = []

    def tracked(path):
        paths.append(path)
        return path

    audit = json.loads(tracked(pilot / "pedigree/inputs/eligible_heldout_founders.tsv.summary.json").read_text())
    if (audit.get("strict_reference_disjoint") is not True
            or audit.get("labels_validated") is not True
            or audit.get("mosaic_training_exclusion_checked") is not True
            or audit.get("flare_reference_overlap_count") != 0):
        raise ValueError("Corrected founder audit is required")
    manifest = json.loads(tracked(pilot / "pedigree/plans/plans/manifest.json").read_text())
    if manifest.get("design") != "reduced":
        raise ValueError("This screen expects the explicitly reduced pilot")
    pilot_audit = json.loads(tracked(pilot / "pilot_inputs/pilot_input_audit.json").read_text())
    hashes = pilot_audit["frozen_artifact_sha256"][REPRESENTATION]
    for name in ("model.joblib", "candidate.json", "predictions.tsv", "metrics.json"):
        if digest(tracked(candidate / name)) != hashes[name]:
            raise ValueError(f"Frozen {name} changed since pilot")
    for name, suffix in (("flare_model", "final.model"), ("model_config", "config.json")):
        if digest(tracked(historical / f"models/{REPRESENTATION}.{suffix}")) != hashes[name]:
            raise ValueError("FLARE representation changed since pilot")
    saved = joblib.load(candidate / "model.joblib")
    if tuple(saved["sources"]) != SOURCES or saved["candidate"]["candidate_id"] != REPRESENTATION:
        raise ValueError("Wrong frozen baseline")
    saved["decoder"].set_params(n_jobs=2)
    cols = saved["feature_columns"]
    feature_paths = sorted((pilot / "pedigree/features").glob(f"{REPRESENTATION}__*.features.tsv"))
    paths.extend(feature_paths)
    truth = pd.read_csv(tracked(pilot / "pedigree/truth/all_samples_ancestry.tsv"), sep="\t")
    truth = truth[truth.is_founder.astype(str).str.lower().isin(["false", "0"])].copy()
    ped = exact_join(feature_frame(feature_paths, cols), truth)
    plan = pd.read_csv(tracked(pilot / "pedigree/plans/plans/expected_ancestry.tsv"), sep="\t")
    plan_desc = plan[plan.is_founder.astype(str).str.lower().isin(["false", "0"])]
    if set(plan_desc.sample_id) != set(ped.sample_id):
        raise ValueError("Plan and realized truth descendants differ")
    comparison = ped[["sample_id", "family_id", "class"]].merge(
        plan_desc[["sample_id", "family_id", "class"]], on="sample_id", validate="one_to_one")
    if not (comparison.family_id_x.eq(comparison.family_id_y)
            & comparison.class_x.eq(comparison.class_y)).all():
        raise ValueError("Plan and truth family/class disagree")
    donors = family_donors(plan)
    all_donors = set().union(*donors.values())
    eligible = pd.read_csv(tracked(pilot / "pedigree/inputs/eligible_heldout_founders.tsv"), sep="\t")
    if not all_donors.issubset(set(eligible["sample"])):
        raise ValueError("Plan uses ineligible founders")
    authority = pd.read_csv(tracked(source / "ancestry_role_audit.tsv"), sep="\t")
    labels = authority[["sample", "fine", "source_group"]].drop_duplicates()
    if labels["sample"].duplicated().any():
        raise ValueError("Ambiguous authoritative labels")
    label_map = labels.set_index("sample").source_group.to_dict()
    founders = plan[plan.is_founder.astype(str).str.lower().isin(["true", "1"])]
    if any(label_map.get(row.founder_donor) != row.founder_source for row in founders.itertuples()):
        raise ValueError("Plan founders differ from authoritative labels")
    reference_ids = set(read_reference_panel(tracked(source / "source_reference.ref.panel")))
    mosaic_ids = mosaic_donor_ids(tracked(
        source / "training_mosaics8/training_mosaics8_metadata.tsv"))
    if all_donors & (reference_ids | mosaic_ids):
        raise ValueError("Plan founders overlap the current FLARE reference or mosaic training")

    def dataset(name):
        fp = tracked(historical / f"features/dev/{REPRESENTATION}__{name}.features.tsv")
        tp = tracked(source / f"{name}_metadata.tsv")
        meta = pd.read_csv(tp, sep="\t")
        if any(label_map.get(row.sample) != row.source_group for row in meta.itertuples()):
            raise ValueError(f"{name}: authoritative ancestry labels differ")
        return exact_join(feature_frame([fp], cols), load_truth(tp, "source_group"))

    pure = dataset("source_train")
    validation = dataset("source_validation")
    if set(pure.sample_id) & (all_donors | set(validation.sample_id)):
        raise ValueError("Pure decoder training overlaps pedigree/reference validation donors")
    reference = validation
    if set(reference.source_group) != set(SOURCES):
        raise ValueError("Pure-validation data does not cover all sources")
    real = feature_frame([tracked(historical / f"features/dev/{REPRESENTATION}__lonni_real.features.tsv")], cols)
    if set(real.sample_id) & (all_donors | set(pure.sample_id) | set(validation.sample_id)):
        raise ValueError("Real regression sample overlaps decoder development donors")
    baseline = scored(ped, frozen_prediction(saved, ped[cols].to_numpy(float)))
    published = pd.read_csv(tracked(pilot / "pedigree/evaluations/v1_frozen_k18.pedigree_predictions.tsv"), sep="\t")
    check = baseline[["sample_id", *[f"pred_{s}" for s in SOURCES]]].merge(
        published[["sample_id", *[f"pred_{s}" for s in SOURCES]]],
        on="sample_id", how="outer", validate="one_to_one", indicator=True)
    if (not check._merge.eq("both").all()
            or not np.allclose(check[[f"pred_{s}_x" for s in SOURCES]],
                               check[[f"pred_{s}_y" for s in SOURCES]], atol=1e-10, rtol=0)):
        raise ValueError("Frozen pedigree predictions cannot be reproduced")
    historical_predictions = pd.read_csv(candidate / "predictions.tsv", sep="\t")
    old = historical_predictions[historical_predictions.dataset_id.eq("source_validation")]
    refcheck = validation[["sample_id"]].merge(old, on="sample_id", how="outer", validate="one_to_one")
    if len(refcheck) != len(validation) or refcheck[list(SOURCES)].isna().any().any():
        raise ValueError("Historical validation IDs differ")
    refpred = frozen_prediction(saved, validation[cols].to_numpy(float))
    old = old.set_index("sample_id").loc[validation.sample_id]
    if not np.allclose(refpred, old[list(SOURCES)], atol=1e-10, rtol=0):
        raise ValueError("Reused reference features do not reproduce frozen predictions")
    folds = purged_folds(ped, donors)
    tracked(root / "config/real_benchmarks.tsv")
    tracked(Path(__file__).resolve())
    tracked(Path(__file__).resolve().with_name("evaluate_finalist.py"))
    tracked(root / "scripts/fit_evaluate_candidate.py")
    provenance = {
        "representation": REPRESENTATION, "pedigree_design": "reduced",
        "reference_disjoint": True, "development_screen_only": True,
        "production_domain_validated": False, "independent_final_test": False,
        "evaluation": "leave-one-family-out with all shared-founder families purged",
        "n_families": len(donors), "n_descendants": len(ped),
        "n_pure_train": len(pure), "n_reference": len(reference),
        "reference_ids_used_as_pedigree_founders": sorted(set(validation.sample_id) & all_donors),
        "reference_evaluation": "average only fold models not trained on that validation donor",
        "reference_per_source": reference.source_group.value_counts().to_dict(),
        "input_sha256": {str(p): digest(p) for p in sorted(set(paths))},
    }
    return saved, cols, ped, pure, reference, real, donors, folds, baseline, provenance


def audit_run(data, out):
    saved, cols, ped, pure, reference, real, donors, folds, baseline, provenance = data
    write_json(out / "input_audit.json", provenance)
    write_json(out / "folds.json", [fold[3] for fold in folds])
    write_json(out / "baseline_metrics.json", report(baseline))
    baseline.to_csv(out / "baseline_predictions.tsv", sep="\t", index=False)
    raw = scored(ped, frozen_prediction(saved, ped[cols].to_numpy(float), calibrate=False))
    write_json(out / "uncalibrated_baseline_metrics.json", report(raw))
    raw.to_csv(out / "uncalibrated_baseline_predictions.tsv", sep="\t", index=False)
    families = baseline.groupby(["family_id", "class"]).agg(
        n=("sample_id", "size"), mean_tv=("TV_distance", "mean"),
        max_tv=("TV_distance", "max")).reset_index()
    families["donors"] = families.family_id.map(lambda f: ",".join(sorted(donors[f])))
    families.sort_values("mean_tv", ascending=False).to_csv(out / "family_errors.tsv", sep="\t", index=False)
    rows = []
    for donor in sorted(set().union(*donors.values())):
        associated = [f for f in donors if donor in donors[f]]
        sub = baseline[baseline.family_id.isin(associated)]
        other = baseline[~baseline.family_id.isin(associated)]
        rows.append({
            "donor": donor, "associated_families": len(associated),
            "associated_descendants": len(sub), "mean_tv": float(sub.TV_distance.mean()),
            "family_macro_tv": float(sub.groupby("family_id").TV_distance.mean().mean()),
            "other_mean_tv": float(other.TV_distance.mean()) if len(other) else None,
            "interpretation": "family association, not realized donor contribution or causation",
        })
    pd.DataFrame(rows).sort_values("family_macro_tv", ascending=False).to_csv(
        out / "donor_associations.tsv", sep="\t", index=False)
    baseline.sort_values("TV_distance", ascending=False).head(25).to_csv(out / "worst_descendants.tsv", sep="\t", index=False)
    p = baseline[[f"pred_{s}" for s in SOURCES]].to_numpy(float)
    y = baseline[[f"realized_{s}" for s in SOURCES]].to_numpy(float)
    e, s = SOURCES.index("EAS"), SOURCES.index("SEA")
    merged_error = (np.abs(p - y).sum(axis=1) - np.abs(p[:, e] - y[:, e])
                    - np.abs(p[:, s] - y[:, s])
                    + np.abs(p[:, e] + p[:, s] - y[:, e] - y[:, s])) / 2
    write_json(out / "residual_summary.json", {
        "baseline_mean_tv": float(baseline.TV_distance.mean()),
        "eas_sea_collapsed_diagnostic_mean_tv": float(merged_error.mean()),
        "within_eas_sea_error_mean_tv": float((baseline.TV_distance - merged_error).mean()),
        "ontology_unchanged": True,
        "notes": [
            "Collapsing EAS/SEA is a diagnostic, not a proposed scoring or label change.",
            "Family donor associations include donors whose contribution may be zero in some descendants.",
            "Cross-family founder reuse prevents interpreting families as independent replicates.",
        ],
    })
    plan_path = next(Path(p) for p in provenance["input_sha256"]
                     if p.endswith("/expected_ancestry.tsv"))
    plan = pd.read_csv(plan_path, sep="\t").sort_values("generation")
    feature_by_donor = pd.concat([pure, reference], ignore_index=True).set_index("sample_id")
    if not set().union(*donors.values()).issubset(feature_by_donor.index):
        raise ValueError("Missing independently inferred founder features for mixture diagnostic")
    node_features = {}
    for row in plan.itertuples():
        if str(row.is_founder).lower() in ("true", "1"):
            node_features[row.sample_id] = feature_by_donor.loc[row.founder_donor, cols].to_numpy(float)
        else:
            node_features[row.sample_id] = (
                node_features[row.father_id] + node_features[row.mother_id]) / 2
    ideal_x = np.vstack([node_features[s] for s in ped.sample_id])
    ideal_pred = frozen_prediction(saved, ideal_x)
    ideal = scored(ped, ideal_pred)
    ideal["feature_tv_vs_expected_parent_average"] = (
        np.abs(ped[cols].to_numpy(float) - ideal_x).sum(axis=1) / 2)
    expected = ped[[f"expected_{s}" for s in SOURCES]].to_numpy(float)
    ideal["ideal_decoder_tv_vs_expected_ancestry"] = np.abs(ideal_pred - expected).sum(axis=1) / 2
    ideal.to_csv(out / "ideal_mixture_diagnostic.tsv", sep="\t", index=False)
    write_json(out / "ideal_mixture_summary.json", {
        "actual_feature_vs_expected_parent_average_mean_tv":
            float(ideal.feature_tv_vs_expected_parent_average.mean()),
        "ideal_feature_decoder_vs_expected_ancestry_mean_tv":
            float(ideal.ideal_decoder_tv_vs_expected_ancestry.mean()),
        "ideal_feature_decoder_vs_realized_ancestry_mean_tv":
            float(ideal.TV_distance.mean()),
        "interpretation": "Expected parent averages ignore realized recombination; diagnostic, not a predictor.",
    })


def training_weights(pure, pedigree=None):
    pure_weights = 1 / pure.groupby("source_group").sample_id.transform("size").to_numpy(float)
    pure_weights /= pure_weights.sum()
    if pedigree is None:
        return pure_weights
    counts = pedigree.groupby("family_id").sample_id.transform("size").to_numpy(float)
    ped_weights = 1 / counts
    ped_weights /= ped_weights.sum()
    return np.concatenate([.5 * pure_weights, .5 * ped_weights])


def fit_stochastic(x, y, weights, alpha=1e-5):
    x = fractions(x, "Stochastic training features")
    y = fractions(y, "Stochastic training truth")
    w = np.asarray(weights, dtype=float)
    if w.shape != (len(x),) or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("Invalid training weights")
    w = w / w.sum()
    k, sources = x.shape[1], y.shape[1]
    h = x.T @ (w[:, None] * x) + alpha * np.eye(k)
    b = x.T @ (w[:, None] * y)

    def objective(flat):
        a = flat.reshape(k, sources)
        return .5 * np.sum(a * (h @ a)) - np.sum(a * b)

    def gradient(flat):
        return (h @ flat.reshape(k, sources) - b).ravel()

    constraint = np.kron(np.eye(k), np.ones((1, sources)))
    result = minimize(objective, np.full(k * sources, 1 / sources),
                      jac=gradient, method="SLSQP", bounds=[(0, 1)] * (k * sources),
                      constraints={"type": "eq", "fun": lambda a: constraint @ a - 1,
                                   "jac": lambda a: constraint},
                      options={"maxiter": 2000, "ftol": 1e-10})
    if not result.success:
        raise RuntimeError(f"Stochastic optimization failed: {result.message}")
    matrix = result.x.reshape(k, sources)
    fractions(matrix, "Stochastic decoder rows")
    return {"kind": "row_stochastic", "matrix": matrix, "alpha": alpha,
            "iterations": int(result.nit), "objective": float(result.fun)}


def fit_decoder(experiment, x, y, weights):
    if experiment.startswith("stochastic"):
        return fit_stochastic(x, y, weights)
    if experiment.startswith("ridge"):
        # Weights sum to one, so alpha is relative to mean squared error.
        model = Ridge(alpha=.001, fit_intercept=False)
    elif experiment == "extratrees_mixed":
        model = ExtraTreesRegressor(n_estimators=200, min_samples_leaf=3,
                                    max_features=1., random_state=SEED, n_jobs=2)
    else:
        raise ValueError(f"Unknown experiment: {experiment}")
    return model.fit(x, y, sample_weight=weights)


def predict_decoder(model, x):
    if isinstance(model, dict) and model.get("kind") == "row_stochastic":
        return fractions(x @ model["matrix"], "Stochastic predictions")
    return simplex_rows(model.predict(x))


def train_parts(experiment, cols, pure, pedigree):
    mixed = experiment.endswith("_mixed")
    x = pure[cols].to_numpy(float)
    y = pure[list(SOURCES)].to_numpy(float)
    if mixed:
        x = np.vstack([x, pedigree[cols].to_numpy(float)])
        y = np.vstack([y, pedigree[[f"realized_{s}" for s in SOURCES]].to_numpy(float)])
    return x, y, training_weights(pure, pedigree if mixed else None)


def reference_metrics(frame, pred):
    clean = frame[["sample_id", "source_group"]].copy()
    clean["family_id"] = clean.sample_id
    clean["class"] = "reference_source"
    for source in SOURCES:
        clean[f"realized_{source}"] = frame[source].to_numpy(float)
    return scored(clean, pred)


def experiment_run(data, experiment, root, out):
    saved, cols, ped, pure, reference, real, donors, folds, baseline, provenance = data
    write_json(out / "input_audit.json", provenance)
    write_json(out / "folds.json", [fold[3] for fold in folds])
    pred = np.full((len(ped), len(SOURCES)), np.nan)
    ref_sum = np.zeros((len(reference), len(SOURCES)))
    ref_count = np.zeros(len(reference), dtype=int)
    (out / "fold_models").mkdir()
    for family, train, test, fold in folds:
        x, y, weights = train_parts(experiment, cols, pure, ped.iloc[train])
        model = fit_decoder(experiment, x, y, weights)
        pred[test] = predict_decoder(model, ped.iloc[test][cols].to_numpy(float))
        ref_indices = reference_indices(reference, fold["training_donors"],
                                        experiment.endswith("_mixed"))
        if len(ref_indices):
            ref_sum[ref_indices] += predict_decoder(
                model, reference.iloc[ref_indices][cols].to_numpy(float))
            ref_count[ref_indices] += 1
        fold = {**fold, "reference_test_ids": reference.iloc[ref_indices].sample_id.tolist()}
        joblib.dump({"decoder": model, "calibrator": None, "sources": SOURCES,
                     "feature_columns": cols, "fold": fold, "development_screen_only": True},
                    out / f"fold_models/{family}.joblib")
        print(f"{experiment}: {family}, {len(train)} purged training descendants", flush=True)
    scored_ped = scored(ped, pred)
    scored_ped["evaluation"] = "founder_purged_leave_one_family_out"
    scored_ped["pedigree_design"] = "reduced"
    scored_ped["reference_disjoint"] = True
    scored_ped.to_csv(out / "oof_predictions.tsv", sep="\t", index=False)
    baseline.to_csv(out / "baseline_predictions.tsv", sep="\t", index=False)
    x, y, weights = train_parts(experiment, cols, pure, ped)
    final = fit_decoder(experiment, x, y, weights)
    config = {
        "experiment": experiment, "representation": REPRESENTATION,
        "seed": SEED, "training": "pure_plus_pedigree" if experiment.endswith("_mixed") else "pure_only",
        "pure_source_balanced": True, "pedigree_family_balanced": True,
        "pure_pedigree_weight": [.5, .5] if experiment.endswith("_mixed") else [1., 0.],
        "calibrator": None, "development_screen_only": True, "promoted": False,
        "hyperparameters": ({"alpha": final["alpha"], "row_simplex_constraints": True}
                            if isinstance(final, dict) else final.get_params()),
    }
    write_json(out / "experiment.json", config)
    joblib.dump({"decoder": final, "calibrator": None, "sources": SOURCES,
                 "feature_columns": cols, "experiment": config,
                 "input_sha256": provenance["input_sha256"]}, out / "development_model.joblib")
    if np.any(ref_count == 0):
        raise ValueError("Reference sample has no donor-disjoint fold prediction")
    ref = reference_metrics(reference, ref_sum / ref_count[:, None])
    ref["n_eligible_fold_models"] = ref_count
    ref["evaluation"] = "donor_disjoint_fold_model_average"
    base_ref = reference_metrics(reference, frozen_prediction(saved, reference[cols].to_numpy(float)))
    ref.to_csv(out / "clean_reference_predictions.tsv", sep="\t", index=False)
    real_pred = predict_decoder(final, real[cols].to_numpy(float))
    rb = pd.read_csv(root / "config/real_benchmarks.tsv", sep="\t")
    rules = rb[rb.dataset_id.eq("lonni_real")]
    if len(rules) != 1 or len(real) != 1:
        raise ValueError("Expected exactly one known real regression case")
    passed, details = evaluate_constraints(dict(zip(SOURCES, real_pred[0])),
                                          json.loads(rules.iloc[0].constraints_json))
    family_delta = (scored_ped.groupby("family_id").TV_distance.mean()
                    - baseline.groupby("family_id").TV_distance.mean())
    write_json(out / "metrics.json", {
        "experiment": experiment, "development_screen_only": True, "promoted": False,
        "pedigree_oof": report(scored_ped), "baseline": report(baseline),
        "paired_family_mean_tv_delta": float(family_delta.mean()),
        "families_improved": int((family_delta < -1e-8).sum()),
        "families_worse": int((family_delta > 1e-8).sum()),
        "clean_reference": report(ref), "baseline_clean_reference": report(base_ref),
        "reference_evaluation": "donor-disjoint fold-model average, not the full-data model",
        "real_regression": {"pass": bool(passed), "details": details,
                            "prediction": dict(zip(SOURCES, real_pred[0]))},
    })
    family_delta.rename("mean_tv_delta_vs_v1").to_csv(out / "paired_family_deltas.tsv", sep="\t")
    print(f"{experiment}: COMPLETE", flush=True)


def summarize(out):
    rows = []
    audit = json.loads((out / "audit/input_audit.json").read_text())
    for experiment in EXPERIMENTS:
        metrics_path = out / experiment / "metrics.json"
        if not metrics_path.is_file():
            raise ValueError(f"Missing completed experiment: {experiment}")
        inputs = json.loads((out / experiment / "input_audit.json").read_text())
        if inputs != audit:
            raise ValueError(f"{experiment}: inputs differ from residual audit")
        value = json.loads(metrics_path.read_text())
        ped = value["pedigree_oof"]
        rows.append({
            "experiment": experiment, "mean_tv": ped["mean_tv"], "p95_tv": ped["p95_tv"],
            "family_macro_tv": ped["family_macro_mean_tv"],
            "class_macro_tv": ped["class_macro_mean_tv"],
            "family_delta_vs_v1": value["paired_family_mean_tv_delta"],
            "families_improved": value["families_improved"], "families_worse": value["families_worse"],
            "off_target_mass": ped["off_target_mass"], "minor_2_5_recall": ped["minor_2_5_recall"],
            "clean_reference_mean_tv": value["clean_reference"]["mean_tv"],
            "real_regression_pass": value["real_regression"]["pass"],
            "development_screen_only": True, "promoted": False,
        })
    pd.DataFrame(rows).sort_values("family_macro_tv").to_csv(
        out / "decoder_screen.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mode", choices=["audit", "summarize", *EXPERIMENTS], required=True)
    args = parser.parse_args()
    if args.mode == "summarize":
        summarize(args.outdir)
        return
    args.outdir.mkdir(parents=True, exist_ok=False)
    data = load_inputs(args.root, args.source)
    if args.mode == "audit":
        audit_run(data, args.outdir)
    else:
        experiment_run(data, args.mode, args.root, args.outdir)


if __name__ == "__main__":
    main()
