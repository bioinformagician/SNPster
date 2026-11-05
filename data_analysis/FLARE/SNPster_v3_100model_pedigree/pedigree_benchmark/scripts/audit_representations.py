#!/usr/bin/env python3
"""Diagnose source information, population coverage and mixture decoding."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

from prepare_founders import load_label_metadata, parse_tsv, validate_labels
from screen_v1_decoders import (
    SOURCES, digest, exact_join, feature_frame, fractions,
    frozen_prediction, load_inputs, reference_metrics, report, scored, write_json,
)
from fit_evaluate_candidate import load_truth

PROBE_C = 10.
PROBE_MAX_ITER = 2000


def reproduce(pred, ids, previous, columns, label):
    ids = pd.Index(ids)
    if (ids.has_duplicates or previous.sample_id.duplicated().any()
            or set(ids) != set(previous.sample_id)):
        raise ValueError(f"{label}: saved and recomputed sample IDs differ")
    older = previous.set_index("sample_id").loc[ids, list(columns)].to_numpy(float)
    if not np.allclose(pred, older, atol=1e-10, rtol=0):
        raise ValueError(f"{label}: saved predictions cannot be reproduced")


def collapsed_eas_sea_tv(pred, truth):
    pred = fractions(pred, "Diagnostic predictions")
    truth = fractions(truth, "Diagnostic truth")
    e, s = SOURCES.index("EAS"), SOURCES.index("SEA")
    return (np.abs(pred - truth).sum(axis=1) - np.abs(pred[:, e] - truth[:, e])
            - np.abs(pred[:, s] - truth[:, s])
            + np.abs(pred[:, e] + pred[:, s] - truth[:, e] - truth[:, s])) / 2


def fit_probe(x, labels):
    if set(labels) != {0, 1}:
        raise ValueError("EAS/SEA probe requires both sources in training")
    probe = LogisticRegression(C=PROBE_C, class_weight="balanced", solver="lbfgs",
                               max_iter=PROBE_MAX_ITER, tol=1e-6)
    probe.fit(x, labels)
    if probe.n_iter_.max() >= PROBE_MAX_ITER:
        raise RuntimeError("Source-separation probe did not converge")
    return probe


def binary_metrics(labels, probability):
    labels = np.asarray(labels, dtype=int)
    if set(labels) != {0, 1}:
        raise ValueError("Aggregate EAS/SEA probe evaluation requires both sources")
    pred = probability >= .5
    return {
        "n": len(labels), "n_eas": int((labels == 0).sum()), "n_sea": int((labels == 1).sum()),
        "balanced_accuracy": float(balanced_accuracy_score(labels, pred)),
        "roc_auc": float(roc_auc_score(labels, probability)),
        "eas_recall": float((~pred[labels == 0]).mean()),
        "sea_recall": float(pred[labels == 1].mean()),
    }


def probe_sources(train, validation, cols, eligible_ids, out):
    train = train[train.source_group.isin(["EAS", "SEA"])].copy()
    validation = validation[validation.source_group.isin(["EAS", "SEA"])].copy()
    if set(train.sample_id) & set(validation.sample_id):
        raise ValueError("Source probe training and validation IDs overlap")
    train_y = train.source_group.eq("SEA").to_numpy(int)
    val_y = validation.source_group.eq("SEA").to_numpy(int)
    model = fit_probe(train[cols].to_numpy(float), train_y)
    probability = model.predict_proba(validation[cols].to_numpy(float))[:, 1]
    population_p = np.full(len(validation), np.nan)
    val_populations = validation.fine.to_numpy()
    folds = []
    for population in sorted(set(val_populations)):
        tr = train.fine.ne(population).to_numpy()
        va = np.flatnonzero(val_populations == population)
        fitted = fit_probe(train.loc[tr, cols].to_numpy(float), train_y[tr])
        population_p[va] = fitted.predict_proba(validation.iloc[va][cols].to_numpy(float))[:, 1]
        folds.append({
            "heldout_population": population,
            "n_same_population_training_donors_removed": int((~tr).sum()),
            "training_donor_ids": train.loc[tr, "sample_id"].tolist(),
            "training_populations": sorted(set(train.loc[tr, "fine"])),
            "validation_donor_ids": validation.iloc[va].sample_id.tolist(),
        })
    if not np.isfinite(population_p).all():
        raise ValueError("Missing leave-population-out probe predictions")
    validation["probe_p_sea"] = probability
    validation["population_holdout_p_sea"] = population_p
    validation["probe_correct"] = (probability >= .5) == val_y
    validation["population_holdout_correct"] = (population_p >= .5) == val_y
    validation["pedigree_eligible"] = validation.sample_id.isin(eligible_ids)
    validation[["sample_id", "source_group", "fine", "pedigree_eligible",
                "probe_p_sea", "population_holdout_p_sea", "probe_correct",
                "population_holdout_correct"]].to_csv(out / "eas_sea_probe_predictions.tsv", sep="\t", index=False)
    write_json(out / "population_probe_folds.json", folds)
    joblib.dump({"classifier": model, "feature_columns": cols, "positive_class": "SEA",
                 "diagnostic_only": True, "not_an_ancestry_proportion_decoder": True},
                out / "source_probe.joblib")
    eligible = validation.pedigree_eligible.to_numpy()
    if set(val_y[eligible]) != {0, 1}:
        raise ValueError("Pedigree-eligible EAS/SEA probe subset must cover both sources")
    return {
        "hyperparameters": {"C": PROBE_C, "class_weight": "balanced",
                            "solver": "lbfgs", "max_iter": PROBE_MAX_ITER},
        "n_train_eas": int((train_y == 0).sum()), "n_train_sea": int((train_y == 1).sum()),
        "validation": binary_metrics(val_y, probability),
        "validation_population_holdout": binary_metrics(val_y, population_p),
        "pedigree_eligible_validation": binary_metrics(val_y[eligible], probability[eligible]),
        "pedigree_eligible_population_holdout": binary_metrics(val_y[eligible], population_p[eligible]),
        "interpretation": (
            "Fixed pure-source classification probe, not admixture accuracy or a proof of "
            "information absence. Population holdout removes labels from that population, "
            "but the already-frozen representation was not retrained."
        ),
    }


def expected_mixtures(plan, founder_features, cols, saved, sample_ids):
    if founder_features.sample_id.duplicated().any():
        raise ValueError("Duplicate founder feature IDs")
    by_donor = founder_features.set_index("sample_id")
    donor_predictions = dict(zip(
        founder_features.sample_id,
        frozen_prediction(saved, founder_features[cols].to_numpy(float)),
    ))
    x_nodes, pred_nodes = {}, {}
    for node in plan.sort_values("generation", kind="stable").itertuples():
        if str(node.is_founder).lower() in ("true", "1"):
            if node.founder_donor not in by_donor.index:
                raise ValueError(f"Missing founder feature vector: {node.founder_donor}")
            x = by_donor.loc[node.founder_donor, cols].to_numpy(float)
            x_nodes[node.sample_id] = x
            pred_nodes[node.sample_id] = donor_predictions[node.founder_donor]
        else:
            if node.father_id not in x_nodes or node.mother_id not in x_nodes:
                raise ValueError(f"Missing ancestor features: {node.sample_id}")
            x_nodes[node.sample_id] = (x_nodes[node.father_id] + x_nodes[node.mother_id]) / 2
            pred_nodes[node.sample_id] = (pred_nodes[node.father_id] + pred_nodes[node.mother_id]) / 2
    return (np.vstack([x_nodes[s] for s in sample_ids]),
            np.vstack([pred_nodes[s] for s in sample_ids]))


def population_errors(dataset, saved, cols, name, eligible_ids, used_ids, train_counts):
    pred = frozen_prediction(saved, dataset[cols].to_numpy(float))
    frame = reference_metrics(dataset, pred)
    frame["dataset"] = name
    frame["fine"] = dataset.fine.to_numpy()
    frame["source_group"] = dataset.source_group.to_numpy()
    frame["decoded_source"] = np.asarray(SOURCES)[pred.argmax(axis=1)]
    frame["source_correct"] = frame.decoded_source.eq(frame.source_group)
    frame["pedigree_eligible"] = frame.sample_id.isin(eligible_ids)
    frame["used_as_founder"] = frame.sample_id.isin(used_ids)
    frame["n_training_donors_same_population"] = frame.fine.map(train_counts).fillna(0).astype(int)
    return frame


def audit(root, source, index, out):
    # Reuse the validated common pilot inputs without modifying old audit scripts.
    _, _, common_ped, _, _, _, donors, _, _, common = load_inputs(root, source)
    manifest_path = root / "config/frozen_pilot_10.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    if len(manifest) != 10 or manifest.representation.nunique() != 10:
        raise ValueError("Expected the ten distinct frozen pilot representations")
    if not 0 <= index < len(manifest):
        raise ValueError("Representation index outside manifest")
    row = manifest.iloc[index]
    cid, rep = row.candidate_id, row.representation
    pilot = root / "results_v3_pedigree_pilot_10"
    historical = root / "results_v3_sweep_100"
    candidate = Path(row.candidate_dir)
    paths = [manifest_path, Path(__file__)]
    input_audit = json.loads((pilot / "pilot_inputs/pilot_input_audit.json").read_text())
    expected = input_audit["frozen_artifact_sha256"][cid]
    for name in ("candidate.json", "metrics.json", "model.joblib", "predictions.tsv"):
        path = candidate / name
        if digest(path) != expected[name]:
            raise ValueError(f"{cid}: frozen {name} changed")
        paths.append(path)
    for name, path in (("flare_model", Path(row.model_path)), ("model_config", Path(row.model_config))):
        if digest(path) != expected[name]:
            raise ValueError(f"{cid}: FLARE artifact changed")
        paths.append(path)
    saved = joblib.load(candidate / "model.joblib")
    if (tuple(saved["sources"]) != SOURCES
            or saved["candidate"]["candidate_id"] != cid
            or saved["candidate"]["representation"] != rep):
        raise ValueError("Frozen model identity/source ordering mismatch")
    if "n_jobs" in saved["decoder"].get_params():
        saved["decoder"].set_params(n_jobs=2)
    cols = saved["feature_columns"]
    labels = load_label_metadata(source / "ancestry_role_audit.tsv")

    def dataset(name):
        fp = historical / f"features/dev/{rep}__{name}.features.tsv"
        tp = source / f"{name}_metadata.tsv"
        validate_labels(parse_tsv(tp), labels, name)
        paths.extend([fp, tp])
        return exact_join(feature_frame([fp], cols), load_truth(tp, "source_group"))

    pure, validation = dataset("source_train"), dataset("source_validation")
    if set(pure.sample_id) & set(validation.sample_id):
        raise ValueError("Pure training and validation donors overlap")
    feature_paths = sorted((pilot / "pedigree/features").glob(f"{rep}__*.features.tsv"))
    paths.extend(feature_paths)
    truth = common_ped.drop(columns=[c for c in common_ped if c.startswith("anc_")])
    ped = exact_join(feature_frame(feature_paths, cols), truth)
    pred = frozen_prediction(saved, ped[cols].to_numpy(float))
    ped_saved_path = pilot / f"pedigree/evaluations/{cid}.pedigree_predictions.tsv"
    previous = pd.read_csv(ped_saved_path, sep="\t")
    paths.append(ped_saved_path)
    reproduce(pred, ped.sample_id, previous, [f"pred_{s}" for s in SOURCES], cid)
    val_pred = frozen_prediction(saved, validation[cols].to_numpy(float))
    old = pd.read_csv(candidate / "predictions.tsv", sep="\t")
    old = old[old.dataset_id.eq("source_validation")]
    reproduce(val_pred, validation.sample_id, old, SOURCES, cid + " reference")
    eligible_ids = set(pd.read_csv(pilot / "pedigree/inputs/eligible_heldout_founders.tsv", sep="\t")["sample"])
    used_ids = set().union(*donors.values())
    if used_ids & set(pure.sample_id):
        raise ValueError("Probe training overlaps pedigree founder donors")
    counts = pure.fine.value_counts().to_dict()
    population_frames = [
        population_errors(data, saved, cols, name, eligible_ids, used_ids, counts)
        for data, name in ((pure, "source_train"), (validation, "source_validation"))
    ]
    populations = pd.concat(population_frames, ignore_index=True)
    populations.to_csv(out / "pure_predictions.tsv", sep="\t", index=False)
    grouped = populations.groupby(["dataset", "source_group", "fine"]).agg(
        n=("sample_id", "size"), mean_tv=("TV_distance", "mean"),
        mean_true_source_prediction=("TV_distance", lambda x: 1 - x.mean()),
        dominant_accuracy=("source_correct", "mean"),
        n_training_donors=("n_training_donors_same_population", "first"),
        n_pedigree_eligible=("pedigree_eligible", "sum"),
        n_used_founders=("used_as_founder", "sum")).reset_index()
    grouped.to_csv(out / "population_errors.tsv", sep="\t", index=False)
    validation_predictions = population_frames[1]
    confusion = pd.crosstab(validation_predictions.source_group,
                            validation_predictions.decoded_source).reindex(
                                index=SOURCES, columns=SOURCES, fill_value=0)
    confusion.to_csv(out / "validation_confusion.tsv", sep="\t")
    validation_predictions[validation_predictions.used_as_founder].sort_values(
        "TV_distance", ascending=False).to_csv(out / "founder_errors.tsv", sep="\t", index=False)
    probe = probe_sources(pure, validation, cols, eligible_ids, out)
    plan_path = pilot / "pedigree/plans/plans/expected_ancestry.tsv"
    plan = pd.read_csv(plan_path, sep="\t")
    x_ideal, founder_average = expected_mixtures(plan, validation, cols, saved, ped.sample_id)
    ideal_pred = frozen_prediction(saved, x_ideal)
    expected_truth = fractions(ped[[f"expected_{s}" for s in SOURCES]], "Expected ancestry")
    realized = fractions(ped[[f"realized_{s}" for s in SOURCES]], "Realized ancestry")
    mixed = scored(ped, pred)
    mixed["feature_tv_vs_expected_parent_average"] = np.abs(ped[cols].to_numpy(float) - x_ideal).sum(axis=1) / 2
    mixed["ideal_decoder_tv_vs_expected_truth"] = np.abs(ideal_pred - expected_truth).sum(axis=1) / 2
    mixed["founder_prediction_average_tv_vs_expected_truth"] = np.abs(founder_average - expected_truth).sum(axis=1) / 2
    mixed["mixture_consistency_gap_tv"] = np.abs(ideal_pred - founder_average).sum(axis=1) / 2
    mixed["eas_sea_collapsed_tv"] = collapsed_eas_sea_tv(pred, realized)
    for j, src in enumerate(SOURCES):
        mixed[f"ideal_pred_{src}"] = ideal_pred[:, j]
        mixed[f"founder_average_pred_{src}"] = founder_average[:, j]
    mixed.to_csv(out / "mixture_diagnostics.tsv", sep="\t", index=False)
    source_centroids = pure.groupby("source_group")[cols].mean().loc[list(SOURCES)]
    source_centroids.to_csv(out / "training_source_centroids.tsv", sep="\t")
    centroid_tv = np.abs(source_centroids.loc["EAS"].to_numpy()
                         - source_centroids.loc["SEA"].to_numpy()).sum() / 2
    nearest = []
    train_x = pure[cols].to_numpy(float)
    for v in validation[validation.sample_id.isin(used_ids)].itertuples():
        x = validation.loc[validation.sample_id.eq(v.sample_id), cols].to_numpy(float)[0]
        distances = np.abs(train_x - x).sum(axis=1) / 2
        for rank, j in enumerate(np.argsort(distances, kind="stable")[:5], 1):
            donor = pure.iloc[j]
            nearest.append({
                "validation_donor": v.sample_id, "true_source": v.source_group, "fine": v.fine,
                "rank": rank, "training_donor": donor.sample_id,
                "training_source": donor.source_group, "training_fine": donor.fine,
                "feature_tv": float(distances[j]),
            })
    pd.DataFrame(nearest).to_csv(out / "founder_nearest_training_donors.tsv", sep="\t", index=False)
    decoded = np.asarray(SOURCES)[val_pred.argmax(axis=1)]
    binary_frozen_accuracy = np.mean([
        (decoded[validation.source_group.eq(src).to_numpy()] == src).mean()
        for src in ("EAS", "SEA")
    ])
    raw_pred = frozen_prediction(saved, ped[cols].to_numpy(float), calibrate=False)
    summary = {
        "candidate_id": cid, "representation": rep, "n_features": len(cols),
        "diagnostic_only": True, "independent_final_test": False, "promoted": False,
        "pedigree": report(mixed), "uncalibrated_pedigree": report(scored(ped, raw_pred)),
        "pure_validation": report(validation_predictions),
        "pure_validation_by_source": {
            str(src): report(group) for src, group in validation_predictions.groupby("source_group")
        },
        "frozen_eas_sea_validation_balanced_accuracy": float(binary_frozen_accuracy),
        "source_probe": probe, "training_eas_sea_centroid_tv": float(centroid_tv),
        "mean_feature_tv_vs_expected_parent_average": float(mixed.feature_tv_vs_expected_parent_average.mean()),
        "mean_ideal_decoder_tv_vs_expected_truth": float(mixed.ideal_decoder_tv_vs_expected_truth.mean()),
        "mean_founder_prediction_average_tv_vs_expected_truth":
            float(mixed.founder_prediction_average_tv_vs_expected_truth.mean()),
        "mean_mixture_consistency_gap_tv": float(mixed.mixture_consistency_gap_tv.mean()),
        "mean_eas_sea_collapsed_pedigree_tv": float(mixed.eas_sea_collapsed_tv.mean()),
        "notes": [
            "Source classification performance does not establish admixture identifiability.",
            "Failure of this fixed linear probe does not prove absence of feature information.",
            "Expected parent averages ignore realized recombination; they are diagnostics, not predictors.",
            "Mixture consistency gap is not an additive causal decomposition of prediction error.",
            "Population-held-out probes exclude same-population labels, not that population from the frozen representation.",
        ],
    }
    write_json(out / "metrics.json", summary)
    write_json(out / "input_audit.json", {
        "candidate_id": cid, "representation": rep, "diagnostic_only": True,
        "common_input_sha256": {
            **common["input_sha256"], str(manifest_path): digest(manifest_path),
            str(Path(__file__)): digest(Path(__file__)),
            str(Path(__file__).with_name("prepare_founders.py")):
                digest(Path(__file__).with_name("prepare_founders.py")),
        },
        "representation_input_sha256": {str(p): digest(p) for p in sorted(set(paths))},
        "reference_disjoint_pedigree": True, "frozen_predictions_reproduced": True,
    })
    print(f"{index}: {rep} audit COMPLETE", flush=True)


def summarize(directory):
    rows = []
    common = None
    for index in range(10):
        sub = directory / f"rep_{index:02d}"
        audit = json.loads((sub / "input_audit.json").read_text())
        if common is None:
            common = audit["common_input_sha256"]
        elif common != audit["common_input_sha256"]:
            raise ValueError("Common source/pedigree inputs changed between audits")
        m = json.loads((sub / "metrics.json").read_text())
        probe = m["source_probe"]
        rows.append({
            "representation": m["representation"], "n_features": m["n_features"],
            "pedigree_mean_tv": m["pedigree"]["mean_tv"],
            "reference_mean_tv": m["pure_validation"]["mean_tv"],
            "eas_reference_mean_tv": m["pure_validation_by_source"]["EAS"]["mean_tv"],
            "sea_reference_mean_tv": m["pure_validation_by_source"]["SEA"]["mean_tv"],
            "frozen_eas_sea_balanced_accuracy": m["frozen_eas_sea_validation_balanced_accuracy"],
            "probe_eas_sea_auc": probe["validation"]["roc_auc"],
            "probe_eas_sea_balanced_accuracy": probe["validation"]["balanced_accuracy"],
            "population_holdout_probe_auc": probe["validation_population_holdout"]["roc_auc"],
            "eligible_probe_balanced_accuracy": probe["pedigree_eligible_validation"]["balanced_accuracy"],
            "eligible_population_holdout_probe_balanced_accuracy":
                probe["pedigree_eligible_population_holdout"]["balanced_accuracy"],
            "feature_vs_parent_average_tv": m["mean_feature_tv_vs_expected_parent_average"],
            "ideal_decoder_tv": m["mean_ideal_decoder_tv_vs_expected_truth"],
            "averaged_founder_prediction_tv": m["mean_founder_prediction_average_tv_vs_expected_truth"],
            "mixture_consistency_gap_tv": m["mean_mixture_consistency_gap_tv"],
            "collapsed_eas_sea_pedigree_tv": m["mean_eas_sea_collapsed_pedigree_tv"],
            "diagnostic_only": True, "promoted": False,
        })
    table = pd.DataFrame(rows)
    if table.representation.nunique() != 10:
        raise ValueError("Audit does not cover ten distinct representations")
    table.to_csv(directory / "representation_diagnostics.tsv", sep="\t", index=False)
    print("Ten representation audits COMPLETE; no models promoted", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--outdir", required=True, type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--index", type=int)
    group.add_argument("--summarize", action="store_true")
    args = parser.parse_args()
    if args.summarize:
        summarize(args.outdir)
    else:
        args.outdir.mkdir(parents=True, exist_ok=False)
        audit(args.root, args.source, args.index, args.outdir)


if __name__ == "__main__":
    main()
