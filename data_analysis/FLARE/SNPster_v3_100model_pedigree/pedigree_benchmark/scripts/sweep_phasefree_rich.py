#!/usr/bin/env python3
"""Bounded richer-feature sweep using only the existing reference and clean training people."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from phasefree_validation import FEATURES, SOURCES, TRUTH, fold_weights, ranking_key, stable_features, summarize
from report_phasefree_accuracy import accuracy_summary
from screen_dosage_models import standardize
from screen_fine_profiles import matrices, predict as frequency_predict, validate_dosages
from screen_model_families import score
from screen_v1_decoders import predict_decoder, simplex_rows
from sweep_phasefree_calibration import (
    MODEL_SHA, checkpoint, configuration, hgb, predict as calibration_predict, read_table,
)


PROJECT = Path(__file__).resolve().parents[2]
CALIBRATION = Path("/calibration")
STACK = Path("/stack")
DOSAGE = Path("/dosage")
PROTOCOL = PROJECT / "config/phasefree_rich_protocol.json"
CONTROL = "stable_control"


def calibration_outputs(hashes):
    return {str(CALIBRATION / Path(path).relative_to("/sweep")): digest
            for path, digest in hashes.items()}


def canonical_features(x, protocol):
    x = np.asarray(x, float)
    settings = protocol["features"]
    a, b, c = (settings[key] for key in ("base", "pca", "fine"))
    if (a != 14 or x.ndim != 2 or x.shape[1] != a + b + c or not len(x)
            or not np.isfinite(x).all()):
        raise ValueError("Invalid finite rich-feature axes")
    base = stable_features(x[:, :a], settings["base_round_decimals"])
    fine = x[:, a + b:]
    if (np.any(fine < -1e-10) or np.any(fine > 1 + 1e-10)
            or not np.allclose(fine.sum(axis=1), 1, atol=1e-8, rtol=0)):
        raise ValueError("Invalid fine-population simplex features")
    return np.column_stack([
        base, np.round(x[:, a:a + b], settings["pca_round_decimals"]),
        np.round(fine, settings["fine_round_decimals"])])


def genotype_features(components, x, protocol):
    validate_dosages(x)
    classifier, frequency = components["classifier"], components["frequency"]
    if classifier["marker_keys"] != frequency["marker_keys"] or x.shape[1] != len(classifier["marker_keys"]):
        raise ValueError("Richer representation requires the frozen complete marker axis")
    scores = classifier["pca"].transform(
        standardize(x, classifier["mean"], classifier["scale"])).astype(float)
    broad, fine, _ = frequency_predict(frequency, x)
    values = np.column_stack([predict_decoder(classifier["head"], scores), broad, scores, fine])
    canonical_features(values, protocol)
    return values


def feature_columns(name, protocol):
    a, b, c = (protocol["features"][key] for key in ("base", "pca", "fine"))
    choices = dict(base=list(range(a)), pca=list(range(a + b)),
                   fine=list(range(a)) + list(range(a + b, a + b + c)),
                   combined=list(range(a + b + c)))
    if name not in choices:
        raise ValueError(f"Unknown rich feature set: {name}")
    return choices[name]


def fit_models(x, y, weights, protocol):
    values = canonical_features(x, protocol)
    y, weights = np.asarray(y, float), np.asarray(weights, float)
    if (y.shape != (len(x), len(SOURCES)) or not np.isfinite(y).all() or np.any(y < 0)
            or not np.allclose(y.sum(axis=1), 1, atol=1e-10, rtol=0)
            or weights.shape != (len(x),) or not np.isfinite(weights).all()
            or np.any(weights <= 0) or not np.isclose(weights.sum(), 1, atol=1e-12, rtol=0)):
        raise ValueError("Invalid rich-model training truth/weights")
    base_columns = feature_columns("base", protocol)
    base = hgb(protocol["head"], protocol).fit(values[:, base_columns], y, sample_weight=weights)
    mean = weights @ values
    scale = np.maximum(np.sqrt(weights @ ((values - mean) ** 2)),
                       protocol["linear"]["minimum_weighted_scale"])
    linear = Ridge(alpha=protocol["linear"]["alpha"], fit_intercept=True).fit(
        (values - mean) / scale, y, sample_weight=weights)
    models = {}
    for name, parameters in protocol["models"].items():
        columns = feature_columns(parameters["features"], protocol)
        model = dict(parameters=parameters, columns=columns, protocol=protocol)
        if parameters["kind"] == "hgb":
            model["head"] = (base if name == CONTROL else
                             hgb(protocol["head"], protocol).fit(values[:, columns], y, sample_weight=weights))
        elif parameters["kind"] == "blend":
            model.update(head=base, columns=base_columns, linear=linear, mean=mean, scale=scale)
        else:
            raise ValueError(f"Unknown rich head: {parameters}")
        models[name] = model
    return models


def predict(model, x):
    values = canonical_features(x, model["protocol"])
    prediction = simplex_rows(model["head"].predict(values[:, model["columns"]]))
    if model["parameters"]["kind"] == "blend":
        other = simplex_rows(model["linear"].predict((values - model["mean"]) / model["scale"]))
        fraction = model["parameters"]["linear_fraction"]
        prediction = (1 - fraction) * prediction + fraction * other
    return prediction


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Use a fresh richer-feature experiment")
    parent = json.loads((CALIBRATION / "experiment.json").read_text())
    provenance.verify_hashes(parent["implementation_sha256"])
    provenance.verify_hashes(parent["input_sha256"])
    provenance.verify_hashes(calibration_outputs(parent["prepared_output_sha256"]))
    checked = json.loads((CALIBRATION / "verification.json").read_text())
    if (checked["status"] != "complete_verified_cv_and_development"
            or checked["cv_completion_sha256"] != adapter.digest(CALIBRATION / "cv_completion.json")
            or checked["development_completion_sha256"] != adapter.digest(CALIBRATION / "development_completion.json")):
        raise ValueError("Require verified standardized calibration folds and control")
    previous_cv = json.loads((CALIBRATION / "cv_completion.json").read_text())
    previous_development = json.loads((CALIBRATION / "development_completion.json").read_text())
    for record in (previous_cv, previous_development):
        if (record["status"] != "complete"
                or record["experiment_sha256"] != adapter.digest(CALIBRATION / "experiment.json")):
            raise ValueError("Incomplete or incompatible calibration source")
    provenance.verify_hashes(calibration_outputs(previous_cv["output_sha256"]))
    control_path = CALIBRATION / "candidates/stable_squared/model.pkl"
    if adapter.digest(control_path) != previous_development["output_sha256"]["/sweep/candidates/stable_squared/model.pkl"]:
        raise ValueError("Verified prior stable control changed")
    protocol = json.loads(PROTOCOL.read_text())
    if protocol["sources"] != SOURCES:
        raise ValueError("Richer sweep source order changed")
    verified = json.loads((DOSAGE / "verification/final_verification.json").read_text())
    if (verified["status"] != "complete_verified"
            or verified.get("n_markers") != 26066):
        raise ValueError("Require independently verified genotype matrices")
    inputs = dict(parent["input_sha256"])
    paths = [CALIBRATION / name for name in (
        "experiment.json", "verification.json", "cv_completion.json", "development_completion.json", "oof_predictions.npz",
        "training.tsv", "training_features.npy", "folds.json", "donors.json",
        "candidates/stable_squared/model.pkl")]
    paths += [DOSAGE / "verification/final_verification.json", DOSAGE / "experiment.json"]
    for chrom in range(1, 23):
        path = DOSAGE / f"chr{chrom}/audit.json"
        audit = json.loads(path.read_text())
        if audit["status"] != "complete" or audit["experiment_sha256"] != adapter.digest(DOSAGE / "experiment.json"):
            raise ValueError(f"Incompatible genotype chromosome: {chrom}")
        provenance.verify_hashes(audit["output_sha256"])
        inputs.update(audit["output_sha256"])
        paths.append(path)
    inputs.update(provenance.fingerprint_files(paths))
    out.mkdir(parents=True, exist_ok=True)
    for name in ("training.tsv", "folds.json", "donors.json", "cohort.json"):
        shutil.copyfile(CALIBRATION / name, out / name)
    (out / "baseline").mkdir()
    shutil.copyfile(CALIBRATION / "baseline/metrics.json", out / "baseline/metrics.json")
    with Path(parent["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    if adapter.digest(parent["original_model_path"]) != MODEL_SHA:
        raise ValueError("Original frozen model changed")
    x, phase_x, axes, phase_rows, markers = matrices(DOSAGE)
    if original["marker_keys"] != markers or not np.array_equal(x[phase_rows], phase_x):
        raise ValueError("Original genotype marker or phase-invariance contract changed")
    table = read_table(STACK / "features.tsv")
    lookup = {(row.group, row.sample_id): i for i, row in enumerate(axes.itertuples())}
    rows = [lookup[key] for key in table[["group", "sample_id"]].itertuples(index=False, name=None)]
    with threadpool_limits(limits=6):
        values = genotype_features(original["feature_models"], x[rows], protocol)
    base_error = float(np.max(np.abs(values[:, :14] - table[FEATURES].to_numpy(float))))
    if base_error > 1e-11:
        raise ValueError(f"Genotype extraction changed the verified base features: {base_error}")
    values[:, :14] = table[FEATURES].to_numpy(float)
    metadata = read_table(out / "training.tsv")
    training_lookup = {row.sample_id: i for i, row in enumerate(table.itertuples())
                       if row.group == "clean_training"}
    if set(training_lookup) != set(metadata.sample_id) or len(metadata) != 799:
        raise ValueError("Richer training features are not the exact original 799 clean people")
    train_x = values[[training_lookup[identity] for identity in metadata.sample_id]]
    np.testing.assert_array_equal(train_x[:, :14], np.load(CALIBRATION / "training_features.npy"))
    reference_ids = set(axes.loc[axes.group.eq("reference"), "sample_id"])
    evaluation_ids = set(table.loc[table.group.isin(["pure", "pedigree", "real"]), "sample_id"])
    if set(metadata.sample_id) & (reference_ids | evaluation_ids):
        raise ValueError("Training/reference/evaluation identity overlap")
    np.save(out / "training_features.npy", train_x)
    np.save(out / "all_features.npy", values)
    table[["group", "sample_id"]].to_csv(out / "feature_axes.tsv", sep="\t", index=False)
    prepared = [out / name for name in ("training.tsv", "training_features.npy", "all_features.npy",
                                      "feature_axes.tsv", "folds.json", "donors.json", "cohort.json",
                                      "baseline/metrics.json")]
    implementations = dict(parent["implementation_sha256"])
    implementations.update(provenance.fingerprint_files([
        Path(__file__), Path(__file__).with_name("verify_phasefree_rich.py"),
        PROJECT / "tests/test_phasefree_rich.py", PROJECT / "run_phasefree_rich_sweep.sh", PROTOCOL]))
    adapter.write_json(out / "experiment.json", dict(
        status="configured", protocol=protocol, targets=parent["targets"], gates=parent["gates"],
        original_model_path=parent["original_model_path"], implementation_sha256=implementations,
        input_sha256=inputs, prepared_output_sha256=provenance.fingerprint_files(prepared),
        feature_base_reproduction_max_error=base_error, n_training_people=799,
        n_reference_people=len(reference_ids), n_features=values.shape[1], n_folds=114,
        current_dataset_only=True, no_new_imputation=True, no_imputed_training=True,
        reference_only_frozen_preprocessing=True, no_evaluation_in_fit=True,
        independent_final_test=False, production_promotion=False))
    print(f"Prepared 6 heads, {values.shape[1]} features, 799 existing training people and114 donor-purged folds.", flush=True)


def run_fold(out, fold_id):
    config = configuration(out)
    split = next(row for row in json.loads((out / "folds.json").read_text()) if row["fold_id"] == fold_id)
    directory = out / "folds" / fold_id
    if checkpoint(directory / "completion.json", out):
        return fold_id
    if directory.exists():
        raise ValueError(f"Inspect incomplete richer-feature fold: {directory}")
    directory.mkdir(parents=True)
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    train = np.flatnonzero(metadata.sample_id.isin(split["train_ids"]))
    test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
    started = time.monotonic()
    predictions, stability = [], {}
    with threadpool_limits(limits=2):
        models = fit_models(x[train], metadata.iloc[train][TRUTH].to_numpy(float),
                            fold_weights(metadata.iloc[train]), config["protocol"])
        for name, model in models.items():
            prediction = predict(model, x[test])
            np.testing.assert_array_equal(prediction, predict(pickle.loads(pickle.dumps(model)), x[test]))
            difference = max(float(np.max(np.abs(prediction - predict(
                model, np.nextafter(x[test], direction))))) for direction in (-np.inf, np.inf))
            if difference > 1e-12:
                raise ValueError(f"One-ULP instability in richer model: {name}, {fold_id}, {difference}")
            predictions.append(prediction)
            stability[name] = difference
    with np.load(CALIBRATION / "oof_predictions.npz") as source:
        if source["sample_ids"].tolist() != metadata.sample_id.tolist():
            raise ValueError("Prior control OOF identities changed")
        prior = source["predictions"][source["candidates"].tolist().index("stable_squared"), test]
    np.testing.assert_allclose(predictions[list(models).index(CONTROL)], prior, atol=1e-12, rtol=0)
    np.savez_compressed(directory / "predictions.npz", predictions=np.asarray(predictions),
                        sample_ids=metadata.iloc[test].sample_id.to_numpy(str), candidates=np.asarray(list(models)))
    adapter.write_json(directory / "completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"), fold=split,
        elapsed_seconds=time.monotonic() - started, max_one_ulp_prediction_difference=stability,
        prior_control_reproduced=True,
        output_sha256=provenance.fingerprint_files([directory / "predictions.npz"])))
    return fold_id


def cross_validate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "cv_completion.json", out):
        return
    splits = json.loads((out / "folds.json").read_text())
    with ProcessPoolExecutor(max_workers=2) as pool:
        tasks = [pool.submit(run_fold, out, split["fold_id"]) for split in splits]
        for index, future in enumerate(as_completed(tasks), 1):
            print(f"Rich CV {index}/{len(splits)}: {future.result()} verified", flush=True)
    metadata = read_table(out / "training.tsv")
    candidates = list(config["protocol"]["models"])
    predictions = np.empty((len(candidates), len(metadata), len(SOURCES)))
    index = {identity: i for i, identity in enumerate(metadata.sample_id)}
    seen = set()
    for split in splits:
        directory = out / "folds" / split["fold_id"]
        if not checkpoint(directory / "completion.json", out):
            raise ValueError("Incomplete richer CV predictions")
        with np.load(directory / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if data["candidates"].tolist() != candidates or seen & set(ids):
                raise ValueError("Duplicate/changed richer OOF coverage")
            predictions[:, [index[identity] for identity in ids]] = data["predictions"]
            seen.update(ids)
    if seen != set(index):
        raise ValueError("Not every training person has one richer OOF prediction")
    reports = {name: summarize(metadata, predictions[i], config["targets"], config["gates"])
               for i, name in enumerate(candidates)}
    order = sorted(candidates, key=lambda name: ranking_key(name, reports[name]))
    selection = dict(ranking=order, finalists=[name for name in order if name != CONTROL][:2],
                     development_metrics_used=False, criteria=config["protocol"]["ranking"])
    np.savez_compressed(out / "oof_predictions.npz", predictions=predictions,
                        sample_ids=metadata.sample_id.to_numpy(str), candidates=np.asarray(candidates))
    adapter.write_json(out / "cv_metrics.json", reports)
    adapter.write_json(out / "selection.json", selection)
    comparison = [dict(candidate=name, cv_rank=order.index(name) + 1,
                       balanced_mean_tv=reports[name]["unified"]["balanced_mean_tv"],
                       pure_mean_tv=reports[name]["pure"]["mean_tv"],
                       pedigree_mean_tv=reports[name]["pedigree"]["mean_tv"],
                       pedigree_within_10pp=reports[name]["pedigree"]["fraction_all_components_within_10pp"],
                       pedigree_within_5pp=reports[name]["pedigree"]["fraction_all_components_within_5pp"])
                  for name in order]
    pd.DataFrame(comparison).to_csv(out / "cv_comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "cv_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        output_sha256=provenance.fingerprint_files([out / name for name in (
            "oof_predictions.npz", "cv_metrics.json", "selection.json", "cv_comparison.tsv")])))
    print(pd.DataFrame(comparison).to_string(index=False), flush=True)


def evaluate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "development_completion.json", out):
        return
    if not checkpoint(out / "cv_completion.json", out):
        raise ValueError("Lock training-only richer CV before development evaluation")
    selection = json.loads((out / "selection.json").read_text())
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    axes, values = read_table(out / "feature_axes.tsv"), np.load(out / "all_features.npy")
    evaluation = axes.group.isin(["pure", "pedigree", "real"]).to_numpy()
    old_array = axes.group.eq("array_imputed").to_numpy()
    if old_array.sum() != 5:
        raise ValueError("Require the same original five array controls")
    with Path(config["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    with (CALIBRATION / "candidates/stable_squared/model.pkl").open("rb") as stream:
        previous_control = pickle.load(stream)
    paths, comparison = [], []
    with threadpool_limits(limits=2):
        models = fit_models(x, metadata[TRUTH].to_numpy(float), metadata.weight_clean.to_numpy(float), config["protocol"])
        for name, model in models.items():
            directory = out / "candidates" / name
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / "model.pkl"
            if path.exists():
                raise ValueError(f"Inspect incomplete richer full fit: {path}")
            bundle = dict(kind="phasefree_rich", sources=SOURCES, marker_keys=original["marker_keys"],
                          feature_models=original["feature_models"], decoder=model,
                          protocol=config["protocol"], experiment_sha256=adapter.digest(out / "experiment.json"),
                          input_contract="Complete pre-QC hard-GT counts in frozen marker order; no missing/DS padding")
            with path.open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            paths.append(path)
            if name not in [CONTROL, *selection["finalists"]]:
                continue
            prediction = predict(model, values[evaluation])
            np.testing.assert_array_equal(prediction, predict(pickle.loads(pickle.dumps(model)), values[evaluation]))
            if name == CONTROL:
                np.testing.assert_allclose(prediction, calibration_predict(
                    previous_control["decoder"], values[evaluation, :14]), atol=1e-12, rtol=0)
            frame = pd.DataFrame(prediction, columns=[f"pred_{source}" for source in SOURCES])
            frame.insert(0, "sample_id", axes.loc[evaluation, "sample_id"].to_numpy())
            frame.insert(0, "dataset", axes.loc[evaluation, "group"].to_numpy())
            frame.to_csv(directory / "predictions.tsv", sep="\t", index=False)
            score(out, name)
            parts = []
            for track in ("pure", "pedigree"):
                part = read_table(directory / f"{track}_predictions_scored.tsv").assign(track=track)
                if track == "pure":
                    part["source_group"] = part.source_group.replace({"EAS": "EAS_SEA", "SEA": "EAS_SEA"})
                parts.append(part)
            combined = pd.concat(parts, ignore_index=True)
            report = summarize(combined, combined[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                               config["targets"], config["gates"])
            array_report, _ = accuracy_summary(
                np.tile([.5, .5, 0, 0, 0, 0, 0], (5, 1)), predict(model, values[old_array]), config["targets"])
            report["old_five_eur_asian_array_controls"] = array_report
            adapter.write_json(directory / "standard_validation.json", report)
            metrics = json.loads((directory / "metrics.json").read_text())
            comparison.append(dict(
                candidate=name, pure_mean_tv=report["pure"]["mean_tv"],
                pedigree_mean_tv=report["pedigree"]["mean_tv"],
                pedigree_within_10pp=report["pedigree"]["fraction_all_components_within_10pp"],
                pedigree_within_5pp=report["pedigree"]["fraction_all_components_within_5pp"],
                array5_mean_tv=array_report["mean_tv"], unchanged_gates_pass=metrics["all_gates_pass"]))
            paths += [directory / name for name in ("predictions.tsv", "metrics.json", "standard_validation.json")]
            paths += list(directory.glob("*scored.tsv"))
    pd.DataFrame(comparison).to_csv(out / "development_comparison.tsv", sep="\t", index=False)
    paths.append(out / "development_comparison.tsv")
    adapter.write_json(out / "development_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        selection_sha256=adapter.digest(out / "selection.json"),
        output_sha256=provenance.fingerprint_files(paths), production_promotion=False))
    print(pd.DataFrame(comparison).to_string(index=False), flush=True)


def array_evaluate(out, arrays):
    config = configuration(out, full=True)
    if checkpoint(out / "array_completion.json", out):
        return
    checked = json.loads((arrays / "verification.json").read_text())
    if (checked["status"] != "complete_verified_seen_donor_stress_test"
            or checked["model_sha256"] != MODEL_SHA
            or checked["experiment_sha256"] != adapter.digest(arrays / "experiment.json")
            or checked["completion_sha256"] != adapter.digest(arrays / "completion.json")):
        raise ValueError("Require matching independently verified shared-array genotypes")
    plan = json.loads((arrays / "experiment.json").read_text())
    with Path(config["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    keys = original["marker_keys"]
    offsets = {key: i for i, key in enumerate(keys)}
    clean, imputed = np.empty((14, len(keys)), dtype=np.int8), np.empty((28, len(keys)), dtype=np.int8)
    seen = set()
    for chrom in range(1, 23):
        directory = arrays / f"chr{chrom}"
        audit = json.loads((directory / "completion.json").read_text())
        if audit["status"] != "complete" or audit["experiment_sha256"] != adapter.digest(arrays / "experiment.json"):
            raise ValueError("Invalid shared-array chromosome provenance")
        provenance.verify_hashes(audit["output_sha256"])
        current = [tuple(key) for key in json.loads((directory / "marker_keys.json").read_text())]
        if seen & set(current):
            raise ValueError("Duplicate shared-array markers")
        seen.update(current)
        indices = [offsets[key] for key in current]
        with np.load(directory / "dosages.npz") as data:
            clean[:, indices], imputed[:, indices] = data["original"], data["array_imputed"]
    metadata = read_table(out / "training.tsv")
    if (seen != set(keys) or len(plan["samples"]) != 14 or len(plan["cases"]) != 28
            or set(plan["masks"]) != {"lonni", "frederik"}
            or set(plan["donor_ids"]) & set(metadata.sample_id)):
        raise ValueError("Shared-array marker/sample/donor contract failed")
    names = [row["sample_id"] for row in plan["samples"]]
    selection = json.loads((out / "selection.json").read_text())
    rows = []
    with threadpool_limits(limits=6):
        clean_features = genotype_features(original["feature_models"], clean, config["protocol"])
        imputed_features = genotype_features(original["feature_models"], imputed, config["protocol"])
    np.savez_compressed(out / "array_features.npz", clean=clean_features, imputed=imputed_features,
                        sample_ids=np.asarray(names), case_ids=np.asarray([row["case_id"] for row in plan["cases"]]))
    with threadpool_limits(limits=2):
        for name in [CONTROL, *selection["finalists"]]:
            with (out / f"candidates/{name}/model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            base, recovered = predict(model["decoder"], clean_features), predict(model["decoder"], imputed_features)
            for i, case in enumerate(plan["cases"]):
                for stage, prediction in (("clean", base[names.index(case["sample_id"])]), ("array_imputed", recovered[i])):
                    rows.append(dict(candidate=name, sample_id=case["sample_id"], mask=case["mask"],
                                     track=case["track"], stage=stage,
                                     **dict(zip(TRUTH, case["truth"])), **dict(zip(SOURCES, prediction))))
    frame = pd.DataFrame(rows)
    comparisons = []
    for (name, mask, track, stage), group in frame.groupby(["candidate", "mask", "track", "stage"]):
        summary, _ = accuracy_summary(group[TRUTH].to_numpy(float), group[SOURCES].to_numpy(float), config["targets"])
        comparisons.append(dict(candidate=name, mask=mask, track=track, stage=stage, **summary))
    frame.to_csv(out / "array_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(comparisons).to_csv(out / "array_comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "array_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        selection_sha256=adapter.digest(out / "selection.json"),
        source_verification_sha256=adapter.digest(arrays / "verification.json"),
        n_unique_people=14, n_conditions=28, no_new_imputation=True,
        output_sha256=provenance.fingerprint_files([out / name for name in (
            "array_predictions.tsv", "array_comparison.tsv", "array_features.npz")]),
        independent_final_test=False, production_promotion=False))


def run(out):
    adapter.write_json(out / "sweep_state.json", dict(status="cross_validating"))
    cross_validate(out)
    adapter.write_json(out / "sweep_state.json", dict(status="fitting_locked_finalists"))
    evaluate(out)
    from verify_phasefree_rich import verify, verify_arrays
    verify(out)
    arrays, deadline = Path("/validation"), time.monotonic() + 48 * 3600
    while not (arrays / "verification.json").exists():
        state = json.loads((arrays / "durable_scheduler.json").read_text())
        if state["status"] == "failed":
            raise RuntimeError("Shared array dependency failed; recover it without changing locked richer models")
        if time.monotonic() >= deadline:
            raise TimeoutError("Shared array verification unavailable after48 hours")
        adapter.write_json(out / "sweep_state.json", dict(
            status="waiting_for_shared_array_verification", development_verified=True, production_promotion=False))
        time.sleep(60)
    adapter.write_json(out / "sweep_state.json", dict(status="evaluating_shared_arrays"))
    array_evaluate(out, arrays)
    verify_arrays(out, arrays)
    adapter.write_json(out / "sweep_state.json", dict(
        status="complete_verified", recommendation_sha256=adapter.digest(out / "recommendation.json"),
        production_promotion=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "smoke", "run"), required=True)
    parser.add_argument("--output", type=Path, default=Path("/sweep"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.output)
    elif args.stage == "smoke":
        split = json.loads((args.output / "folds.json").read_text())[0]
        run_fold(args.output, split["fold_id"])
        path = args.output / "folds" / split["fold_id"] / "completion.json"
        if not checkpoint(path, args.output):
            raise ValueError("Richer smoke fold did not produce a verified checkpoint")
        adapter.write_json(args.output / "smoke_passed.json", dict(
            status="complete", experiment_sha256=adapter.digest(args.output / "experiment.json"),
            fold_completion_sha256=adapter.digest(path)))
    else:
        try:
            run(args.output)
        except (ValueError, RuntimeError, OSError, AssertionError) as error:
            adapter.write_json(args.output / "sweep_state.json", dict(
                status="failed", error_type=type(error).__name__, error=str(error), production_promotion=False))
            raise
