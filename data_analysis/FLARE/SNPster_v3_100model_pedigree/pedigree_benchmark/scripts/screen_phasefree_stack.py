#!/usr/bin/env python3
"""Training-only calibration of reference-only classification and frequency fits."""

import argparse
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
from screen_dosage_models import predict as pca_predict, standardize
from screen_fine_profiles import matrices, predict as frequency_predict
from screen_model_families import score
from screen_reference_trees import model_for
from screen_supervised_genotypes import training
from screen_v1_decoders import predict_decoder


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
MODES = {
    "phasefree_stack_ridge_clean": dict(kind="ridge", policy="clean"),
    "phasefree_stack_ridge_matched": dict(kind="ridge", policy="matched"),
    "phasefree_stack_boost_clean": dict(kind="boosting", policy="clean"),
}


def features(bundle, x):
    return np.column_stack([
        pca_predict(bundle["classifier"], x),
        frequency_predict(bundle["frequency"], x)[0]])


def decoder_for(kind):
    if kind == "ridge":
        return Ridge(alpha=.001, fit_intercept=True)
    if kind == "boosting":
        return model_for(dict(kind="boosting", iterations=300))
    raise ValueError("Unknown phase-independent stacking decoder")


def fit_decoders(values, feature_rows, axes, cohort):
    rows, y, policies, table = training(axes, cohort)
    policies["matched"] = .75 * policies["clean"]
    policies["matched"][axes.iloc[rows].group.eq("imputed_training").to_numpy()] = .25 / 7
    table["weight_matched"] = policies["matched"]
    if len(feature_rows) != len(values) or len(set(feature_rows)) != len(feature_rows):
        raise ValueError("Stacking feature/sample axes mismatch")
    lookup = np.full(len(axes), -1, dtype=int)
    lookup[feature_rows] = np.arange(len(feature_rows))
    models = {}
    for mode, parameters in MODES.items():
        chosen = np.flatnonzero(policies[parameters["policy"]] > 0)
        selected = lookup[rows[chosen]]
        if np.any(selected < 0) or not axes.iloc[rows[chosen]].group.isin(
                ["clean_training", "imputed_training"]).all():
            raise ValueError("Stacking decoder attempted to use a nontraining feature")
        models[mode] = decoder_for(parameters["kind"]).fit(
            values[selected], y[chosen], sample_weight=policies[parameters["policy"]][chosen])
    return models, table


def prepare(dosage, trees, out):
    if (out / "experiment.json").exists():
        raise ValueError("Phase-independent stacking requires a fresh output")
    parent = json.loads((trees / "experiment.json").read_text())
    previous = json.loads((trees / "completion.json").read_text())
    checked = json.loads((trees / "verification/final_verification.json").read_text())
    if (previous["status"] != "complete" or checked["status"] != "complete_verified"
            or checked["completion_sha256"] != adapter.digest(trees / "completion.json")):
        raise ValueError("Reference-tree contrast must be independently verified")
    inputs = {**parent["input_sha256"], **previous["input_sha256"], **previous["output_sha256"]}
    inputs.update(provenance.fingerprint_files([
        trees / "completion.json", trees / "verification/final_verification.json"]))
    provenance.verify_hashes(inputs)
    implementations = [
        Path(__file__), Path(adapter.__file__), Path(provenance.__file__),
        PROJECT / "pedigree_benchmark/scripts/screen_dosage_models.py",
        PROJECT / "pedigree_benchmark/scripts/screen_fine_profiles.py",
        PROJECT / "pedigree_benchmark/scripts/screen_supervised_genotypes.py",
        PROJECT / "pedigree_benchmark/scripts/screen_reference_trees.py",
        PROJECT / "pedigree_benchmark/scripts/train_matched_global.py",
        PROJECT / "pedigree_benchmark/scripts/screen_v1_decoders.py",
        PROJECT / "pedigree_benchmark/scripts/evaluate_finalist.py",
        PROJECT / "pedigree_benchmark/scripts/screen_model_families.py",
        PROJECT / "run_phasefree_stack.sh", PROJECT / "tests/test_phasefree_stack.py",
        PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    adapter.write_json(out / "experiment.json", dict(
        status="configured", models=MODES, sources=SOURCES,
        classifier="400 ExtraTrees, leaf1; source/fine-balanced reference only, never clean/imputed training",
        features="7 reference-only classifier fractions +7 frozen shrink50 frequency fractions",
        ridge=dict(alpha=.001, fit_intercept=True),
        boosting=dict(iterations=300, max_leaf_nodes=15, min_samples_leaf=10,
                      learning_rate=.05, early_stopping=False, l2_regularization=0),
        decoder_policies=["clean", "matched"],
        matched_decoder_weights="37.5% pure +37.5% pedigree +25% imputed; no reference people",
        no_evaluation_in_fit=True,
        base_training_disjoint_from_all_decoder_training_identities=True,
        no_crossfit_required="Base classifier sees reference people only, disjoint from decoder people and founders",
        fixed_array_diagnostic_mean_tv_target=.10, production_promotion=False,
        array_target_is_development_diagnostic_not_certification=True,
        implementation_sha256=provenance.fingerprint_files(implementations), input_sha256=inputs))
    shutil.copyfile(dosage / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(dosage / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared3 stacking heads over a genuinely disjoint reference-only base classifier", flush=True)


def run(dosage, fine, out):
    settings = json.loads((out / "experiment.json").read_text())
    if settings["models"] != MODES or settings["sources"] != SOURCES:
        raise ValueError("Frozen stacking configuration changed")
    provenance.verify_hashes(settings["implementation_sha256"])
    provenance.verify_hashes(settings["input_sha256"])
    if (out / "completion.json").exists():
        raise ValueError("Completed stacking results must not be overwritten")
    x, phase_x, axes, phase_rows, markers = matrices(dosage)
    cohort = json.loads((out / "cohort.json").read_text())
    rows, y, policies, table = training(axes, cohort)
    ref_rows = np.flatnonzero(axes.group.eq("reference"))
    clean_rows = np.flatnonzero(axes.group.eq("clean_training"))
    extra_rows = np.flatnonzero(axes.group.eq("imputed_training"))
    with (dosage / "candidates/dosage_pca64_ridge_clean/model.pkl").open("rb") as stream:
        original = pickle.load(stream)
    with (fine / "candidates/fine_simplex_shrink50/model.pkl").open("rb") as stream:
        frequency = pickle.load(stream)
    if original["marker_keys"] != markers or frequency["marker_keys"] != markers:
        raise ValueError("Stacking component marker identities changed")
    start = time.monotonic()
    with threadpool_limits(limits=6):
        scores = original["pca"].transform(standardize(
            x[ref_rows], original["mean"], original["scale"])).astype(float)
        weights = policies["reference"][:len(ref_rows)]
        weights = weights / weights.sum()
        classifier = dict(
            kind="pca", sources=SOURCES, marker_keys=markers,
            mean=original["mean"], scale=original["scale"], pca=original["pca"],
            reference_only_preprocessing=True, reference_only_classifier=True,
            experiment_sha256=adapter.digest(out / "experiment.json"),
            head=model_for(dict(kind="extratrees", leaf=1)).fit(
                scores, y[:len(ref_rows)], sample_weight=weights))
        components = dict(classifier=classifier, frequency=frequency)
        feature_rows = np.r_[clean_rows, extra_rows, phase_rows]
        parts = [features(components, x[selected]) for selected in (clean_rows, extra_rows, phase_rows)]
        values = np.concatenate(parts)
        heldout = parts[-1]
        phase_features = features(components, phase_x)
        feature_phase_error = float(np.abs(heldout - phase_features).max())
        if feature_phase_error > 1e-10:
            raise ValueError("Stacking features depend on phase")
        frame = axes.iloc[feature_rows].reset_index(drop=True).copy()
        columns = [f"{kind}_{source}" for kind in ("classifier", "frequency") for source in SOURCES]
        frame[columns] = values
        frame.to_csv(out / "features.tsv", sep="\t", index=False)
        with (out / "feature_models.pkl").open("wb") as stream:
            pickle.dump(components, stream, protocol=pickle.HIGHEST_PROTOCOL)
        models, table = fit_decoders(values, feature_rows, axes, cohort)
        results, outputs = {}, [out / "features.tsv", out / "feature_models.pkl"]
        for mode, decoder in models.items():
            directory = out / "candidates" / mode
            if directory.exists():
                raise ValueError(f"Candidate exists; inspect before retry: {directory}")
            directory.mkdir(parents=True)
            bundle = dict(
                kind="phasefree_stack", sources=SOURCES, marker_keys=markers,
                feature_models=components, decoder=decoder, phase_independent=True,
                experiment_sha256=adapter.digest(out / "experiment.json"))
            path = directory / "model.pkl"
            with path.with_suffix(".pending").open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            path.with_suffix(".pending").replace(path)
            with path.open("rb") as stream:
                restored = pickle.load(stream)
            predictions = predict_decoder(decoder, heldout)
            restored_features = features(restored["feature_models"], x[phase_rows])
            reloaded = predict_decoder(restored["decoder"], restored_features)
            altered = predict_decoder(restored["decoder"], phase_features)
            reload_error = float(np.abs(predictions - reloaded).max())
            phase_error = float(np.abs(reloaded - altered).max())
            if reload_error > 1e-10 or phase_error > 1e-10:
                raise ValueError("Stacking model serialization or phase invariance failed")
            predictions_frame = axes.iloc[phase_rows].rename(columns={"group": "dataset"}).copy()
            for i, source in enumerate(SOURCES):
                predictions_frame["pred_" + source] = predictions[:, i]
            array = predictions_frame[predictions_frame.dataset.eq("array_imputed")].copy()
            if len(array) != 5:
                raise ValueError("Incomplete array-imputed evaluation cohort")
            array["TV_distance"] = np.abs(
                array[["pred_" + source for source in SOURCES]].to_numpy()
                - [.5, .5, 0, 0, 0, 0, 0]).sum(axis=1) / 2
            array.to_csv(directory / "array_imputed_predictions.tsv", sep="\t", index=False)
            predictions_frame[~predictions_frame.dataset.eq("array_imputed")].to_csv(
                directory / "predictions.tsv", sep="\t", index=False)
            score(out, mode)
            metrics = json.loads((directory / "metrics.json").read_text())
            array_tv = float(array.TV_distance.mean())
            results[mode] = dict(
                pure_mean_tv=metrics["pure"]["mean_tv"], pedigree_mean_tv=metrics["pedigree"]["mean_tv"],
                array_imputed_mean_tv=array_tv, real=metrics["real"],
                unchanged_gates_pass=metrics["all_gates_pass"],
                failed_gates=[key for key, value in metrics["gates"].items() if not value],
                array_diagnostic_target_pass=array_tv <= settings["fixed_array_diagnostic_mean_tv_target"],
                max_roundtrip_prediction_error=reload_error, max_phase_prediction_error=phase_error)
            outputs.extend(path for path in directory.iterdir() if path.is_file())
    table.to_csv(out / "training.tsv", sep="\t", index=False)
    pd.DataFrame([dict(candidate=mode, **{key: result[key] for key in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv",
        "unchanged_gates_pass", "array_diagnostic_target_pass")})
        for mode, result in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    adapter.write_json(out / "completion.json", dict(
        status="complete", n_candidates=3, n_markers=len(markers), n_chromosomes=22,
        n_reference_base_training=3447, n_clean_decoder_training=799, n_imputed_decoder_training=7,
        results=results, feature_phase_prediction_error=feature_phase_error,
        elapsed_seconds=time.monotonic() - start, production_promotion=False, independent_final_test=False,
        evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "training.tsv", out / "comparison.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "run"))
    parser.add_argument("--dosage", type=Path, default=Path("/dosage"))
    parser.add_argument("--fine", type=Path, default=Path("/fine"))
    parser.add_argument("--trees", type=Path, default=Path("/trees"))
    parser.add_argument("--outdir", type=Path, default=Path("/stack"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.dosage, args.trees, args.outdir)
    else:
        run(args.dosage, args.fine, args.outdir)


if __name__ == "__main__":
    main()
