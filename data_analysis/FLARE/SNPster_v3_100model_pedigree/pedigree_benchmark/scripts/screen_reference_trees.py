#!/usr/bin/env python3
"""Reference-augmented nonlinear heads on the frozen reference-only PCA."""

import argparse
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from screen_dosage_models import predict
from screen_fine_profiles import matrices
from screen_model_families import score
from screen_supervised_genotypes import training


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
MODES = {
    "pca64_reference_et_leaf1": dict(kind="extratrees", leaf=1),
    "pca64_reference_et_leaf5": dict(kind="extratrees", leaf=5),
    "pca64_reference_hgb150": dict(kind="boosting", iterations=150),
    "pca64_reference_hgb400": dict(kind="boosting", iterations=400),
}


def model_for(parameters):
    if parameters["kind"] == "extratrees":
        return ExtraTreesRegressor(
            n_estimators=400, min_samples_leaf=parameters["leaf"], max_features=1.,
            random_state=adapter.SEED, n_jobs=1)
    if parameters["kind"] != "boosting":
        raise ValueError("Unknown reference-tree model")
    return MultiOutputRegressor(HistGradientBoostingRegressor(
        max_iter=parameters["iterations"], learning_rate=.05, max_leaf_nodes=15,
        min_samples_leaf=10, l2_regularization=0, early_stopping=False,
        random_state=adapter.SEED), n_jobs=1)


def fit_heads(scores, axes, cohort):
    rows, y, policies, table = training(axes, cohort)
    chosen = np.flatnonzero(policies["reference"] > 0)
    x, truth, weights = scores[rows[chosen]], y[chosen], policies["reference"][chosen]
    models = {}
    for mode, parameters in MODES.items():
        models[mode] = model_for(parameters).fit(x, truth, sample_weight=weights)
        print(f"Trained {mode} on{len(x)} reference/clean rows, no held-out fitting", flush=True)
    return models, table.iloc[chosen]


def prepare(dosage, kernel, out):
    if (out / "experiment.json").exists():
        raise ValueError("Reference-tree comparison requires a fresh output")
    parent = json.loads((kernel / "experiment.json").read_text())
    previous = json.loads((kernel / "completion.json").read_text())
    if previous["status"] != "complete" or any(
            value["unchanged_gates_pass"] for value in previous["results"].values()):
        raise ValueError("Expected completed negative kernel contrast before this final batch")
    inputs = {**parent["input_sha256"], **previous["input_sha256"], **previous["output_sha256"]}
    inputs.update(provenance.fingerprint_files([kernel / "completion.json"]))
    provenance.verify_hashes(inputs)
    paths = [
        Path(__file__), Path(adapter.__file__), Path(provenance.__file__),
        PROJECT / "pedigree_benchmark/scripts/screen_dosage_models.py",
        PROJECT / "pedigree_benchmark/scripts/screen_fine_profiles.py",
        PROJECT / "pedigree_benchmark/scripts/screen_supervised_genotypes.py",
        PROJECT / "pedigree_benchmark/scripts/train_matched_global.py",
        PROJECT / "pedigree_benchmark/scripts/screen_v1_decoders.py",
        PROJECT / "pedigree_benchmark/scripts/evaluate_finalist.py",
        PROJECT / "pedigree_benchmark/scripts/screen_model_families.py",
        PROJECT / "run_reference_tree_screen.sh", PROJECT / "tests/test_reference_trees.py",
        PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    adapter.write_json(out / "experiment.json", dict(
        status="configured", models=MODES, sources=SOURCES,
        pca_bundle=str(dosage / "candidates/dosage_pca64_ridge_clean/model.pkl"),
        pca="Exact frozen reference-only64-PC scaling/whitening; no new PCA fitting",
        training_policy="reference",
        training_weights="50% source/fine-balanced reference +25% pure +25% family-balanced pedigree",
        trees=dict(n_estimators=400, max_features=1., min_samples_leaf=[1, 5]),
        boosting=dict(max_iter=[150, 400], max_leaf_nodes=15, min_samples_leaf=10,
                      learning_rate=.05, l2_regularization=0, early_stopping=False),
        no_evaluation_in_fit=True, no_imputed_training=True,
        fixed_array_diagnostic_mean_tv_target=.10,
        array_target_is_development_diagnostic_not_certification=True, production_promotion=False,
        implementation_sha256=provenance.fingerprint_files(paths), input_sha256=inputs))
    shutil.copyfile(dosage / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(dosage / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared4 fixed reference-augmented tree/boosting models", flush=True)


def run(dosage, out):
    settings = json.loads((out / "experiment.json").read_text())
    if settings["models"] != MODES or settings["sources"] != SOURCES:
        raise ValueError("Frozen reference-tree configuration changed")
    provenance.verify_hashes(settings["implementation_sha256"])
    provenance.verify_hashes(settings["input_sha256"])
    if (out / "completion.json").exists():
        raise ValueError("Completed reference-tree results must not be overwritten")
    x, phase_x, axes, phase_rows, markers = matrices(dosage)
    with Path(settings["pca_bundle"]).open("rb") as stream:
        frozen = pickle.load(stream)
    if frozen["marker_keys"] != markers or not frozen["reference_only_preprocessing"]:
        raise ValueError("Frozen reference-only PCA or marker axis changed")
    from screen_dosage_models import standardize
    with threadpool_limits(limits=6):
        scores = frozen["pca"].transform(standardize(x, frozen["mean"], frozen["scale"])).astype(float)
        start = time.monotonic()
        models, table = fit_heads(scores, axes, json.loads((out / "cohort.json").read_text()))
        results, outputs = {}, []
        for mode, head in models.items():
            directory = out / "candidates" / mode
            if directory.exists():
                raise ValueError(f"Candidate directory exists; inspect before retry: {directory}")
            directory.mkdir(parents=True)
            bundle = dict(
                kind="pca", sources=SOURCES, marker_keys=markers, head=head,
                mean=frozen["mean"], scale=frozen["scale"], pca=frozen["pca"],
                reference_only_preprocessing=True, phase_independent=True,
                training_policy="reference", experiment_sha256=adapter.digest(out / "experiment.json"))
            path = directory / "model.pkl"
            with path.with_suffix(".pending").open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            path.with_suffix(".pending").replace(path)
            with path.open("rb") as stream:
                restored = pickle.load(stream)
            predictions = predict(bundle, x[phase_rows])
            reloaded, altered = predict(restored, x[phase_rows]), predict(restored, phase_x)
            reload_error = float(np.max(np.abs(predictions - reloaded)))
            phase_error = float(np.max(np.abs(reloaded - altered)))
            if reload_error > 1e-10 or phase_error > 1e-10:
                raise ValueError("Reference-tree serialization or phase-invariance failure")
            frame = axes.iloc[phase_rows].rename(columns={"group": "dataset"}).copy()
            for i, source in enumerate(SOURCES):
                frame["pred_" + source] = predictions[:, i]
            array = frame[frame.dataset.eq("array_imputed")].copy()
            if len(array) != 5:
                raise ValueError("Incomplete array-imputed evaluation cohort")
            array["TV_distance"] = np.abs(
                array[["pred_" + source for source in SOURCES]].to_numpy()
                - [.5, .5, 0, 0, 0, 0, 0]).sum(axis=1) / 2
            array.to_csv(directory / "array_imputed_predictions.tsv", sep="\t", index=False)
            frame[~frame.dataset.eq("array_imputed")].to_csv(
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
                max_phase_prediction_error=phase_error, max_roundtrip_prediction_error=reload_error)
            outputs.extend(path for path in directory.iterdir() if path.is_file())
    table.to_csv(out / "training.tsv", sep="\t", index=False)
    pd.DataFrame([dict(candidate=mode, **{key: result[key] for key in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv",
        "unchanged_gates_pass", "array_diagnostic_target_pass")})
        for mode, result in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    adapter.write_json(out / "completion.json", dict(
        status="complete", n_candidates=4, n_markers=len(markers), n_chromosomes=22,
        n_reference=3447, n_clean_training=799, n_imputed_training=0, results=results,
        elapsed_seconds=time.monotonic() - start, production_promotion=False, independent_final_test=False,
        evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "training.tsv", out / "comparison.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "run"))
    parser.add_argument("--dosage", type=Path, default=Path("/dosage"))
    parser.add_argument("--kernel", type=Path, default=Path("/kernel"))
    parser.add_argument("--outdir", type=Path, default=Path("/trees"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.dosage, args.kernel, args.outdir)
    else:
        run(args.dosage, args.outdir)


if __name__ == "__main__":
    main()
