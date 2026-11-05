#!/usr/bin/env python3
"""Full-marker, phase-independent supervised global ancestry comparison."""

import argparse
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from screen_dosage_models import standardize
from screen_fine_profiles import matrices, validate_dosages
from screen_model_families import score
from screen_v1_decoders import simplex_rows, training_weights
from train_matched_global import training_table


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
ALPHAS = (.1, 1., 10.)
POLICIES = ("clean", "reference", "matched")
MODES = {f"genotype_ridge_{policy}_a{str(alpha).replace('.', 'p')}":
         dict(policy=policy, alpha=alpha) for policy in POLICIES for alpha in ALPHAS}


def training(axes, cohort):
    pure, ped = training_table(cohort)
    metadata = pd.concat([pure, ped], ignore_index=True)
    ref = pd.read_csv(cohort["reference"]["labels_tsv"], sep="\t")
    ref_rows = np.flatnonzero(axes.group.eq("reference"))
    clean_rows = np.flatnonzero(axes.group.eq("clean_training"))
    extra_rows = np.flatnonzero(axes.group.eq("imputed_training"))
    clean_ids = axes.iloc[clean_rows].sample_id
    extra_ids = axes.iloc[extra_rows].sample_id
    if (len(clean_ids) != 799 or metadata.sample_id.duplicated().any()
            or len(set(clean_ids)) != 799 or set(clean_ids) != set(metadata.sample_id)
            or len(extra_ids) != 7 or len(set(extra_ids)) != 7
            or not set(extra_ids) <= set(ped.sample_id)):
        raise ValueError("Supervised training identities are incomplete or include a nontraining sample")
    if (axes.iloc[ref_rows].sample_id.tolist() != ref.sample_id.tolist()
            or ref.sample_id.duplicated().any() or set(ref.source_group) != set(SOURCES)
            or ref.fine.isna().any() or ref.fine.astype(str).str.strip().eq("").any()
            or ref.groupby("fine").source_group.nunique().ne(1).any()
            or set(ref.sample_id) & set(metadata.sample_id)):
        raise ValueError("Reference identities/fine labels are invalid or overlap clean training")
    train_ids = set(ref.sample_id) | set(clean_ids)
    if train_ids & set(axes.loc[axes.group.isin(["pure", "pedigree", "real", "array_imputed"]), "sample_id"]):
        raise ValueError("Held-out identity entered supervised training")
    truth = ["truth_" + source for source in SOURCES]
    ref_y = np.column_stack([ref.source_group.eq(source).to_numpy(float) for source in SOURCES])
    y = np.concatenate([ref_y, metadata.set_index("sample_id").loc[clean_ids, truth].to_numpy(float),
                        ped.set_index("sample_id").loc[extra_ids, truth].to_numpy(float)])
    if (not np.isfinite(y).all() or np.any(y < 0)
            or not np.allclose(y.sum(axis=1), 1, atol=1e-12, rtol=0)):
        raise ValueError("Invalid supervised global truth")
    clean_weights = pd.Series(training_weights(pure, ped), index=metadata.sample_id).loc[clean_ids].to_numpy()
    fine_counts = ref.groupby(["source_group", "fine"]).sample_id.transform("size").to_numpy(float)
    source_populations = ref.groupby("source_group").fine.transform("nunique").to_numpy(float)
    ref_weights = 1 / (len(SOURCES) * fine_counts * source_populations)
    if not np.isclose(ref_weights.sum(), 1, atol=1e-12):
        raise ValueError("Invalid source/fine-balanced reference weights")
    n_ref = len(ref)
    clean = np.r_[np.zeros(n_ref), clean_weights, np.zeros(7)]
    reference = np.r_[.5 * ref_weights, .5 * clean_weights, np.zeros(7)]
    matched = .75 * reference
    matched[-7:] = .25 / 7
    rows = np.r_[ref_rows, clean_rows, extra_rows]
    table = axes.iloc[rows].reset_index(drop=True).copy()
    for i, source in enumerate(SOURCES):
        table["truth_" + source] = y[:, i]
    for policy, weights in (("clean", clean), ("reference", reference), ("matched", matched)):
        table["weight_" + policy] = weights
    return rows, y, {policy: table["weight_" + policy].to_numpy() for policy in POLICIES}, table


def weighted_ridges(x, y, weights, alphas, gram=None):
    x, y, weights = np.asarray(x, float), np.asarray(y, float), np.asarray(weights, float)
    if (x.ndim != 2 or y.ndim != 2 or len(x) != len(y) or not len(x)
            or weights.shape != (len(x),) or not np.isfinite(x).all()
            or not np.isfinite(y).all() or not np.isfinite(weights).all()
            or np.any(weights <= 0) or not np.isclose(weights.sum(), 1, atol=1e-12)
            or not alphas or any(not np.isfinite(alpha) or alpha <= 0 for alpha in alphas)):
        raise ValueError("Invalid weighted full-marker Ridge inputs")
    gram = x @ x.T if gram is None else np.asarray(gram, float)
    if gram.shape != (len(x), len(x)) or not np.isfinite(gram).all():
        raise ValueError("Invalid full-marker Gram matrix")
    row_mean = gram @ weights
    centered = gram - row_mean[:, None] - row_mean[None, :] + weights @ row_mean
    root = np.sqrt(weights)
    centered *= root[:, None] * root[None, :]
    x_mean, y_mean = weights @ x, weights @ y
    rhs = root[:, None] * (y - y_mean)
    models = {}
    for alpha in alphas:
        kernel = centered.copy()
        kernel.flat[::len(x) + 1] += alpha
        dual = root[:, None] * cho_solve(cho_factor(kernel, lower=True), rhs)
        coefficients = x.T @ dual
        intercept = y_mean - x_mean @ coefficients
        residual = x @ coefficients + intercept - y
        gradient = x.T @ (weights[:, None] * residual) + alpha * coefficients
        error = float(np.max(np.abs(gradient)))
        intercept_error = float(np.max(np.abs(weights @ residual)))
        if error > 1e-8 or intercept_error > 1e-10:
            raise ValueError(f"Full-marker Ridge solve failed normal equations: {error}, {intercept_error}")
        models[alpha] = dict(coefficients=coefficients, intercept=intercept, alpha=alpha,
                             max_normal_equation_error=error,
                             max_intercept_equation_error=intercept_error)
    return models


def predict(bundle, x):
    validate_dosages(x)
    transformed = standardize(x, bundle["mean"], bundle["scale"]).astype(float)
    return simplex_rows(transformed @ bundle["coefficients"] + bundle["intercept"])


def prepare(dosage, profiles, out):
    if (out / "experiment.json").exists():
        raise ValueError("Supervised genotype preparation requires a new output")
    old = json.loads((profiles / "completion.json").read_text())
    checked = json.loads((profiles / "verification/final_verification.json").read_text())
    if (old["status"] != "complete" or checked["status"] != "complete_verified"
            or checked["completion_sha256"] != adapter.digest(profiles / "completion.json")):
        raise ValueError("Fine-profile comparison must already be independently verified")
    parent = json.loads((profiles / "experiment.json").read_text())
    cohort = json.loads((dosage / "cohort.json").read_text())
    adapter.validate_splits(cohort)
    inputs = {**parent["input_sha256"], **old["input_sha256"], **old["output_sha256"]}
    original = json.loads((dosage / "experiment.json").read_text())
    metadata = [Path(cohort["pure_train"]["labels_tsv"])]
    metadata += [Path(batch["samples_tsv"]) for batch in cohort["pedigree_train"]["batches"]]
    for path in metadata:
        inputs[str(path)] = original["input_sha256"][str(path)]
    inputs.update(provenance.fingerprint_files([
        profiles / "completion.json", profiles / "verification/final_verification.json"]))
    provenance.verify_hashes(inputs)
    implementation = [
        Path(__file__), Path(adapter.__file__), Path(provenance.__file__),
        PROJECT / "pedigree_benchmark/scripts/screen_dosage_models.py",
        PROJECT / "pedigree_benchmark/scripts/screen_fine_profiles.py",
        PROJECT / "pedigree_benchmark/scripts/train_matched_global.py",
        PROJECT / "pedigree_benchmark/scripts/screen_v1_decoders.py",
        PROJECT / "pedigree_benchmark/scripts/evaluate_finalist.py",
        PROJECT / "pedigree_benchmark/scripts/screen_model_families.py",
        PROJECT / "run_supervised_genotype_screen.sh", PROJECT / "tests/test_supervised_genotypes.py",
        PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    settings = dict(
        status="configured", models=MODES, alphas=list(ALPHAS), policies=list(POLICIES),
        sources=SOURCES, n_features=26066, reference_only_scaling=True,
        features="Every frozen thinned hard-GT marker, reference binomial scaling, float32 then float64",
        regression="Weighted squared global truth error + alpha*sum(coefficients**2); unpenalized intercept",
        clean_weights="Half source-balanced pure, half family-balanced pedigree",
        reference_weights="Half source/fine-balanced reference, half clean policy",
        matched_weights="75% reference policy +25% equally weighted7 imputed training descendants",
        fixed_array_diagnostic_mean_tv_target=.10,
        array_target_is_development_diagnostic_not_certification=True,
        no_evaluation_in_fit=True, no_marker_or_pca_selection=True, no_new_imputation=True,
        baselines_immutable=True, production_promotion=False,
        implementation_sha256=provenance.fingerprint_files(implementation), input_sha256=inputs)
    adapter.write_json(out / "experiment.json", settings)
    shutil.copyfile(dosage / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(dosage / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared9 full-marker Ridge models:3 training policies x3 fixed regularizers", flush=True)


def configuration(out):
    settings = json.loads((out / "experiment.json").read_text())
    if settings["models"] != MODES or settings["sources"] != SOURCES:
        raise ValueError("Frozen supervised genotype configuration changed")
    provenance.verify_hashes(settings["implementation_sha256"])
    provenance.verify_hashes(settings["input_sha256"])
    return settings


def run(dosage, profiles, out, smoke=False):
    settings = configuration(out)
    if not smoke and (out / "completion.json").exists():
        raise ValueError("Completed supervised models must not be overwritten")
    x, phase_x, axes, phase_rows, markers = matrices(dosage)
    cohort = json.loads((out / "cohort.json").read_text())
    rows, y, weights, table = training(axes, cohort)
    ref = x[np.flatnonzero(axes.group.eq("reference"))]
    mean = ref.mean(axis=0, dtype=np.float64)
    p = mean / 2
    scale = np.sqrt(2 * p * (1 - p))
    train_x = standardize(x[rows], mean, scale).astype(float)
    if smoke:
        selected = np.r_[np.arange(0, 3447, 30), np.arange(3447, len(rows))]
        train_x, y = train_x[selected, :256], y[selected]
        weights = {policy: value[selected] / value[selected].sum() for policy, value in weights.items()}
    results, outputs, start = {}, [], time.monotonic()
    with threadpool_limits(limits=6):
        gram = train_x @ train_x.T
        for policy in POLICIES:
            chosen = np.flatnonzero(weights[policy] > 0)
            fitted = weighted_ridges(train_x[chosen], y[chosen], weights[policy][chosen],
                                     ALPHAS, gram[np.ix_(chosen, chosen)])
            if smoke:
                print(f"{policy}: training-only256-marker solve passed for all3 alphas", flush=True)
                continue
            for mode, parameters in MODES.items():
                if parameters["policy"] != policy:
                    continue
                directory = out / "candidates" / mode
                if directory.exists():
                    raise ValueError(f"Candidate directory exists; inspect before retry: {directory}")
                directory.mkdir(parents=True)
                bundle = dict(
                    kind="full_marker_ridge", sources=SOURCES, marker_keys=markers,
                    mean=mean, scale=scale, training_policy=policy,
                    phase_independent=True, reference_only_scaling=True,
                    experiment_sha256=adapter.digest(out / "experiment.json"),
                    **fitted[parameters["alpha"]])
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
                    raise ValueError("Full-marker model persistence or phase invariance failed")
                frame = axes.iloc[phase_rows].rename(columns={"group": "dataset"}).copy()
                for i, source in enumerate(SOURCES):
                    frame["pred_" + source] = predictions[:, i]
                array = frame[frame.dataset.eq("array_imputed")].copy()
                if len(array) != 5:
                    raise ValueError("Missing array-imputed evaluation case")
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
                    max_phase_prediction_error=phase_error, max_roundtrip_prediction_error=reload_error,
                    max_normal_equation_error=bundle["max_normal_equation_error"])
                outputs.extend(path for path in directory.iterdir() if path.is_file())
    if smoke:
        adapter.write_json(out / "smoke.json", dict(
            status="complete", training_only=True, n_markers=256,
            experiment_sha256=adapter.digest(out / "experiment.json")))
        return
    table.to_csv(out / "training.tsv", sep="\t", index=False)
    pd.DataFrame([dict(candidate=mode, **{key: result[key] for key in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv",
        "unchanged_gates_pass", "array_diagnostic_target_pass")})
        for mode, result in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    adapter.write_json(out / "completion.json", dict(
        status="complete", n_candidates=9, n_markers=len(markers), n_chromosomes=22,
        n_reference=3447, n_clean_training=799, n_imputed_training=7,
        results=results, retained_fine_profile_results=json.loads(
            (profiles / "completion.json").read_text())["results"],
        retained_dosage_results=json.loads((dosage / "completion.json").read_text())["results"],
        evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        elapsed_seconds=time.monotonic() - start, production_promotion=False,
        independent_final_test=False,
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "training.tsv", out / "comparison.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "smoke", "run"))
    parser.add_argument("--dosage", type=Path, default=Path("/dosage"))
    parser.add_argument("--profiles", type=Path, default=Path("/fine"))
    parser.add_argument("--outdir", type=Path, default=Path("/supervised"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.dosage, args.profiles, args.outdir)
    else:
        run(args.dosage, args.profiles, args.outdir, smoke=args.stage == "smoke")


if __name__ == "__main__":
    main()
