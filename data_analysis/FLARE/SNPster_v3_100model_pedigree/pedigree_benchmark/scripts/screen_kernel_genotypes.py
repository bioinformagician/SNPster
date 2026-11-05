#!/usr/bin/env python3
"""Nonlinear global kernels on the frozen phase-independent genotype panel."""

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
from screen_supervised_genotypes import training
from screen_v1_decoders import simplex_rows


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
ALPHAS = (1e-5, 1e-4, 1e-3)
KINDS = ("poly2", "rbf05", "rbf2")
MODES = {f"kernel_{kind}_a{alpha:.0e}": dict(kind=kind, alpha=alpha)
         for kind in KINDS for alpha in ALPHAS}


def kernel_values(dot, left_norm, right_norm, n_features, kind):
    if (dot.ndim != 2 or left_norm.shape != (dot.shape[0],)
            or right_norm.shape != (dot.shape[1],) or n_features <= 0
            or not np.isfinite(dot).all() or not np.isfinite(left_norm).all()
            or not np.isfinite(right_norm).all() or np.any(left_norm < 0)
            or np.any(right_norm < 0)):
        raise ValueError("Invalid global kernel geometry")
    if kind == "poly2":
        return (1 + dot / n_features) ** 2
    if kind not in ("rbf05", "rbf2"):
        raise ValueError(f"Unknown global genotype kernel: {kind}")
    distance = (left_norm[:, None] + right_norm[None, :] - 2 * dot) / n_features
    if np.any(distance < -1e-10):
        raise ValueError("Inconsistent kernel norms/dot products")
    gamma = .5 if kind == "rbf05" else 2.
    return np.exp(-gamma * np.maximum(distance, 0))


def kernel_heads(kernel, y, weights, alphas):
    if (kernel.ndim != 2 or kernel.shape != (len(y), len(y)) or y.ndim != 2
            or weights.shape != (len(y),) or not len(y)
            or not np.isfinite(kernel).all() or not np.isfinite(y).all()
            or not np.isfinite(weights).all() or np.any(weights <= 0)
            or not np.isclose(weights.sum(), 1, atol=1e-12)
            or not np.allclose(kernel, kernel.T, atol=1e-12, rtol=0)
            or not alphas or any(not np.isfinite(alpha) or alpha <= 0 for alpha in alphas)):
        raise ValueError("Invalid weighted global-kernel training input")
    root = np.sqrt(weights)
    weighted = kernel * root[:, None] * root[None, :]
    rhs, models = np.column_stack([root[:, None] * y, root]), {}
    for alpha in alphas:
        system = weighted.copy()
        system.flat[::len(y) + 1] += alpha
        solved = cho_solve(cho_factor(system, lower=True), rhs)
        denominator = root @ solved[:, -1]
        if not np.isfinite(denominator) or denominator <= 0:
            raise ValueError("Invalid unpenalized kernel-intercept solve")
        intercept = (root @ solved[:, :-1]) / denominator
        dual = root[:, None] * (solved[:, :-1] - solved[:, -1, None] * intercept)
        residual = kernel @ dual + intercept - y
        error = float(np.max(np.abs(residual + alpha * dual / weights[:, None])))
        intercept_error = float(np.max(np.abs(weights @ residual)))
        if error > 1e-8 or intercept_error > 1e-10:
            raise ValueError(f"Global kernel solve failed normal equations: {error}, {intercept_error}")
        models[alpha] = dict(dual=dual, intercept=intercept, alpha=alpha,
                             max_normal_equation_error=error,
                             max_intercept_equation_error=intercept_error)
    return models


def predict(bundle, x):
    validate_dosages(x)
    validate_dosages(bundle["anchor_dosages"])
    z = standardize(x, bundle["mean"], bundle["scale"]).astype(float)
    anchors = standardize(bundle["anchor_dosages"], bundle["mean"], bundle["scale"]).astype(float)
    kernel = kernel_values(z @ anchors.T, (z * z).sum(axis=1), (anchors * anchors).sum(axis=1),
                           z.shape[1], bundle["kernel_kind"])
    return simplex_rows(kernel @ bundle["dual"] + bundle["intercept"])


def prepare(linear, out):
    if (out / "experiment.json").exists():
        raise ValueError("Global-kernel preparation requires a fresh output")
    parent = json.loads((linear / "experiment.json").read_text())
    previous = json.loads((linear / "completion.json").read_text())
    if previous["status"] != "complete" or any(
            result["unchanged_gates_pass"] for result in previous["results"].values()):
        raise ValueError("Expected the completed negative linear-model contrast before this batch")
    inputs = {**parent["input_sha256"], **previous["input_sha256"], **previous["output_sha256"]}
    inputs.update(provenance.fingerprint_files([linear / "completion.json"]))
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
        PROJECT / "run_kernel_genotype_screen.sh", PROJECT / "tests/test_kernel_genotypes.py",
        PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    adapter.write_json(out / "experiment.json", dict(
        status="configured", models=MODES, sources=SOURCES,
        kernels=dict(poly2="(1 + z@z'/26066)**2",
                     rbf05="exp(-0.5*||z-z'||**2/26066)",
                     rbf2="exp(-2*||z-z'||**2/26066)"),
        alphas=list(ALPHAS), training_policy="reference",
        reference_weights="50% source/fine-balanced reference +25% pure +25% family-balanced pedigree",
        no_imputed_training="Prior linear contrast showed no improvement from the seven matched cases",
        reference_only_scaling=True, n_features=26066, no_evaluation_in_fit=True,
        fixed_array_diagnostic_mean_tv_target=.10, production_promotion=False,
        array_target_is_development_diagnostic_not_certification=True,
        implementation_sha256=provenance.fingerprint_files(paths), input_sha256=inputs))
    shutil.copyfile(linear / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(linear / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared9 global kernel heads; no chromosome fits, haplotypes or new imputation", flush=True)


def run(dosage, linear, out, smoke=False):
    settings = json.loads((out / "experiment.json").read_text())
    if settings["models"] != MODES or settings["sources"] != SOURCES:
        raise ValueError("Frozen kernel batch configuration changed")
    provenance.verify_hashes(settings["implementation_sha256"])
    provenance.verify_hashes(settings["input_sha256"])
    if not smoke and (out / "completion.json").exists():
        raise ValueError("Completed global-kernel results must not be overwritten")
    x, phase_x, axes, phase_rows, markers = matrices(dosage)
    cohort = json.loads((out / "cohort.json").read_text())
    rows, truth, policies, table = training(axes, cohort)
    chosen = np.flatnonzero(policies["reference"] > 0)
    weights, y = policies["reference"][chosen], truth[chosen]
    anchors = x[rows[chosen]]
    reference = x[np.flatnonzero(axes.group.eq("reference"))]
    mean = reference.mean(axis=0, dtype=np.float64)
    p = mean / 2
    scale = np.sqrt(2 * p * (1 - p))
    z = standardize(anchors, mean, scale).astype(float)
    if smoke:
        selected = np.r_[np.arange(0, 3447, 30), np.arange(3447, len(z))]
        z, y, weights = z[selected, :256], y[selected], weights[selected]
        weights = weights / weights.sum()
    outputs, results, start = [], {}, time.monotonic()
    with threadpool_limits(limits=6):
        dot, norms = z @ z.T, (z * z).sum(axis=1)
        if not smoke:
            eval_z = standardize(x[phase_rows], mean, scale).astype(float)
            eval_dot, eval_norms = eval_z @ z.T, (eval_z * eval_z).sum(axis=1)
        for kind in KINDS:
            kernel = kernel_values(dot, norms, norms, z.shape[1], kind)
            models = kernel_heads(kernel, y, weights, ALPHAS)
            if smoke:
                print(f"{kind}: training-only compatibility solve passed for3 alphas", flush=True)
                continue
            cross = kernel_values(eval_dot, eval_norms, norms, z.shape[1], kind)
            for mode, parameters in MODES.items():
                if parameters["kind"] != kind:
                    continue
                directory = out / "candidates" / mode
                if directory.exists():
                    raise ValueError(f"Candidate exists; inspect before retry: {directory}")
                directory.mkdir(parents=True)
                bundle = dict(
                    kind="global_genotype_kernel", kernel_kind=kind, marker_keys=markers,
                    sources=SOURCES, mean=mean, scale=scale, anchor_dosages=anchors,
                    anchor_ids=table.iloc[chosen][["group", "sample_id"]].to_dict("records"),
                    reference_only_scaling=True, phase_independent=True, training_policy="reference",
                    experiment_sha256=adapter.digest(out / "experiment.json"),
                    **models[parameters["alpha"]])
                path = directory / "model.pkl"
                with path.with_suffix(".pending").open("wb") as stream:
                    pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
                path.with_suffix(".pending").replace(path)
                with path.open("rb") as stream:
                    restored = pickle.load(stream)
                predictions = simplex_rows(cross @ bundle["dual"] + bundle["intercept"])
                reloaded, altered = predict(restored, x[phase_rows]), predict(restored, phase_x)
                reload_error = float(np.max(np.abs(predictions - reloaded)))
                phase_error = float(np.max(np.abs(reloaded - altered)))
                if reload_error > 1e-10 or phase_error > 1e-10:
                    raise ValueError("Kernel serialization or phase invariance failed")
                frame = axes.iloc[phase_rows].rename(columns={"group": "dataset"}).copy()
                for i, source in enumerate(SOURCES):
                    frame["pred_" + source] = predictions[:, i]
                array = frame[frame.dataset.eq("array_imputed")].copy()
                if len(array) != 5:
                    raise ValueError("Incomplete array-imputed cohort")
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
    table.iloc[chosen].to_csv(out / "training.tsv", sep="\t", index=False)
    pd.DataFrame([dict(candidate=mode, **{key: result[key] for key in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv",
        "unchanged_gates_pass", "array_diagnostic_target_pass")})
        for mode, result in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    adapter.write_json(out / "completion.json", dict(
        status="complete", n_candidates=9, n_markers=len(markers), n_chromosomes=22,
        n_reference=3447, n_clean_training=799, n_imputed_training=0,
        results=results, retained_linear_results=json.loads(
            (linear / "completion.json").read_text())["results"],
        elapsed_seconds=time.monotonic() - start, production_promotion=False,
        independent_final_test=False, evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "training.tsv", out / "comparison.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "smoke", "run"))
    parser.add_argument("--linear", type=Path, default=Path("/supervised"))
    parser.add_argument("--dosage", type=Path, default=Path("/dosage"))
    parser.add_argument("--outdir", type=Path, default=Path("/kernel"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.linear, args.outdir)
    else:
        run(args.dosage, args.linear, args.outdir, smoke=args.stage == "smoke")


if __name__ == "__main__":
    main()
