#!/usr/bin/env python3
"""Compare reference-only fine-population mixtures on frozen genotype matrices."""

import argparse
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from screen_model_families import score


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
MODES = {"fine_simplex_raw": 0, "fine_simplex_shrink10": 10, "fine_simplex_shrink50": 50}


def validate_dosages(x):
    if (x.ndim != 2 or not np.all(np.isfinite(x))
            or not np.isin(x, [0, 1, 2]).all()):
        raise ValueError("Complete hard-GT allele counts 0/1/2 are required")


def build_profiles(x, axes, labels, strength):
    if not np.isfinite(strength) or strength < 0:
        raise ValueError("Shrinkage strength must be nonnegative")
    if x.ndim != 2 or len(x) != len(axes):
        raise ValueError("Genotype sample axis mismatch")
    reference_rows = np.flatnonzero(axes.group.eq("reference").to_numpy())
    if (labels.sample_id.duplicated().any() or labels.fine.isna().any()
            or labels.fine.astype(str).str.strip().eq("").any()
            or set(labels.source_group) != set(SOURCES)
            or axes.iloc[reference_rows].sample_id.tolist() != labels.sample_id.tolist()
            or labels.groupby("fine").source_group.nunique().ne(1).any()):
        raise ValueError("Reference identities, fine labels or source assignments are inconsistent")
    reference = x[reference_rows]
    validate_dosages(reference)
    broad = np.column_stack([
        reference[labels.source_group.eq(source)].mean(axis=0, dtype=np.float64) / 2
        for source in SOURCES])
    informative = np.ptp(broad, axis=1) >= .03
    if informative.sum() < 100:
        raise ValueError("Insufficient informative markers on the frozen broad-source mask")
    profiles, frequencies = [], []
    for source_index, source in enumerate(SOURCES):
        for fine in sorted(labels.loc[labels.source_group.eq(source), "fine"].unique()):
            rows = labels.source_group.eq(source) & labels.fine.eq(fine)
            count = int(rows.sum())
            weight = count / (count + strength)
            fine_af = reference[rows].mean(axis=0, dtype=np.float64) / 2
            frequencies.append(weight * fine_af + (1 - weight) * broad[:, source_index])
            profiles.append(dict(profile_id=f"{source}:{fine}", source_group=source,
                                 fine=fine, n_reference=count, fine_weight=weight))
    source_indices = np.array([SOURCES.index(row["source_group"]) for row in profiles])
    source_counts = np.bincount(source_indices, minlength=len(SOURCES))
    return dict(
        kind="fine_population_simplex", sources=SOURCES, profiles=profiles,
        frequencies=np.column_stack(frequencies), informative=informative,
        source_indices=source_indices, shrinkage_pseudoindividuals=strength,
        initial_weights=1 / (len(SOURCES) * source_counts[source_indices]),
        phase_independent=True, reference_only_preprocessing=True)


def predict(bundle, x):
    validate_dosages(x)
    frequencies = bundle["frequencies"]
    informative = bundle["informative"]
    if (frequencies.ndim != 2 or x.shape[1] != len(frequencies)
            or informative.shape != (len(frequencies),)
            or informative.dtype != bool or informative.sum() < 100
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies < 0) or np.any(frequencies > 1)):
        raise ValueError("Invalid fine-population frequency model or marker axis")
    f = frequencies[informative]
    # With sum(q)=1, subtracting a shared marker center leaves squared error unchanged.
    center = f.mean(axis=1)
    f = f - center[:, None]
    gram = f.T @ f / len(f)
    target = (x[:, informative] / 2 - center) @ f / len(f)
    weights, diagnostics = [], []
    for vector in target:
        fit = minimize(
            lambda q: q @ gram @ q - 2 * vector @ q, bundle["initial_weights"].copy(),
            jac=lambda q: 2 * (gram @ q - vector),
            constraints={"type": "eq", "fun": lambda q: q.sum() - 1,
                         "jac": lambda q: np.ones_like(q)},
            bounds=[(0, 1)] * len(vector), method="SLSQP",
            options={"ftol": 1e-12, "maxiter": 1000})
        if not fit.success:
            raise ValueError(f"Fine-profile optimization failed: {fit.message}")
        gradient = 2 * (gram @ fit.x - vector)
        gap = float(fit.x @ gradient - gradient.min())
        if (not np.all(np.isfinite(fit.x)) or np.any(fit.x < -1e-10)
                or not np.isclose(fit.x.sum(), 1, atol=1e-8, rtol=0) or gap > 1e-6):
            raise ValueError(f"Invalid or insufficiently optimized fine-profile weights; gap={gap}")
        weights.append(fit.x)
        diagnostics.append(dict(iterations=int(fit.nit), simplex_dual_gap=gap))
    weights = np.asarray(weights)
    membership = np.eye(len(SOURCES))[bundle["source_indices"]]
    return weights @ membership, weights, pd.DataFrame(diagnostics)


def matrices(dosage):
    arrays, altered, markers, axes, phase_rows = [], [], [], None, None
    for chrom in range(1, 23):
        directory = dosage / f"chr{chrom}"
        current = pd.read_csv(directory / "samples.tsv", sep="\t")
        with np.load(directory / "dosages.npz") as data:
            if axes is not None and (not axes.equals(current)
                                      or not np.array_equal(phase_rows, data["phase_rows"])):
                raise ValueError("Frozen chromosome axes disagree")
            arrays.append(data["X"])
            altered.append(data["phase_X"])
            phase_rows = data["phase_rows"]
        axes = current
        frame = pd.read_csv(directory / "markers.tsv", sep="\t", dtype={"chrom": str})
        if not frame.chrom.eq(str(chrom)).all():
            raise ValueError("Frozen marker chromosome mismatch")
        markers.extend(list(frame.itertuples(index=False, name=None)))
    x, phase_x = np.concatenate(arrays, axis=1), np.concatenate(altered, axis=1)
    if (x.shape != (4541, 26066) or len(set(markers)) != 26066 or len(phase_rows) != 288
            or axes.duplicated(["group", "sample_id"]).any()
            or not np.array_equal(x[phase_rows], phase_x)
            or set(axes.iloc[phase_rows].group) != {"pure", "pedigree", "real", "array_imputed"}):
        raise ValueError("Frozen dosage coverage or phase-invariance check failed")
    validate_dosages(x)
    return x, phase_x, axes, phase_rows, markers


def prepare(dosage, out):
    if (out / "experiment.json").exists():
        raise ValueError("Fine-profile preparation requires a new experiment")
    previous = json.loads((dosage / "completion.json").read_text())
    verified = json.loads((dosage / "verification/final_verification.json").read_text())
    if (previous["status"] != "complete" or verified["status"] != "complete_verified"
            or verified["completion_sha256"] != adapter.digest(dosage / "completion.json")):
        raise ValueError("Dosage baselines must be complete and independently verified")
    cohort = json.loads((dosage / "cohort.json").read_text())
    adapter.validate_splits(cohort)
    inputs = {**previous["input_sha256"], **previous["output_sha256"]}
    files = [dosage / name for name in
             ("completion.json", "verification/final_verification.json", "cohort.json", "image.txt")]
    files += [Path(cohort["reference"]["labels_tsv"]), dosage / "baseline/metrics.json"]
    files += [Path(batch["samples_tsv"]) for group in cohort["evaluation"].values()
              for batch in group["batches"]]
    files.append(Path(cohort["pure_train"]["labels_tsv"]))
    for chrom in range(1, 23):
        path = dosage / f"chr{chrom}/audit.json"
        audit = json.loads(path.read_text())
        if (audit["status"] != "complete"
                or audit["experiment_sha256"] != adapter.digest(dosage / "experiment.json")):
            raise ValueError(f"Invalid parent dosage checkpoint: {path}")
        inputs.update(audit["output_sha256"])
        files.append(path)
    inputs.update(provenance.fingerprint_files(files))
    provenance.verify_hashes(inputs)
    implementation = [
        Path(__file__), Path(adapter.__file__), Path(provenance.__file__),
        PROJECT / "pedigree_benchmark/scripts/screen_model_families.py",
        PROJECT / "run_fine_profile_screen.sh", PROJECT / "tests/test_fine_profiles.py",
        PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    settings = dict(
        status="configured", models=MODES, sources=SOURCES,
        profile_labels="All100 existing reference fine populations; never evaluation labels",
        shrinkage="n/(n+k)*fine_AF + k/(n+k)*broad_AF; k in0,10,50 individuals",
        informative_mask="Identical to broad baseline: reference broad AF range>=0.03",
        marker_panel="Frozen26066 first markers per100kb physical bin; not r-squared LD pruning",
        optimization=dict(method="SLSQP", ftol=1e-12, maxiter=1000, max_simplex_dual_gap=1e-6),
        initial_weights="Equal total mass per broad source, equal profiles within source",
        objective="Unweighted squared allele-frequency error, sum(q)=1, q>=0; no new intercept",
        no_supervised_training=True, no_evaluation_in_fit=True,
        implementation_sha256=provenance.fingerprint_files(implementation), input_sha256=inputs,
        baselines_immutable=True, production_promotion=False)
    adapter.write_json(out / "experiment.json", settings)
    shutil.copyfile(dosage / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(dosage / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared three fixed fine-profile models; four dosage baselines retained read-only", flush=True)


def configuration(out):
    settings = json.loads((out / "experiment.json").read_text())
    if settings["models"] != MODES or settings["sources"] != SOURCES:
        raise ValueError("Frozen fine-profile settings changed")
    provenance.verify_hashes(settings["implementation_sha256"])
    provenance.verify_hashes(settings["input_sha256"])
    return settings


def run(dosage, out, smoke=False):
    settings = configuration(out)
    if not smoke and (out / "completion.json").exists():
        raise ValueError("Completed fine-profile experiment must not be overwritten")
    x, phase_x, axes, phase_rows, markers = matrices(dosage)
    cohort = json.loads((out / "cohort.json").read_text())
    labels = pd.read_csv(cohort["reference"]["labels_tsv"], sep="\t")
    training = pd.read_csv(cohort["pure_train"]["labels_tsv"], sep="\t")
    if smoke:
        ids = training.groupby("source_group", sort=True).first().sample_id.tolist()
        selected = axes.group.eq("clean_training") & axes.sample_id.isin(ids)
        if selected.sum() != 7:
            raise ValueError("Training-only compatibility sample coverage changed")
        samples = x[selected]
    else:
        samples = x[phase_rows]
    results, outputs, start = {}, [], time.monotonic()
    with threadpool_limits(limits=6):
        for mode, strength in MODES.items():
            directory = out / "candidates" / mode
            checkpoint = directory / "audit.json"
            if not smoke and checkpoint.exists():
                audit = json.loads(checkpoint.read_text())
                if (audit["status"] != "complete"
                        or audit["experiment_sha256"] != adapter.digest(out / "experiment.json")):
                    raise ValueError(f"Incompatible fine-profile checkpoint: {checkpoint}")
                provenance.verify_hashes(audit["output_sha256"])
                results[mode] = audit["result"]
                outputs += [Path(path) for path in audit["output_sha256"]] + [checkpoint]
                continue
            bundle = build_profiles(x, axes, labels, strength)
            if len(bundle["profiles"]) != 100:
                raise ValueError("Expected all100 frozen fine populations")
            bundle.update(marker_keys=markers, experiment_sha256=adapter.digest(out / "experiment.json"))
            predictions, weights, diagnostics = predict(bundle, samples)
            if smoke:
                print(f"{mode}: training-only smoke passed, maxgap="
                      f"{diagnostics.simplex_dual_gap.max():.3g}", flush=True)
                continue
            if directory.exists():
                raise ValueError(f"Incomplete candidate directory requires inspection: {directory}")
            directory.mkdir(parents=True)
            path = directory / "model.pkl"
            with path.with_suffix(".pending").open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            path.with_suffix(".pending").replace(path)
            with path.open("rb") as stream:
                restored = pickle.load(stream)
            reloaded, _, _ = predict(restored, samples)
            changed, _, _ = predict(restored, phase_x)
            reload_delta = float(np.max(np.abs(reloaded - predictions)))
            phase_delta = float(np.max(np.abs(reloaded - changed)))
            if reload_delta > 1e-10 or phase_delta > 1e-10:
                raise ValueError("Saved model or phase-invariance prediction check failed")
            meta = axes.iloc[phase_rows].rename(columns={"group": "dataset"}).reset_index(drop=True)
            populations = pd.DataFrame(weights, columns=[row["profile_id"] for row in bundle["profiles"]])
            pd.concat([meta, populations], axis=1).to_csv(
                directory / "population_weights.tsv", sep="\t", index=False)
            pd.concat([meta, diagnostics], axis=1).to_csv(
                directory / "solver_diagnostics.tsv", sep="\t", index=False)
            pd.DataFrame(bundle["profiles"]).to_csv(directory / "profiles.tsv", sep="\t", index=False)
            for index, source in enumerate(SOURCES):
                meta["pred_" + source] = predictions[:, index]
            array = meta[meta.dataset.eq("array_imputed")].copy()
            if len(array) != 5:
                raise ValueError("Incomplete imputed evaluation coverage")
            array["TV_distance"] = np.abs(
                array[["pred_" + source for source in SOURCES]].to_numpy()
                - [.5, .5, 0, 0, 0, 0, 0]).sum(axis=1) / 2
            array.to_csv(directory / "array_imputed_predictions.tsv", sep="\t", index=False)
            meta[~meta.dataset.eq("array_imputed")].to_csv(
                directory / "predictions.tsv", sep="\t", index=False)
            score(out, mode)
            metrics = json.loads((directory / "metrics.json").read_text())
            results[mode] = dict(
                pure_mean_tv=metrics["pure"]["mean_tv"], pedigree_mean_tv=metrics["pedigree"]["mean_tv"],
                array_imputed_mean_tv=float(array.TV_distance.mean()), real=metrics["real"],
                unchanged_gates_pass=metrics["all_gates_pass"],
                failed_gates=[key for key, value in metrics["gates"].items() if not value],
                max_roundtrip_prediction_error=reload_delta, max_phase_prediction_error=phase_delta,
                max_simplex_dual_gap=float(diagnostics.simplex_dual_gap.max()))
            paths = [path for path in directory.iterdir() if path.is_file()]
            adapter.write_json(checkpoint, dict(
                status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
                result=results[mode], output_sha256=provenance.fingerprint_files(paths)))
            outputs += paths + [checkpoint]
    if smoke:
        adapter.write_json(out / "smoke.json", dict(
            status="complete", training_only=True, n_samples=7,
            experiment_sha256=adapter.digest(out / "experiment.json")))
        return
    previous = json.loads((dosage / "completion.json").read_text())
    rows = [dict(role=role, candidate=name, **{key: value[key] for key in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv", "unchanged_gates_pass")})
        for role, group in (("retained_baseline", previous["results"]), ("candidate", results))
        for name, value in group.items()]
    pd.DataFrame(rows).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    adapter.write_json(out / "completion.json", dict(
        status="complete", results=results, retained_dosage_baselines=previous["results"],
        previous_haplotype_model_results=previous["previous_haplotype_model_results"],
        n_candidates=3, n_profiles=100, n_markers=len(markers), n_chromosomes=22,
        n_reference=3447, evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        no_supervised_training=True, production_promotion=False, independent_final_test=False,
        elapsed_seconds=time.monotonic() - start,
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "comparison.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "smoke", "run"))
    parser.add_argument("--dosage", type=Path, default=Path("/dosage"))
    parser.add_argument("--outdir", type=Path, default=Path("/fine"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.dosage, args.outdir)
    else:
        run(args.dosage, args.outdir, smoke=args.stage == "smoke")


if __name__ == "__main__":
    main()
