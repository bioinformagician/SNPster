#!/usr/bin/env python3
"""Two fixed training-only phase-augmentation arms for frozen LR+CRF bases."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from pilot_matched_imputation import fingerprint_files, verify_hashes
from screen_model_families import score, verify


PROJECT = Path(__file__).resolve().parents[2]
BASE = "gnomix_fast_050"
ARMS = {"phase_blocks_050": .5, "phase_blocks_200": 2.}


def phase_view(x, labels, cm, mean_block_cm, seed):
    if (x.shape != labels.shape or len(x) % 2 or x.shape[1] != len(cm)
            or mean_block_cm <= 0 or not np.all(np.isfinite(cm))
            or np.any(np.diff(cm) < 0) or np.any(x < 0) or np.any(x > 1)):
        raise ValueError("Invalid complete paired genotypes, labels or genetic coordinates")
    paired = x.reshape(-1, 2, x.shape[1])
    ancestry = labels.reshape(paired.shape)
    rng = np.random.default_rng(seed)
    # Parity of Poisson switches, rather than at most one event across a long marker gap.
    probability = -np.expm1(-2 * np.diff(cm) / mean_block_cm) / 2
    transitions = rng.random((len(paired), len(cm) - 1)) < probability
    orientation = np.concatenate(
        [rng.integers(0, 2, (len(paired), 1)), transitions], axis=1).cumsum(axis=1) % 2
    mask = orientation[:, None, :].astype(bool)
    altered = np.where(mask, paired[:, ::-1], paired).reshape(x.shape)
    truth = np.where(mask, ancestry[:, ::-1], ancestry).reshape(labels.shape)
    if (not np.array_equal(altered.reshape(paired.shape).sum(axis=1), paired.sum(axis=1))
            or not np.array_equal(np.sort(truth.reshape(paired.shape), axis=1),
                                  np.sort(ancestry, axis=1))):
        raise ValueError("Phase augmentation changed diploid genotype or ancestry truth")
    return altered, truth


def stable_seed(sample, chrom):
    text = f"{adapter.SEED}:{sample}:{chrom}".encode()
    return int.from_bytes(hashlib.sha256(text).digest()[:8], "little")


def load_experiment(out):
    config = json.loads((out / "experiment.json").read_text())
    if config["arms"] != ARMS or config["script_sha256"] != adapter.digest(__file__):
        raise ValueError("Frozen experiment implementation or settings changed")
    return config


def verified_checkpoint(path, out):
    if not path.exists():
        return False
    audit = json.loads(path.read_text())
    if (audit["status"] != "complete"
            or audit["experiment_sha256"] != adapter.digest(out / "experiment.json")):
        raise ValueError(f"Incomplete or incompatible checkpoint: {path}")
    verify_hashes(audit["input_sha256"])
    verify_hashes(audit["output_sha256"])
    print(f"Verified checkpoint {path}; no repeated work", flush=True)
    return True


def prepare(screen, pilot, out):
    if (out / "experiment.json").exists():
        raise ValueError("Preparation requires a new isolated result directory")
    cohort = verify(screen)
    counts = adapter.validate_splits(cohort)
    original = json.loads((screen / f"candidates/{BASE}/output_audit.json").read_text())
    if original["status"] != "complete" or adapter.input_hashes(cohort) != original["input_sha256"]:
        raise ValueError("Original candidate provenance changed")
    matched = json.loads((pilot / "summary.json").read_text())
    if matched["status"] != "complete" or matched["n_chromosomes"] != 22:
        raise ValueError("Matched imputation validation must be complete")
    train_root = PROJECT / "results_v3_pedigree_training_screen_01/cohort/pedigree/plans/plans"
    development_root = PROJECT / "results_v3_pedigree_pilot_10/pedigree/plans/plans"
    train_donors, dev_donors, provenance = set(), set(), []
    for root, target in ((train_root, train_donors), (development_root, dev_donors)):
        plans = sorted(root.glob("batch_*/pedigrees.json"))
        if not plans:
            raise ValueError(f"Original founder plans missing: {root}")
        for path in plans:
            plan = json.loads(path.read_text())
            target.update(item["vcf_sample"] for item in plan["founders"].values())
            provenance.append(path)
    pure_train = {r["sample_id"] for r in adapter.sample_records(cohort["pure_train"]["labels_tsv"])}
    reference = {r["sample_id"] for r in adapter.sample_records(cohort["reference"]["labels_tsv"])}
    validation = {r["sample_id"] for b in cohort["evaluation"]["pure"]["batches"]
                  for r in adapter.sample_records(b["samples_tsv"])}
    if (not train_donors <= pure_train or train_donors & (dev_donors | reference | validation)
            or dev_donors & (pure_train | reference)):
        raise ValueError("Training founder pools are not reference/evaluation disjoint")
    out.mkdir(parents=True, exist_ok=True)
    config = dict(
        status="configured", baseline=BASE, arms=ARMS, seed=adapter.SEED,
        sample_counts=counts, n_training_founder_donors=len(train_donors),
        n_development_founder_donors=len(dev_donors),
        script_sha256=adapter.digest(__file__), adapter_sha256=adapter.digest(adapter.__file__),
        input_sha256={**adapter.input_hashes(cohort),
                      **fingerprint_files(provenance + [screen / "cohort.json",
                                                        pilot / "pilot.json", pilot / "summary.json"])},
        policy="Original training sequences plus one equally weighted switched view per individual; separate CRF for each fixed phase-block length.",
        frozen_base=True, smoother_only=True, no_new_synthetic_genomes=True,
        augmentation_preserves_diploid_genotypes=True,
        ancestry_labels_swapped_with_genotypes=True,
        phase_rates_fit_to_evaluation=False, training_uses_evaluation_genotypes=False,
        crf_hyperparameters="Unchanged upstream defaults; doubled sequence count changes effective data/regularization balance.",
        augmentation_limitation="Genotype-preserving phase perturbations, not actual array-imputed training or a simulation of genotype-call errors.",
        evaluation_tracks=["pure", "pedigree", "real", "array_imputed"],
        evaluation_caveat="Previously inspected development cohorts, including five shared-donor controls and one real case; not independent certification.",
        production_promotion=False, ancestry_sources=list(adapter.SOURCES))
    adapter.write_json(out / "experiment.json", config)
    shutil.copyfile(screen / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(screen / "baseline/metrics.json", out / "baseline/metrics.json")
    print(f"Configured two fixed phase-augmentation arms; {len(train_donors)} training donors disjoint",
          flush=True)


def profiles(screen, out, chrom):
    config = load_experiment(out)
    directory = out / f"profiles/chr{chrom}"
    if verified_checkpoint(directory / "audit.json", out):
        return
    if directory.exists():
        raise ValueError("Uncommitted profile directory exists; inspect before rerunning")
    directory.mkdir(parents=True)
    cohort = json.loads((screen / "cohort.json").read_text())
    adapter.load_upstream()
    source_model = screen / f"candidates/{BASE}/models/chr{chrom}.pkl"
    bundle, _ = adapter.saved_chromosome(cohort, adapter.config_for(BASE), chrom,
                                         screen / "candidates" / BASE)
    model, keys, layout = bundle["model"], bundle["keys"], bundle["layout"]
    points = adapter.map_points(cohort["genetic_map"], chrom)
    cm = np.interp(layout["positions"], points[:, 0], points[:, 1])
    training = [(cohort["pure_train"],
                 adapter.sample_records(cohort["pure_train"]["labels_tsv"], True), True)]
    training += [(batch, adapter.sample_records(batch["samples_tsv"]), False)
                 for batch in cohort["pedigree_train"]["batches"]]
    probabilities = {arm: [] for arm in ("original", *ARMS)}
    labels = {arm: [] for arm in probabilities}
    ids_seen, inputs = [], [source_model, source_model.with_suffix(".json"),
                           Path(cohort["genetic_map"]), screen / "cohort.json"]
    start = time.monotonic()
    with threadpool_limits(limits=1):
        for batch, rows, pure in training:
            label_file = batch["labels_tsv"] if pure else batch["samples_tsv"]
            inputs += [Path(batch["vcf"]), Path(label_file)]
            if not pure:
                inputs.append(Path(batch["segments_tsv"]))
            for part in adapter.chunks(rows, 64):
                ids = [row["sample_id"] for row in part]
                ids_seen.extend(ids)
                x = adapter.load_haplotypes(batch["vcf"], ids, chrom, keys)
                if np.any(x > 1):
                    raise ValueError("Augmentation requires complete training genotypes")
                if pure:
                    ancestry = np.repeat(
                        np.repeat([adapter.SOURCES.index(row["source_group"]) for row in part], 2)
                        [:, None], len(keys), axis=1).astype(np.int8)
                else:
                    ancestry = adapter.segment_labels(batch["segments_tsv"], ids, chrom,
                                                       layout["positions"])
                probabilities["original"].append(model.base.predict_proba(x).astype(np.float32))
                labels["original"].append(adapter.window_labels(ancestry, layout))
                for arm, mean in ARMS.items():
                    switched, target = [], []
                    for index, sample in enumerate(ids):
                        view, truth = phase_view(x[2 * index:2 * index + 2],
                                                 ancestry[2 * index:2 * index + 2],
                                                 cm, mean, stable_seed(sample, chrom))
                        switched.append(view)
                        target.append(truth)
                    altered = np.concatenate(switched)
                    probabilities[arm].append(model.base.predict_proba(altered).astype(np.float32))
                    labels[arm].append(adapter.window_labels(np.concatenate(target), layout))
    if len(set(ids_seen)) != len(ids_seen) or len(ids_seen) != 799:
        raise ValueError("Training sample-axis coverage differs from the frozen experiment")
    for arm in probabilities:
        np.savez_compressed(directory / f"{arm}.npz",
                            probabilities=np.concatenate(probabilities[arm]),
                            labels=np.concatenate(labels[arm]))
    adapter.write_json(directory / "samples.json", ids_seen)
    hashes = fingerprint_files(inputs)
    if any(str(path) in config["input_sha256"]
           and expected != config["input_sha256"][str(path)] for path, expected in hashes.items()):
        raise ValueError("Training input differs from the prepared cohort")
    audit = dict(status="complete", chromosome=chrom, n_samples=len(ids_seen),
                 n_markers=len(keys), n_windows=layout["W"], seconds=time.monotonic() - start,
                 experiment_sha256=adapter.digest(out / "experiment.json"),
                 input_sha256=hashes,
                 output_sha256=fingerprint_files(list(directory.glob("*.npz")) +
                                                  [directory / "samples.json"]))
    adapter.write_json(directory / "audit.json", audit)
    print(f"chr{chrom}: training-only original and two perturbed profiles saved", flush=True)


def predict_batch(model, x):
    return adapter.ordered_probabilities(model, model.base.predict_proba(x))


def fit_evaluate(screen, pilot, out, chrom, arm):
    config = load_experiment(out)
    directory = out / f"candidates/{arm}/chr{chrom}"
    if verified_checkpoint(directory / "audit.json", out):
        return
    if directory.exists():
        raise ValueError("Uncommitted fit exists; inspect rather than overwrite")
    directory.mkdir(parents=True)
    source = out / f"profiles/chr{chrom}"
    if not verified_checkpoint(source / "audit.json", out):
        raise ValueError("Verified training profiles are missing")
    adapter.load_upstream()
    cohort = json.loads((screen / "cohort.json").read_text())
    bundle, original_metadata = adapter.saved_chromosome(
        cohort, adapter.config_for(BASE), chrom, screen / "candidates" / BASE)
    layout, keys, model = bundle["layout"], bundle["keys"], bundle["model"]
    base_hash = hashlib.sha256(pickle.dumps(model.base)).hexdigest()
    with np.load(source / "original.npz") as clean, np.load(source / f"{arm}.npz") as noisy:
        b = np.concatenate([clean["probabilities"], noisy["probabilities"]])
        y = np.concatenate([clean["labels"], noisy["labels"]])
    if set(np.unique(y)) != set(range(len(adapter.SOURCES))):
        raise ValueError("Augmented smoother training lacks a source")
    from src.Smooth.models import CRF_Smoother
    model.smooth = CRF_Smoother(layout["W"], len(adapter.SOURCES), model.S,
                                calibrate=False, seed=adapter.SEED, n_jobs=1, verbose=False)
    start = time.monotonic()
    with threadpool_limits(limits=1):
        model.smooth.train(b, y)
    del b, y
    seconds = time.monotonic() - start
    if hashlib.sha256(pickle.dumps(model.base)).hexdigest() != base_hash:
        raise ValueError("Smoother retraining changed the frozen base classifier")
    bundle["robustness"] = dict(arm=arm, mean_phase_block_cm=ARMS[arm],
                                original_model_sha256=original_metadata["model_sha256"])
    model_path = directory / "model.pkl"
    with model_path.with_suffix(".pending").open("wb") as stream:
        pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
    model_path.with_suffix(".pending").replace(model_path)
    with model_path.open("rb") as stream:
        reloaded = pickle.load(stream)
    predictions, max_error, inputs = [], 0., [source / "original.npz", source / f"{arm}.npz",
                                            source / "audit.json", screen / "cohort.json",
                                            out / "cohort.json", out / "baseline/metrics.json"]
    pilot_audit_path = pilot / f"chr{chrom}/audit.json"
    pilot_audit = json.loads(pilot_audit_path.read_text())
    if pilot_audit["status"] != "complete":
        raise ValueError("Array-imputed validation checkpoint is incomplete")
    verify_hashes(pilot_audit["output_sha256"])
    inputs.append(pilot_audit_path)
    tracks = [(name, batch, adapter.sample_records(batch["samples_tsv"]))
              for name, group in cohort["evaluation"].items() for batch in group["batches"]]
    pilot_config = json.loads((pilot / "pilot.json").read_text())
    tracks += [("array_imputed",
                dict(vcf=str(pilot / f"chr{chrom}/{row['sample_id']}/truth_panel.vcf.gz")),
                [dict(sample_id=row["sample_id"])]) for row in pilot_config["samples"]]
    with threadpool_limits(limits=1):
        for track, batch, rows in tracks:
            inputs.append(Path(batch["vcf"]))
            if "samples_tsv" in batch:
                inputs.append(Path(batch["samples_tsv"]))
            for part in adapter.chunks(rows, 64):
                ids = [row["sample_id"] for row in part]
                x = adapter.load_haplotypes(batch["vcf"], ids, chrom, keys)
                before = predict_batch(model, x)
                after = predict_batch(reloaded["model"], x)
                delta = float(np.max(np.abs(before - after)))
                if delta > 1e-10:
                    raise ValueError("Saved model reload changed evaluation posterior tensor")
                max_error = max(max_error, delta)
                marker, physical = adapter.aggregate_probabilities(after, layout)
                for index, sample in enumerate(ids):
                    for estimand, values, denominator in (
                            ("marker", marker[index], len(keys)),
                            ("physical", physical[index], layout["physical_weights"].sum())):
                        predictions.append(dict(
                            dataset=track, sample_id=sample, estimand=estimand, n_markers=len(keys),
                            physical_span=float(layout["physical_weights"].sum()),
                            **dict(zip(adapter.SOURCES, values / denominator))))
    pd.DataFrame(predictions).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    inputs.append(screen / f"candidates/{BASE}/models/chr{chrom}.pkl")
    hashes = fingerprint_files(inputs)
    profile_audit = json.loads((source / "audit.json").read_text())
    verify_hashes(profile_audit["input_sha256"])
    verify_hashes(profile_audit["output_sha256"])
    audit = dict(status="complete", arm=arm, chromosome=chrom, fit_seconds=seconds,
                 base_unchanged=True, model_roundtrip_verified=True,
                 max_posterior_reload_error=max_error, n_individual_predictions=len(predictions) // 2,
                 n_training_sequences=3196, n_validation_array_controls=5,
                 experiment_sha256=adapter.digest(out / "experiment.json"),
                 input_sha256=hashes, output_sha256=fingerprint_files(
                     [model_path, directory / "predictions.tsv"]))
    adapter.write_json(directory / "audit.json", audit)
    print(f"{arm} chr{chrom}: smoother fit {seconds:.1f}s; all evaluation tracks reload-verified",
          flush=True)


def summarize(screen, pilot, out):
    config = load_experiment(out)
    verify_hashes(config["input_sha256"])
    original_crf = json.loads((screen / f"candidates/{BASE}/metrics.json").read_text())
    rows, result = [], {}
    pilot_baseline = json.loads((pilot / "summary.json").read_text())
    for arm in ARMS:
        directory = out / f"candidates/{arm}"
        chromosome_predictions = []
        for chrom in range(1, 23):
            folder = directory / f"chr{chrom}"
            if not verified_checkpoint(folder / "audit.json", out):
                raise ValueError(f"Missing {arm} chr{chrom}")
            frame = pd.read_csv(folder / "predictions.tsv", sep="\t")
            frame["chromosome"] = chrom
            chromosome_predictions.append(frame)
        predictions = pd.concat(chromosome_predictions, ignore_index=True)
        if predictions.duplicated(["dataset", "sample_id", "estimand", "chromosome"]).any():
            raise ValueError("Duplicate chromosome predictions")
        for _, group in predictions.groupby(["dataset", "sample_id", "estimand"]):
            if set(group.chromosome) != set(range(1, 23)) or group.n_markers.sum() != 197260:
                raise ValueError("Incomplete chromosome or marker coverage")
        array_mean = None
        for estimand in ("marker", "physical"):
            frame = predictions[predictions.estimand.eq(estimand)].copy()
            weights = frame.n_markers if estimand == "marker" else frame.physical_span
            for source in adapter.SOURCES:
                frame["pred_" + source] = frame[source] * weights
            frame["denominator"] = weights
            summed = frame.groupby(["dataset", "sample_id"])[[
                *["pred_" + s for s in adapter.SOURCES], "denominator"]].sum()
            cols = ["pred_" + s for s in adapter.SOURCES]
            summed[cols] = summed[cols].div(summed.denominator, axis=0)
            summed = summed.drop(columns="denominator").reset_index()
            array = summed[summed.dataset.eq("array_imputed")].copy()
            if len(array) != 5:
                raise ValueError("Array validation axis is incomplete")
            array["TV_distance"] = np.abs(
                array[cols].to_numpy() - np.array([.5, .5, 0, 0, 0, 0, 0])).sum(axis=1) / 2
            suffix = "" if estimand == "marker" else "_physical"
            array.to_csv(directory / f"array_imputed_predictions{suffix}.tsv", sep="\t", index=False)
            if estimand == "marker":
                array_mean = float(array.TV_distance.mean())
            summed[~summed.dataset.eq("array_imputed")].to_csv(
                directory / f"predictions{suffix}.tsv", sep="\t", index=False)
        score(out, arm)
        metrics = json.loads((directory / "metrics.json").read_text())
        comparisons = dict(
            pure_vs_original_crf=metrics["pure"]["mean_tv"] - original_crf["pure"]["mean_tv"],
            pedigree_vs_original_crf=metrics["pedigree"]["mean_tv"] - original_crf["pedigree"]["mean_tv"],
            array_vs_original_crf=array_mean - pilot_baseline["stage_mean_tv"]["array_imputed"])
        result[arm] = dict(
            pure_mean_tv=metrics["pure"]["mean_tv"], pedigree_mean_tv=metrics["pedigree"]["mean_tv"],
            array_imputed_mean_tv=array_mean, n_array_controls=5,
            unchanged_gates_pass=metrics["all_gates_pass"],
            failed_unchanged_gates=[key for key, value in metrics["gates"].items() if not value],
            real=metrics["real"], real_constraints=metrics["real_constraints"],
            comparisons_to_original_crf=comparisons,
            all_three_means_no_original_crf_regression=all(value <= 1e-10 for value in comparisons.values()))
        rows.append(dict(candidate=arm, **{key: result[arm][key] for key in
                                          ("pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv",
                                           "unchanged_gates_pass")}))
    pd.DataFrame(rows).to_csv(out / "comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "completion.json", dict(
        status="complete", n_candidates=2, n_chromosomes=22, results=result,
        original_crf=dict(pure_mean_tv=original_crf["pure"]["mean_tv"],
                          pedigree_mean_tv=original_crf["pedigree"]["mean_tv"],
                          array_imputed_mean_tv=pilot_baseline["stage_mean_tv"]["array_imputed"]),
        production_promotion=False, independent_final_test=False,
        limitations=config["augmentation_limitation"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "profiles", "fit", "summarize"), required=True)
    parser.add_argument("--screen", type=Path, default=Path("/results"))
    parser.add_argument("--pilot", type=Path, default=Path("/matched"))
    parser.add_argument("--outdir", type=Path, default=Path("/robust"))
    parser.add_argument("--chromosome", type=int, choices=range(1, 23))
    parser.add_argument("--arm", choices=tuple(ARMS))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.screen, args.pilot, args.outdir)
    elif args.stage == "summarize":
        summarize(args.screen, args.pilot, args.outdir)
    elif args.chromosome is None:
        parser.error("--chromosome is required")
    elif args.stage == "profiles":
        profiles(args.screen, args.outdir, args.chromosome)
    elif args.arm is None:
        parser.error("--arm is required for fitting")
    else:
        fit_evaluate(args.screen, args.pilot, args.outdir, args.chromosome, args.arm)


if __name__ == "__main__":
    main()
