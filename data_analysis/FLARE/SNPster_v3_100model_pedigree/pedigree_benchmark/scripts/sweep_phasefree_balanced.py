#!/usr/bin/env python3
"""Training-only founder-mosaic augmentation with unchanged donor-purged validation."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations, permutations
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from phasefree_validation import SOURCES, TRUTH, fold_weights, ranking_key, summarize
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import score
from sweep_phasefree_calibration import checkpoint, configuration, hgb, read_table
from sweep_phasefree_rich import (
    CALIBRATION, CONTROL, canonical_features, genotype_features, predict, array_evaluate,
)


PROJECT = Path(__file__).resolve().parents[2]
RICH = Path("/rich")
DOSAGE = Path("/dosage")
PROTOCOL = PROJECT / "config/phasefree_balanced_protocol.json"


def templates(settings):
    result = []
    for scenario, size, ordered in (("f1", 2, False), ("asymmetric", 2, True),
                                   ("minor", 2, True), ("three_way", 3, False),
                                   ("four_way", 4, False)):
        groups = permutations(SOURCES, size) if ordered else combinations(SOURCES, size)
        for sources in groups:
            variants = [(list(sources), [.5, .5], "")]
            if scenario == "asymmetric":
                variants = [(list(sources), [.75, .25], "")]
            elif scenario == "minor":
                variants = [(list(sources), [1 - fraction, fraction], f":{fraction:.3f}")
                            for fraction in settings["minor_fractions"]]
            elif scenario == "three_way":
                variants = [(list(sources[i:] + sources[:i]), [.5, .25, .25], "") for i in range(3)]
            elif scenario == "four_way":
                variants = [(list(sources), [.25] * 4, "")]
            count = settings["initial_replicates"][scenario]
            if count < len(variants) or count % len(variants):
                raise ValueError("Initial replicates must balance every role/fraction within a stratum")
            for ordered_sources, weights, suffix in variants:
                result.append(dict(
                    scenario=scenario, sources=ordered_sources, weights=weights,
                    stratum=scenario + ":" + "+".join(ordered_sources) + suffix,
                    replicates=count // len(variants)))
    return result


def new_plan(template, labels, rng, index):
    sources, donors = template["sources"], []
    for source in sources:
        eligible = labels[labels.source_group.eq(source)]
        if eligible.empty or eligible.fine.isna().any():
            raise ValueError(f"No labelled training founders for {source}")
        fine = rng.choice(sorted(eligible.fine.unique()))
        donors.append(str(rng.choice(sorted(eligible.loc[eligible.fine.eq(fine), "sample_id"]))))
    return dict(sample_id=f"aug{index:06d}", scenario=template["scenario"], stratum=template["stratum"],
                sources=sources, donors=donors, target_weights=list(template["weights"]),
                segment_seed=int(rng.integers(0, 2 ** 32)))


def safe_indices(plans, held):
    held = set(held)
    return np.asarray([i for i, plan in enumerate(plans) if not set(plan["donors"]) & held], dtype=int)


def make_plans(labels, splits, settings):
    if (labels.sample_id.duplicated().any() or set(labels.source_group) != set(SOURCES)
            or labels.fine.isna().any() or labels.fine.astype(str).str.strip().eq("").any()):
        raise ValueError("Incomplete or duplicate pure-training donor labels")
    rng, plans = np.random.default_rng(settings["seed"]), []
    specifications = templates(settings)
    for template in specifications:
        for _ in range(template["replicates"]):
            plans.append(new_plan(template, labels, rng, len(plans)))
    initial = len(plans)
    if initial > settings["maximum_mosaics"]:
        raise ValueError("Initial augmentation exceeds the fixed mosaic budget")
    for split in splits:
        held = set(split["heldout_donors"])
        available = labels[~labels.sample_id.isin(held)]
        counts = pd.Series([plans[i]["stratum"] for i in safe_indices(plans, held)]).value_counts()
        for template in specifications:
            needed = max(0, settings["minimum_safe_per_stratum_per_fold"]
                         - int(counts.get(template["stratum"], 0)))
            for _ in range(needed):
                if len(plans) >= settings["maximum_mosaics"]:
                    raise ValueError("Donor-safe augmentation coverage exceeds the fixed mosaic budget")
                plans.append(new_plan(template, available, rng, len(plans)))
    return plans, initial


def segments_for(plan, labels, spans, block_bp):
    if block_bp <= 0 or len(set(plan["donors"])) != len(plan["donors"]):
        raise ValueError("Invalid mosaic block size or reused founder")
    lookup = {row.sample_id: (i, row.source_group) for i, row in enumerate(labels.itertuples())}
    donors = []
    for donor, source in zip(plan["donors"], plan["sources"], strict=True):
        if donor not in lookup or lookup[donor][1] != source:
            raise ValueError("Mosaic donor identity/source mismatch")
        donors.append(lookup[donor][0])
    weights = np.asarray(plan["target_weights"], float)
    if (weights.shape != (len(donors),) or not np.isfinite(weights).all()
            or np.any(weights <= 0) or not np.isclose(weights.sum(), 1, atol=1e-12, rtol=0)):
        raise ValueError("Invalid mosaic source weights")
    chromosomes = sorted(map(int, spans))
    if any(spans[str(chrom)][1] <= spans[str(chrom)][0] for chrom in chromosomes):
        raise ValueError("Invalid physical chromosome span")
    rng, segments = np.random.default_rng(plan["segment_seed"]), []
    if plan["scenario"] == "f1":
        if len(donors) != 2 or not np.array_equal(weights, [.5, .5]):
            raise ValueError("F1 mosaics require two equally weighted founders")
        for chrom in chromosomes:
            start, end = spans[str(chrom)]
            for hap in range(2):
                segments.append((chrom, hap, start, end, donors[hap], int(rng.integers(2)),
                                 SOURCES.index(plan["sources"][hap])))
        return np.asarray(segments, dtype=np.int64)
    total = 2 * sum(spans[str(chrom)][1] - spans[str(chrom)][0] for chrom in chromosomes)
    boundaries = np.rint(np.r_[0., np.cumsum(weights)] * total).astype(np.int64)
    boundaries[-1] = total
    tiles = [(i, min(block_bp, count - start))
             for i, count in enumerate(np.diff(boundaries))
             for start in range(0, int(count), block_bp)]
    rng.shuffle(tiles)
    tile_index, tile_used = 0, 0
    orientations = rng.integers(2, size=(len(chromosomes), len(donors)))
    for chrom_index, chrom in enumerate(chromosomes):
        start, end = spans[str(chrom)]
        for hap in range(2):
            position = start
            while position < end:
                source, size = tiles[tile_index]
                length = min(end - position, size - tile_used)
                # Complementary copies avoid artificial autozygosity when a source uses one donor.
                donor_hap = int(orientations[chrom_index, source]) ^ hap
                segments.append((chrom, hap, position, position + length, donors[source],
                                 donor_hap, SOURCES.index(plan["sources"][source])))
                position += length
                tile_used += length
                if tile_used == size:
                    tile_index += 1
                    tile_used = 0
    if tile_index != len(tiles) or tile_used:
        raise ValueError("Mosaic tiles do not cover the physical genome exactly")
    return np.asarray(segments, dtype=np.int64)


def render(segments, founders, keys, spans):
    if (segments.ndim != 2 or segments.shape[1] != 7 or not len(segments)
            or not set(segments[:, 0]) <= set(map(int, spans))
            or not np.isin(segments[:, 1], [0, 1]).all()):
        raise ValueError("Invalid mosaic segment axes or chromosome/copy identifiers")
    values = np.full((2, len(keys)), 255, dtype=np.uint8)
    source_bp = np.zeros(len(SOURCES), dtype=np.int64)
    for chrom in sorted(map(int, spans)):
        indices = np.asarray([i for i, key in enumerate(keys) if key[0] == str(chrom)], dtype=int)
        positions = np.asarray([keys[i][1] for i in indices])
        if not len(indices) or np.any(np.diff(positions) <= 0):
            raise ValueError("Mosaic markers must be nonempty and uniquely ordered by chromosome")
        start, end = spans[str(chrom)]
        if positions[0] < start or positions[-1] >= end:
            raise ValueError("Mosaic markers lie outside physical truth spans")
        for hap in range(2):
            selected = segments[(segments[:, 0] == chrom) & (segments[:, 1] == hap)]
            previous = start
            for _, _, left, right, donor, donor_hap, source in selected:
                if (left != previous or right <= left or right > end or donor_hap not in (0, 1)
                        or not 0 <= donor < len(founders) or not 0 <= source < len(SOURCES)):
                    raise ValueError("Mosaic segment gap, overlap or invalid founder/source")
                chosen = indices[np.searchsorted(positions, left):np.searchsorted(positions, right)]
                values[hap, chosen] = founders[donor, donor_hap, chosen]
                source_bp[source] += right - left
                previous = right
            if previous != end:
                raise ValueError("Mosaic segments do not cover the full chromosome copy")
    if not np.isin(values, [0, 1]).all():
        raise ValueError("Mosaic has unassigned/nonbiallelic inherited alleles")
    return values.sum(axis=0, dtype=np.uint8), source_bp / source_bp.sum()


def augmentation_weights(plans, profile):
    frame = pd.DataFrame(plans)
    if set(frame.scenario) != set(profile) or not np.isclose(sum(profile.values()), 1):
        raise ValueError("Incomplete augmentation scenario weights")
    counts = frame.groupby("stratum").sample_id.transform("size").to_numpy(float)
    strata = frame.groupby("scenario").stratum.transform("nunique").to_numpy(float)
    weights = frame.scenario.map(profile).to_numpy(float) / (counts * strata)
    if np.any(weights <= 0) or not np.isclose(weights.sum(), 1, atol=1e-12, rtol=0):
        raise ValueError("Invalid balanced augmentation weights")
    return weights


def fit_models(x, y, weights, aug_x, aug_y, plans, protocol):
    clean, extra = canonical_features(x, protocol)[:, :14], canonical_features(aug_x, protocol)[:, :14]
    if (y.shape != (len(x), 7) or aug_y.shape != (len(aug_x), 7) or len(plans) != len(aug_x)
            or weights.shape != (len(x),) or not np.isfinite(weights).all()
            or np.any(weights <= 0) or not np.isclose(weights.sum(), 1, atol=1e-12, rtol=0)):
        raise ValueError("Invalid original/augmented training axes or weights")
    for truth in (y, aug_y):
        if not np.isfinite(truth).all() or np.any(truth < 0) or not np.allclose(truth.sum(axis=1), 1):
            raise ValueError("Invalid physical ancestry training truth")
    models = {}
    for name, parameters in protocol["models"].items():
        fraction = parameters["augmentation_fraction"]
        if not 0 <= fraction < 1:
            raise ValueError("Invalid augmentation mass")
        if fraction:
            extra_w = augmentation_weights(plans, protocol["augmentation"]["profiles"][parameters["profile"]])
            train_x, train_y = np.r_[clean, extra], np.r_[y, aug_y]
            train_w = np.r_[(1 - fraction) * weights, fraction * extra_w]
        else:
            train_x, train_y, train_w = clean, y, weights
        head = hgb(protocol["head"], protocol).fit(train_x, train_y, sample_weight=train_w)
        models[name] = dict(head=head, parameters=parameters, columns=list(range(14)), protocol=protocol,
                            n_original=len(x), n_augmented=len(aug_x) if fraction else 0)
    return models


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Use a fresh balanced-mixture experiment")
    parent = json.loads((RICH / "experiment.json").read_text())
    provenance.verify_hashes(parent["implementation_sha256"])
    provenance.verify_hashes(parent["input_sha256"])
    remapped = {str(RICH / Path(path).relative_to("/sweep")): digest
                for path, digest in parent["prepared_output_sha256"].items()}
    provenance.verify_hashes(remapped)
    checked = json.loads((RICH / "verification.json").read_text())
    if (checked["status"] != "complete_verified_cv_and_development"
            or checked["cv_completion_sha256"] != adapter.digest(RICH / "cv_completion.json")):
        raise ValueError("Require independently verified parent feature caches")
    settings = json.loads(PROTOCOL.read_text())
    protocol = {**parent["protocol"], "seed": settings["decoder_seed"], "augmentation": settings,
                "models": settings["models"], "n_finalists": settings["n_finalists"]}
    inputs = dict(parent["input_sha256"])
    inputs.update(remapped)
    inputs.update(provenance.fingerprint_files([
        RICH / "experiment.json", RICH / "verification.json", RICH / "cv_completion.json"]))
    out.mkdir(parents=True, exist_ok=True)
    for name in ("training.tsv", "training_features.npy", "feature_axes.tsv", "all_features.npy",
                 "folds.json", "donors.json", "cohort.json"):
        shutil.copyfile(RICH / name, out / name)
    (out / "baseline").mkdir()
    shutil.copyfile(RICH / "baseline/metrics.json", out / "baseline/metrics.json")
    cohort = json.loads((out / "cohort.json").read_text())
    labels = read_table(cohort["pure_train"]["labels_tsv"])
    metadata = read_table(out / "training.tsv")
    if len(labels) != 307 or set(labels.sample_id) != set(metadata.loc[metadata.track.eq("pure"), "sample_id"]):
        raise ValueError("Augmentation must use exactly the original307 pure-training people")
    vcf = Path(cohort["pure_train"]["vcf"])
    original_audit = json.loads((DOSAGE / "chr1/audit.json").read_text())
    if adapter.digest(vcf) != original_audit["input_sha256"][str(vcf)]:
        raise ValueError("Original phased training genotypes changed")
    inputs.update(provenance.fingerprint_files([vcf, Path(cohort["pure_train"]["labels_tsv"])]))
    with Path(parent["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    keys = original["marker_keys"]
    founders = np.empty((len(labels), 2, len(keys)), dtype=np.uint8)
    for chrom in range(1, 23):
        indices = [i for i, key in enumerate(keys) if key[0] == str(chrom)]
        haplotypes = adapter.load_haplotypes(vcf, labels.sample_id.tolist(), chrom, [keys[i] for i in indices])
        if not np.isin(haplotypes, [0, 1]).all():
            raise ValueError("Missing/nonbiallelic founder haplotypes")
        founders[:, :, indices] = haplotypes.reshape(len(labels), 2, len(indices))
        directory = DOSAGE / f"chr{chrom}"
        axes = read_table(directory / "samples.tsv")
        lookup = {row.sample_id: i for i, row in enumerate(axes.itertuples()) if row.group == "clean_training"}
        with np.load(directory / "dosages.npz") as data:
            np.testing.assert_array_equal(founders[:, :, indices].sum(axis=1),
                                          data["X"][[lookup[identity] for identity in labels.sample_id]])
        print(f"Founder chr{chrom}: complete phased genotypes match frozen dosage cache.", flush=True)
    labels.to_csv(out / "founder_labels.tsv", sep="\t", index=False)
    np.savez_compressed(out / "founders.npz", haplotypes=founders, sample_ids=labels.sample_id.to_numpy(str))
    splits = json.loads((out / "folds.json").read_text())
    plans, initial_count = make_plans(labels, splits, settings)
    all_x, all_y, all_segments = [], [], []
    for i, plan in enumerate(plans):
        segments = segments_for(plan, labels, cohort["chromosome_spans"], settings["block_bp"])
        values, truth = render(segments, founders, keys, cohort["chromosome_spans"])
        all_x.append(values)
        all_y.append(truth)
        all_segments.append(np.column_stack([np.full(len(segments), i), segments]))
    x, truth = np.asarray(all_x, dtype=np.uint8), np.asarray(all_y)
    adapter.write_json(out / "augmentation_plans.json", plans)
    np.savez_compressed(out / "augmentation_genotypes.npz", X=x, truth=truth,
                        segments=np.concatenate(all_segments),
                        sample_ids=np.asarray([plan["sample_id"] for plan in plans]))
    with threadpool_limits(limits=6):
        features = genotype_features(original["feature_models"], x, protocol)
    np.save(out / "augmentation_features.npy", features)
    prepared = [out / name for name in (
        "training.tsv", "training_features.npy", "feature_axes.tsv", "all_features.npy", "folds.json",
        "donors.json", "cohort.json", "baseline/metrics.json", "founder_labels.tsv", "founders.npz",
        "augmentation_plans.json", "augmentation_genotypes.npz", "augmentation_features.npy")]
    implementation = dict(parent["implementation_sha256"])
    implementation.update(provenance.fingerprint_files([
        Path(__file__), Path(__file__).with_name("verify_phasefree_balanced.py"), PROTOCOL,
        PROJECT / "run_phasefree_balanced_sweep.sh", PROJECT / "tests/test_phasefree_balanced.py"]))
    adapter.write_json(out / "experiment.json", dict(
        status="configured", protocol=protocol, targets=parent["targets"], gates=parent["gates"],
        original_model_path=parent["original_model_path"], implementation_sha256=implementation,
        input_sha256=inputs, prepared_output_sha256=provenance.fingerprint_files(prepared),
        initial_mosaics=initial_count, n_mosaics=len(plans), n_unique_founders=307,
        n_original_oof_people=799, n_donor_purged_folds=114,
        augmentation_is_training_only=True, independent_final_test=False, production_promotion=False))
    from verify_phasefree_balanced import verify_bank
    verify_bank(out)
    print(f"Prepared {len(plans)} training-only mosaics; no new validation people.", flush=True)


def fold_models(out, split):
    config = configuration(out)
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    plans = json.loads((out / "augmentation_plans.json").read_text())
    kept = safe_indices(plans, split["heldout_donors"])
    if set(plans[i]["stratum"] for i in kept) != {plan["stratum"] for plan in plans}:
        raise ValueError("Founder purging removed an entire augmentation stratum")
    with np.load(out / "augmentation_genotypes.npz") as data:
        aug_y = data["truth"][kept]
    aug_x = np.load(out / "augmentation_features.npy")[kept]
    train = np.flatnonzero(metadata.sample_id.isin(split["train_ids"]))
    models = fit_models(x[train], metadata.iloc[train][TRUTH].to_numpy(float),
                        fold_weights(metadata.iloc[train]), aug_x, aug_y, [plans[i] for i in kept], config["protocol"])
    return models, kept


def run_fold(out, split):
    directory = out / "folds" / split["fold_id"]
    if checkpoint(directory / "completion.json", out):
        return split["fold_id"]
    if directory.exists():
        raise ValueError(f"Inspect incomplete balanced-mixture fold: {directory}")
    directory.mkdir(parents=True)
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
    started, predictions, stability = time.monotonic(), [], {}
    with threadpool_limits(limits=2):
        models, kept = fold_models(out, split)
        for name, model in models.items():
            values = predict(model, x[test])
            np.testing.assert_array_equal(values, predict(pickle.loads(pickle.dumps(model)), x[test]))
            error = max(float(np.max(np.abs(values - predict(
                model, np.nextafter(x[test], direction))))) for direction in (-np.inf, np.inf))
            if error > 1e-12:
                raise ValueError("Augmented decoder violates the stable feature contract")
            predictions.append(values)
            stability[name] = error
    with np.load(CALIBRATION / "oof_predictions.npz") as prior:
        if prior["sample_ids"].tolist() != metadata.sample_id.tolist():
            raise ValueError("Prior stable control OOF identities changed")
        control = prior["predictions"][prior["candidates"].tolist().index("stable_squared"), test]
    np.testing.assert_allclose(predictions[list(models).index(CONTROL)], control, atol=1e-12, rtol=0)
    np.savez_compressed(directory / "predictions.npz", predictions=np.asarray(predictions),
                        sample_ids=metadata.iloc[test].sample_id.to_numpy(str), candidates=np.asarray(list(models)))
    adapter.write_json(directory / "completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"), fold=split,
        augmented_indices=kept.tolist(), elapsed_seconds=time.monotonic() - started,
        max_one_ulp_prediction_difference=stability, prior_control_reproduced=True,
        output_sha256=provenance.fingerprint_files([directory / "predictions.npz"])))
    return split["fold_id"]


def cross_validate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "cv_completion.json", out):
        return
    splits = json.loads((out / "folds.json").read_text())
    with ProcessPoolExecutor(max_workers=2) as pool:
        tasks = [pool.submit(run_fold, out, split) for split in splits]
        for i, future in enumerate(as_completed(tasks), 1):
            print(f"Balanced CV {i}/{len(splits)}: {future.result()} verified", flush=True)
    metadata = read_table(out / "training.tsv")
    candidates = list(config["protocol"]["models"])
    values = np.empty((len(candidates), len(metadata), 7))
    lookup, seen = {identity: i for i, identity in enumerate(metadata.sample_id)}, set()
    for split in splits:
        directory = out / "folds" / split["fold_id"]
        if not checkpoint(directory / "completion.json", out):
            raise ValueError("Missing augmented-training fold")
        with np.load(directory / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if seen & set(ids) or data["candidates"].tolist() != candidates:
                raise ValueError("Balanced OOF identity/candidate overlap")
            values[:, [lookup[identity] for identity in ids]] = data["predictions"]
            seen.update(ids)
    if seen != set(lookup) or len(seen) != 799:
        raise ValueError("Augmentation changed the original OOF validation cohort")
    reports = {name: summarize(metadata, values[i], config["targets"], config["gates"])
               for i, name in enumerate(candidates)}
    order = sorted(candidates, key=lambda name: ranking_key(name, reports[name]))
    adapter.write_json(out / "cv_metrics.json", reports)
    adapter.write_json(out / "selection.json", dict(
        ranking=order, finalists=[name for name in order if name != CONTROL][:2],
        development_metrics_used=False, synthetic_people_in_validation=0))
    np.savez_compressed(out / "oof_predictions.npz", predictions=values,
                        sample_ids=metadata.sample_id.to_numpy(str), candidates=np.asarray(candidates))
    comparison = [dict(
        candidate=name, balanced_mean_tv=reports[name]["unified"]["balanced_mean_tv"],
        pure_mean_tv=reports[name]["pure"]["mean_tv"], pedigree_mean_tv=reports[name]["pedigree"]["mean_tv"],
        pedigree_within_10pp=reports[name]["pedigree"]["fraction_all_components_within_10pp"],
        pedigree_within_5pp=reports[name]["pedigree"]["fraction_all_components_within_5pp"]) for name in order]
    pd.DataFrame(comparison).to_csv(out / "cv_comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "cv_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        output_sha256=provenance.fingerprint_files([out / name for name in (
            "cv_metrics.json", "selection.json", "oof_predictions.npz", "cv_comparison.tsv")])))
    print(pd.DataFrame(comparison).to_string(index=False), flush=True)


def evaluate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "development_completion.json", out):
        return
    if not checkpoint(out / "cv_completion.json", out):
        raise ValueError("Lock balanced training-only selection before development evaluation")
    selection = json.loads((out / "selection.json").read_text())
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    plans, aug_x = json.loads((out / "augmentation_plans.json").read_text()), np.load(out / "augmentation_features.npy")
    with np.load(out / "augmentation_genotypes.npz") as data:
        aug_y = data["truth"]
    with Path(config["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    axes, features = read_table(out / "feature_axes.tsv"), np.load(out / "all_features.npy")
    chosen = axes.group.isin(["pure", "pedigree", "real"]).to_numpy()
    old_array = axes.group.eq("array_imputed").to_numpy()
    if old_array.sum() != 5:
        raise ValueError("Changed original array-control cohort")
    paths, comparison = [], []
    with threadpool_limits(limits=2):
        models = fit_models(x, metadata[TRUTH].to_numpy(float), metadata.weight_clean.to_numpy(float),
                            aug_x, aug_y, plans, config["protocol"])
        for name, model in models.items():
            directory = out / "candidates" / name
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / "model.pkl"
            if path.exists():
                raise ValueError("Inspect incomplete augmented full fit before retrying")
            bundle = dict(kind="phasefree_balanced_training", sources=SOURCES, marker_keys=original["marker_keys"],
                          feature_models=original["feature_models"], decoder=model, protocol=config["protocol"],
                          experiment_sha256=adapter.digest(out / "experiment.json"),
                          input_contract="Complete pre-QC hard GT; same frozen reference transforms and 14 decoder inputs")
            with path.open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            paths.append(path)
            if name not in [CONTROL, *selection["finalists"]]:
                continue
            prediction = predict(model, features[chosen])
            np.testing.assert_array_equal(prediction, predict(pickle.loads(pickle.dumps(model)), features[chosen]))
            frame = pd.DataFrame(prediction, columns=[f"pred_{source}" for source in SOURCES])
            frame.insert(0, "sample_id", axes.loc[chosen, "sample_id"].to_numpy())
            frame.insert(0, "dataset", axes.loc[chosen, "group"].to_numpy())
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
                np.tile([.5, .5, 0, 0, 0, 0, 0], (5, 1)), predict(model, features[old_array]), config["targets"])
            report["old_five_eur_asian_array_controls"] = array_report
            adapter.write_json(directory / "standard_validation.json", report)
            metrics = json.loads((directory / "metrics.json").read_text())
            comparison.append(dict(
                candidate=name, pure_mean_tv=report["pure"]["mean_tv"],
                pedigree_mean_tv=report["pedigree"]["mean_tv"],
                pedigree_within_10pp=report["pedigree"]["fraction_all_components_within_10pp"],
                pedigree_within_5pp=report["pedigree"]["fraction_all_components_within_5pp"],
                array5_mean_tv=array_report["mean_tv"], unchanged_gates_pass=metrics["all_gates_pass"]))
            paths += [directory / file for file in ("predictions.tsv", "metrics.json", "standard_validation.json")]
            paths += list(directory.glob("*scored.tsv"))
    pd.DataFrame(comparison).to_csv(out / "development_comparison.tsv", sep="\t", index=False)
    paths.append(out / "development_comparison.tsv")
    adapter.write_json(out / "development_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        selection_sha256=adapter.digest(out / "selection.json"),
        output_sha256=provenance.fingerprint_files(paths), production_promotion=False))
    print(pd.DataFrame(comparison).to_string(index=False), flush=True)


def run(out):
    config = configuration(out, full=True)
    bank = json.loads((out / "augmentation_verification.json").read_text())
    if bank["experiment_sha256"] != adapter.digest(out / "experiment.json"):
        raise ValueError("Require the independently verified training-only augmentation bank")
    adapter.write_json(out / "sweep_state.json", dict(status="cross_validating"))
    cross_validate(out)
    adapter.write_json(out / "sweep_state.json", dict(status="fitting_locked_finalists"))
    evaluate(out)
    from verify_phasefree_balanced import verify
    verify(out)
    arrays, deadline = Path("/validation"), time.monotonic() + 48 * 3600
    while not (arrays / "verification.json").exists():
        state = json.loads((arrays / "durable_scheduler.json").read_text())
        if state["status"] == "failed":
            raise RuntimeError("Shared array dependency failed; keep locked augmented models and recover it")
        if time.monotonic() >= deadline:
            raise TimeoutError("Shared array verification unavailable after48 hours")
        adapter.write_json(out / "sweep_state.json", dict(
            status="waiting_for_shared_array_verification", development_verified=True, production_promotion=False))
        time.sleep(60)
    adapter.write_json(out / "sweep_state.json", dict(status="evaluating_shared_arrays"))
    array_evaluate(out, arrays)
    from verify_phasefree_rich import verify_arrays
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
        run_fold(args.output, split)
        path = args.output / "folds" / split["fold_id"] / "completion.json"
        if not checkpoint(path, args.output):
            raise ValueError("Balanced smoke fold lacks verified predictions")
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
