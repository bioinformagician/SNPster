#!/usr/bin/env python3
"""Independently verify mosaic inheritance, physical truth and donor-purged model results."""

import json
from pathlib import Path
import pickle

import numpy as np
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from phasefree_validation import SOURCES, TRUTH
from sweep_phasefree_calibration import checkpoint, configuration, predict as calibration_predict, read_table
from verify_phasefree_rich import check_metrics
from verify_phasefree_sweep import independently_balanced_tv


def verify_bank(out):
    config = configuration(out, full=True)
    settings = config["protocol"]["augmentation"]
    labels = read_table(out / "founder_labels.tsv")
    metadata = read_table(out / "training.tsv")
    plans = json.loads((out / "augmentation_plans.json").read_text())
    cohort = json.loads((out / "cohort.json").read_text())
    spans = cohort["chromosome_spans"]
    with Path(config["original_model_path"]).open("rb") as stream:
        frozen = pickle.load(stream)
    keys = frozen["marker_keys"]
    with np.load(out / "founders.npz") as data:
        founders = data["haplotypes"]
        if data["sample_ids"].tolist() != labels.sample_id.tolist():
            raise ValueError("Founder haplotype identities changed")
    if (len(labels) != 307 or len(set(labels.sample_id)) != 307
            or set(labels.sample_id) != set(metadata.loc[metadata.track.eq("pure"), "sample_id"])
            or founders.shape != (307, 2, len(keys)) or not np.isin(founders, [0, 1]).all()):
        raise ValueError("Augmentation did not use exactly the existing complete training founders")
    with np.load(out / "augmentation_genotypes.npz") as data:
        x, truth, segments = data["X"], data["truth"], data["segments"]
        if data["sample_ids"].tolist() != [plan["sample_id"] for plan in plans]:
            raise ValueError("Augmentation genotype/plan identities differ")
    if (x.shape != (len(plans), len(keys)) or truth.shape != (len(plans), 7)
            or segments.ndim != 2 or segments.shape[1] != 8 or not len(segments)
            or np.any(np.diff(segments[:, 0]) < 0)
            or set(segments[:, 0]) != set(range(len(plans)))
            or not set(segments[:, 1]) <= set(map(int, spans))
            or not np.isin(segments[:, 2], [0, 1]).all()):
        raise ValueError("Invalid augmentation genotype/segment axes")
    total_bp = 2 * sum(end - start for start, end in spans.values())
    marker_axes = {
        int(chrom): (np.asarray([i for i, key in enumerate(keys) if key[0] == chrom]),
                     np.asarray([key[1] for key in keys if key[0] == chrom]))
        for chrom in spans}
    source_index = [SOURCES.index(source) for source in labels.source_group]
    for sample_index, plan in enumerate(plans):
        if (len(set(plan["donors"])) != len(plan["donors"])
                or not set(plan["donors"]) <= set(labels.sample_id)):
            raise ValueError("Duplicate or nontraining augmentation donor")
        calls = np.full((2, len(keys)), 255, dtype=np.uint8)
        donor_choices = np.full((2, len(keys)), -1, dtype=np.int16)
        hap_choices = np.full((2, len(keys)), -1, dtype=np.int8)
        assigned_bp = np.zeros(7, dtype=np.int64)
        left, right = np.searchsorted(segments[:, 0], [sample_index, sample_index + 1])
        selected = segments[left:right]
        for chrom, (indices, positions) in marker_axes.items():
            for hap in (0, 1):
                rows = selected[(selected[:, 1] == chrom) & (selected[:, 2] == hap)]
                rows = rows[np.argsort(rows[:, 3])]
                cursor = spans[str(chrom)][0]
                for _, _, _, start, end, donor, donor_hap, source in rows:
                    if (start != cursor or end <= start or end > spans[str(chrom)][1]
                            or donor_hap not in (0, 1) or not 0 <= donor < len(labels)
                            or not 0 <= source < 7 or source_index[donor] != source
                            or labels.iloc[donor].sample_id not in plan["donors"]):
                        raise ValueError("Independent mosaic span/founder/source audit failed")
                    columns = indices[(positions >= start) & (positions < end)]
                    calls[hap, columns] = founders[donor, donor_hap, columns]
                    donor_choices[hap, columns] = donor
                    hap_choices[hap, columns] = donor_hap
                    assigned_bp[source] += end - start
                    cursor = end
                if cursor != spans[str(chrom)][1]:
                    raise ValueError("Independent mosaic audit found an uncovered chromosome copy")
        if not np.isin(calls, [0, 1]).all() or assigned_bp.sum() != total_bp:
            raise ValueError("Incomplete inherited mosaic genotypes or physical truth")
        np.testing.assert_array_equal(calls.sum(axis=0), x[sample_index])
        same = donor_choices[0] == donor_choices[1]
        if np.any(hap_choices[0, same] == hap_choices[1, same]):
            raise ValueError("Mosaic duplicated a founder haplotype across both copies")
        columns = np.flatnonzero(same)
        selected_donors = donor_choices[0, same]
        np.testing.assert_array_equal(
            x[sample_index, same],
            founders[selected_donors, 0, columns] + founders[selected_donors, 1, columns])
        np.testing.assert_allclose(assigned_bp / total_bp, truth[sample_index], atol=1e-12, rtol=0)
        target = np.zeros(7)
        for donor, source, fraction in zip(plan["donors"], plan["sources"], plan["target_weights"], strict=True):
            if labels.set_index("sample_id").loc[donor, "source_group"] != source:
                raise ValueError("Augmentation plan donor/source disagreement")
            target[SOURCES.index(source)] = fraction
        np.testing.assert_allclose(truth[sample_index], target, atol=2 / total_bp + 1e-12, rtol=0)
    splits = json.loads((out / "folds.json").read_text())
    strata = {plan["stratum"] for plan in plans}
    supports = []
    for split in splits:
        held = set(split["heldout_donors"])
        safe = [plan for plan in plans if not set(plan["donors"]) & held]
        counts = {stratum: sum(plan["stratum"] == stratum for plan in safe) for stratum in strata}
        if min(counts.values()) < settings["minimum_safe_per_stratum_per_fold"]:
            raise ValueError("Mosaic bank lacks the predeclared donor-safe stratum coverage")
        supports.append(dict(fold=split["fold_id"], n_safe=len(safe), minimum_per_stratum=min(counts.values())))
    minor = (truth >= config["gates"]["minor_min"]) & (truth <= config["gates"]["minor_max"])
    adapter.write_json(out / "augmentation_verification.json", dict(
        status="complete_verified_training_only_mosaics",
        experiment_sha256=adapter.digest(out / "experiment.json"),
        n_training_founders=307, n_mosaics=len(plans), n_strata=len(strata),
        n_inherited_genotypes_verified=int(x.size), physical_truth_independently_recomputed=True,
        same_founder_diploid_genotypes_preserved=True,
        all_donors_training_only=True, fold_support=supports,
        minor_contributions_by_source=dict(zip(SOURCES, minor.sum(axis=0).astype(int).tolist())),
        synthetic_people_in_validation=0, independent_final_test=False, production_promotion=False))
    print("Independent mosaic genotype/physical-truth/donor-coverage verification passed.", flush=True)


def verify(out):
    from sweep_phasefree_balanced import CALIBRATION, CONTROL, fold_models, predict

    config = configuration(out, full=True)
    if not all(checkpoint(out / name, out) for name in ("cv_completion.json", "development_completion.json")):
        raise ValueError("Require completed balanced CV and locked development results")
    metadata, features = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    plans = json.loads((out / "augmentation_plans.json").read_text())
    donors = {family: set(ids) for family, ids in json.loads((out / "donors.json").read_text()).items()}
    splits = json.loads((out / "folds.json").read_text())
    if splits != json.loads(Path("/rich/folds.json").read_text()):
        raise ValueError("Augmentation changed the established donor-purged validation splits")
    candidates = list(config["protocol"]["models"])
    reports, selection = (json.loads((out / name).read_text()) for name in ("cv_metrics.json", "selection.json"))
    with np.load(out / "oof_predictions.npz") as data:
        oof = data["predictions"]
        if data["sample_ids"].tolist() != metadata.sample_id.tolist() or data["candidates"].tolist() != candidates:
            raise ValueError("Balanced OOF sample/candidate axes changed")
    lookup, seen = {identity: i for i, identity in enumerate(metadata.sample_id)}, set()
    for split in splits:
        held, tested = set(split["heldout_donors"]), set(split["test_ids"])
        train = metadata[metadata.sample_id.isin(split["train_ids"])]
        used = set(train.loc[train.track.eq("pure"), "sample_id"])
        for family in train.loc[train.track.eq("pedigree"), "family_id"].unique():
            used.update(donors[family])
        if used & held or tested & set(train.sample_id) or seen & tested:
            raise ValueError("Independent original-training donor leakage/OOF overlap")
        seen.update(tested)
        directory = out / "folds" / split["fold_id"]
        if not checkpoint(directory / "completion.json", out):
            raise ValueError("Missing balanced fold checkpoint")
        record = json.loads((directory / "completion.json").read_text())
        expected = [i for i, plan in enumerate(plans) if not set(plan["donors"]) & held]
        if (record["fold"] != split or record["augmented_indices"] != expected
                or max(record["max_one_ulp_prediction_difference"].values()) > 1e-12):
            raise ValueError("Independent augmentation donor purge or numeric check failed")
        with np.load(directory / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if set(ids) != tested or data["candidates"].tolist() != candidates:
                raise ValueError("Balanced fold sample/candidate mismatch")
            np.testing.assert_array_equal(data["predictions"], oof[:, [lookup[identity] for identity in ids]])
    if seen != set(lookup) or len(seen) != 799 or len(splits) != 114:
        raise ValueError("Augmentation introduced validation people or incomplete original OOF coverage")
    with np.load(CALIBRATION / "oof_predictions.npz") as prior:
        np.testing.assert_allclose(oof[candidates.index(CONTROL)],
                                   prior["predictions"][prior["candidates"].tolist().index("stable_squared")],
                                   atol=1e-12, rtol=0)
    truth, tolerance = metadata[TRUTH].to_numpy(float), config["targets"]["comparison_tolerance"]
    for i, name in enumerate(candidates):
        np.testing.assert_allclose(
            independently_balanced_tv(metadata, truth, oof[i], config["gates"]["absent_ancestry_max"]),
            reports[name]["unified"]["balanced_mean_tv"], atol=1e-12, rtol=0)
        for track in ("pure", "pedigree"):
            mask = metadata.track.eq(track).to_numpy()
            check_metrics(truth[mask], oof[i, mask], reports[name][track], tolerance)
    order = sorted(candidates, key=lambda name: (
        reports[name]["unified"]["balanced_mean_tv"], reports[name]["pedigree"]["mean_tv"],
        -reports[name]["pedigree"]["fraction_all_components_within_5pp"],
        -reports[name]["pedigree"]["fraction_all_components_within_10pp"], name))
    if (selection["ranking"] != order or selection["development_metrics_used"]
            or selection["synthetic_people_in_validation"] != 0
            or selection["finalists"] != [name for name in order if name != CONTROL][:2]):
        raise ValueError("Balanced model selection violated training-only original-person ranking")
    development = json.loads((out / "development_completion.json").read_text())
    if development["selection_sha256"] != adapter.digest(out / "selection.json"):
        raise ValueError("Balanced selection changed after development evaluation")
    axes, all_x = read_table(out / "feature_axes.tsv"), np.load(out / "all_features.npy")
    chosen = axes.group.isin(["pure", "pedigree", "real"]).to_numpy()
    checks = []
    with threadpool_limits(limits=2):
        for split in (splits[0], next(row for row in splits if row["track"] == "pure")):
            models, _ = fold_models(out, split)
            test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
            for name in [CONTROL, *selection["finalists"]]:
                error = float(np.max(np.abs(predict(models[name], features[test]) - oof[candidates.index(name), test])))
                if error > 1e-12:
                    raise ValueError("Independent augmented-training refit predictions differ")
                checks.append(dict(fold=split["fold_id"], candidate=name, max_error=error))
        with (CALIBRATION / "candidates/stable_squared/model.pkl").open("rb") as stream:
            previous = pickle.load(stream)
        for name in [CONTROL, *selection["finalists"]]:
            directory = out / "candidates" / name
            with (directory / "model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            frame = read_table(directory / "predictions.tsv")
            if (frame.sample_id.tolist() != axes.loc[chosen, "sample_id"].tolist()
                    or frame.dataset.tolist() != axes.loc[chosen, "group"].tolist()
                    or model["sources"] != SOURCES
                    or model["experiment_sha256"] != adapter.digest(out / "experiment.json")):
                raise ValueError("Augmented saved-model metadata mismatch")
            prediction = predict(model["decoder"], all_x[chosen])
            np.testing.assert_allclose(frame[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                                       prediction, atol=1e-12, rtol=0)
            if name == CONTROL:
                np.testing.assert_allclose(prediction, calibration_predict(previous["decoder"], all_x[chosen, :14]),
                                           atol=1e-12, rtol=0)
            standard = json.loads((directory / "standard_validation.json").read_text())
            for track in ("pure", "pedigree"):
                scored = read_table(directory / f"{track}_predictions_scored.tsv")
                check_metrics(scored[TRUTH].to_numpy(float),
                              scored[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                              standard[track], tolerance)
            check_metrics(np.tile([.5, .5, 0, 0, 0, 0, 0], (5, 1)),
                          predict(model["decoder"], all_x[axes.group.eq("array_imputed").to_numpy()]),
                          standard["old_five_eur_asian_array_controls"], tolerance)
    adapter.write_json(out / "verification.json", dict(
        status="complete_verified_cv_and_development",
        cv_completion_sha256=adapter.digest(out / "cv_completion.json"),
        development_completion_sha256=adapter.digest(out / "development_completion.json"),
        augmentation_verification_sha256=adapter.digest(out / "augmentation_verification.json"),
        original_and_augmented_donor_overlaps_zero=True, n_original_oof_people=799,
        n_donor_purged_folds=114, synthetic_people_in_validation=0,
        original_stable_control_reproduced=True, metrics_independently_recomputed=True,
        refit_spot_checks=checks, expanded_array_verified=False,
        independent_final_test=False, production_promotion=False))
    print("Independent augmented donor/OOF/metric/refit/development verification passed.", flush=True)
