#!/usr/bin/env python3
"""Independent donor, metric, refit and transfer checks for richer-feature sweeps."""

import json
from pathlib import Path
import pickle

import numpy as np
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from phasefree_validation import SOURCES, TRUTH, fold_weights
from verify_phasefree_sweep import independently_balanced_tv
from sweep_phasefree_calibration import checkpoint, configuration, read_table


def check_metrics(truth, prediction, report, tolerance):
    error = np.abs(prediction - truth)
    if (prediction.shape != truth.shape or not len(truth) or not np.isfinite(prediction).all()
            or np.any(prediction < -tolerance)
            or not np.allclose(prediction.sum(axis=1), 1, atol=tolerance, rtol=0)):
        raise ValueError("Invalid richer ancestry predictions")
    np.testing.assert_allclose(
        [error.sum(axis=1).mean() / 2, (error.max(axis=1) <= .10 + tolerance).mean(),
         (error.max(axis=1) <= .05 + tolerance).mean()],
        [report["mean_tv"], report["fraction_all_components_within_10pp"],
         report["fraction_all_components_within_5pp"]], atol=1e-12, rtol=0)


def verify(out):
    from sweep_phasefree_rich import CALIBRATION, CONTROL, fit_models, predict

    config = configuration(out, full=True)
    if not all(checkpoint(out / name, out) for name in ("cv_completion.json", "development_completion.json")):
        raise ValueError("Require complete richer CV and locked development results")
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    donors = {family: set(ids) for family, ids in json.loads((out / "donors.json").read_text()).items()}
    splits = json.loads((out / "folds.json").read_text())
    if splits != json.loads((CALIBRATION / "folds.json").read_text()):
        raise ValueError("Richer models changed the established donor-purged splits")
    candidates = list(config["protocol"]["models"])
    reports = json.loads((out / "cv_metrics.json").read_text())
    selection = json.loads((out / "selection.json").read_text())
    with np.load(out / "oof_predictions.npz") as data:
        oof = data["predictions"]
        if data["sample_ids"].tolist() != metadata.sample_id.tolist() or data["candidates"].tolist() != candidates:
            raise ValueError("Richer OOF identity/candidate axis mismatch")
    index, seen = {identity: i for i, identity in enumerate(metadata.sample_id)}, set()
    for split in splits:
        train = metadata[metadata.sample_id.isin(split["train_ids"])]
        used = set(train.loc[train.track.eq("pure"), "sample_id"])
        for family in train.loc[train.track.eq("pedigree"), "family_id"].unique():
            used.update(donors[family])
        test = set(split["test_ids"])
        if used & set(split["heldout_donors"]) or set(split["train_ids"]) & test or seen & test:
            raise ValueError("Independent richer donor/test overlap audit failed")
        seen.update(test)
        directory = out / "folds" / split["fold_id"]
        if not checkpoint(directory / "completion.json", out):
            raise ValueError("Missing verified richer fold")
        record = json.loads((directory / "completion.json").read_text())
        if record["fold"] != split or max(record["max_one_ulp_prediction_difference"].values()) > 1e-12:
            raise ValueError("Richer fold identity or numeric contract changed")
        with np.load(directory / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if set(ids) != test or data["candidates"].tolist() != candidates:
                raise ValueError("Changed richer fold prediction axes")
            np.testing.assert_array_equal(data["predictions"], oof[:, [index[identity] for identity in ids]])
    if seen != set(index) or len(seen) != 799 or len(splits) != 114:
        raise ValueError("Richer OOF validation does not cover the exact established cohort")
    with np.load(CALIBRATION / "oof_predictions.npz") as prior:
        np.testing.assert_allclose(oof[candidates.index(CONTROL)],
                                   prior["predictions"][prior["candidates"].tolist().index("stable_squared")],
                                   atol=1e-12, rtol=0)
    tolerance = config["targets"]["comparison_tolerance"]
    truth = metadata[TRUTH].to_numpy(float)
    for i, name in enumerate(candidates):
        balanced = independently_balanced_tv(metadata, truth, oof[i], config["gates"]["absent_ancestry_max"])
        np.testing.assert_allclose(balanced, reports[name]["unified"]["balanced_mean_tv"], atol=1e-12, rtol=0)
        for track in ("pure", "pedigree"):
            chosen = metadata.track.eq(track).to_numpy()
            check_metrics(truth[chosen], oof[i, chosen], reports[name][track], tolerance)
    order = sorted(candidates, key=lambda name: (
        reports[name]["unified"]["balanced_mean_tv"], reports[name]["pedigree"]["mean_tv"],
        -reports[name]["pedigree"]["fraction_all_components_within_5pp"],
        -reports[name]["pedigree"]["fraction_all_components_within_10pp"], name))
    if (selection["ranking"] != order or selection["development_metrics_used"]
            or selection["finalists"] != [name for name in order if name != CONTROL][:2]):
        raise ValueError("Richer finalist selection does not match the predeclared training-only ranking")
    completion = json.loads((out / "development_completion.json").read_text())
    if completion["selection_sha256"] != adapter.digest(out / "selection.json"):
        raise ValueError("Richer selection changed after development evaluation")
    axes, values = read_table(out / "feature_axes.tsv"), np.load(out / "all_features.npy")
    evaluation = axes.group.isin(["pure", "pedigree", "real"]).to_numpy()
    checks = []
    with threadpool_limits(limits=2):
        for split in (splits[0], next(row for row in splits if row["track"] == "pure")):
            train = np.flatnonzero(metadata.sample_id.isin(split["train_ids"]))
            test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
            models = fit_models(x[train], truth[train], fold_weights(metadata.iloc[train]), config["protocol"])
            for name in [CONTROL, *selection["finalists"]]:
                error = float(np.max(np.abs(predict(models[name], x[test]) - oof[candidates.index(name), test])))
                if error > 1e-12:
                    raise ValueError("Independent richer CV refit predictions differ")
                checks.append(dict(fold=split["fold_id"], candidate=name, max_error=error))
        for name in [CONTROL, *selection["finalists"]]:
            directory = out / "candidates" / name
            with (directory / "model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            if model["sources"] != SOURCES or model["experiment_sha256"] != adapter.digest(out / "experiment.json"):
                raise ValueError("Richer saved-model provenance mismatch")
            frame = read_table(directory / "predictions.tsv")
            if (frame.sample_id.tolist() != axes.loc[evaluation, "sample_id"].tolist()
                    or frame.dataset.tolist() != axes.loc[evaluation, "group"].tolist()):
                raise ValueError("Richer development prediction identities changed")
            np.testing.assert_allclose(frame[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                                       predict(model["decoder"], values[evaluation]), atol=1e-12, rtol=0)
            standard = json.loads((directory / "standard_validation.json").read_text())
            for track in ("pure", "pedigree"):
                scored = read_table(directory / f"{track}_predictions_scored.tsv")
                check_metrics(scored[TRUTH].to_numpy(float),
                              scored[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                              standard[track], tolerance)
    adapter.write_json(out / "verification.json", dict(
        status="complete_verified_cv_and_development", n_donor_purged_folds=114, n_oof_people=799,
        all_fold_donor_overlaps_zero=True, prior_stable_control_oof_reproduced=True,
        metrics_independently_recomputed=True, refit_spot_checks=checks,
        cv_completion_sha256=adapter.digest(out / "cv_completion.json"),
        development_completion_sha256=adapter.digest(out / "development_completion.json"),
        expanded_array_verified=False, current_dataset_only=True,
        independent_final_test=False, production_promotion=False))
    print("Independent richer CV donor/metric/refit and development checks passed.", flush=True)


def verify_arrays(out, arrays):
    from sweep_phasefree_calibration import predict as calibration_predict
    from sweep_phasefree_rich import CALIBRATION, CONTROL, predict

    config = configuration(out, full=True)
    if not checkpoint(out / "array_completion.json", out):
        raise ValueError("Richer shared-array comparison has not completed")
    record = json.loads((out / "array_completion.json").read_text())
    if (record["source_verification_sha256"] != adapter.digest(arrays / "verification.json")
            or record["selection_sha256"] != adapter.digest(out / "selection.json")):
        raise ValueError("Shared-array source or locked richer selection changed")
    selection = json.loads((out / "selection.json").read_text())
    plan = json.loads((arrays / "experiment.json").read_text())
    names = [row["sample_id"] for row in plan["samples"]]
    frame, comparison = read_table(out / "array_predictions.tsv"), read_table(out / "array_comparison.tsv")
    candidates = [CONTROL, *selection["finalists"]]
    if (len(frame) != len(candidates) * 28 * 2 or set(frame.candidate) != set(candidates)
            or frame.duplicated(["candidate", "sample_id", "mask", "stage"]).any()):
        raise ValueError("Incomplete richer shared-array conditions")
    with np.load(out / "array_features.npz") as data:
        clean, imputed = data["clean"], data["imputed"]
        if (data["sample_ids"].tolist() != names
                or data["case_ids"].tolist() != [case["case_id"] for case in plan["cases"]]):
            raise ValueError("Richer shared-array feature identities changed")
    with threadpool_limits(limits=2):
        with (CALIBRATION / "candidates/stable_squared/model.pkl").open("rb") as stream:
            previous = pickle.load(stream)
        for name in candidates:
            with (out / f"candidates/{name}/model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            original, recovered = predict(model["decoder"], clean), predict(model["decoder"], imputed)
            if name == CONTROL:
                np.testing.assert_allclose(original, calibration_predict(previous["decoder"], clean[:, :14]),
                                           atol=1e-12, rtol=0)
                np.testing.assert_allclose(recovered, calibration_predict(previous["decoder"], imputed[:, :14]),
                                           atol=1e-12, rtol=0)
            for i, case in enumerate(plan["cases"]):
                for stage, prediction in (("clean", original[names.index(case["sample_id"])]), ("array_imputed", recovered[i])):
                    row = frame[(frame.candidate == name) & (frame["mask"] == case["mask"])
                                & (frame.sample_id == case["sample_id"]) & (frame.stage == stage)]
                    if len(row) != 1 or row.iloc[0].track != case["track"]:
                        raise ValueError("Richer shared-array prediction metadata changed")
                    np.testing.assert_allclose(row[TRUTH].to_numpy(float)[0], case["truth"], atol=1e-12, rtol=0)
                    np.testing.assert_allclose(row[SOURCES].to_numpy(float)[0], prediction, atol=1e-12, rtol=0)
    for row in comparison.to_dict("records"):
        chosen = frame[(frame.candidate == row["candidate"]) & (frame["mask"] == row["mask"])
                       & (frame.track == row["track"]) & (frame.stage == row["stage"])]
        if len(chosen) != row["n_people"]:
            raise ValueError("Richer array metric support mismatch")
        check_metrics(chosen[TRUTH].to_numpy(float), chosen[SOURCES].to_numpy(float),
                      row, config["targets"]["comparison_tolerance"])
    reports, decisions = json.loads((out / "cv_metrics.json").read_text()), []
    for name in selection["finalists"]:
        directory = out / "candidates" / name
        metrics = json.loads((directory / "metrics.json").read_text())
        standard = json.loads((directory / "standard_validation.json").read_text())
        group = comparison[comparison.candidate.eq(name) & comparison.stage.eq("array_imputed")]
        if len(group) != 4:
            raise ValueError("Require both masks and both pure/mixed tracks")
        checks = dict(
            improves_training_cv=reports[name]["unified"]["balanced_mean_tv"]
            < reports[CONTROL]["unified"]["balanced_mean_tv"] - 1e-12,
            unchanged_development_gates=bool(metrics["all_gates_pass"]),
            original_array5_diagnostic=standard["old_five_eur_asian_array_controls"]["mean_tv"] <= .10,
            expanded_array_diagnostics=bool((group.mean_tv <= config["protocol"]["array_mean_tv_diagnostic_max"]).all()),
            user_pedigree_mean_acceptable=standard["pedigree"]["mean_tv"] <= .10)
        decisions.append(dict(candidate=name, checks=checks, eligible=all(checks.values())))
    eligible = [row["candidate"] for row in decisions if row["eligible"]]
    adapter.write_json(out / "array_verification.json", dict(
        status="complete_verified_shared_array_comparison", metrics_independently_recomputed=True,
        saved_model_predictions_reproduced=True, n_unique_people=14, n_conditions=28,
        array_completion_sha256=adapter.digest(out / "array_completion.json"),
        independent_final_test=False, production_promotion=False))
    adapter.write_json(out / "recommendation.json", dict(
        status="eligible_for_further_validation" if eligible else "retain_frozen_baseline",
        candidate=eligible[0] if eligible else "phasefree_stack_boost_clean", decisions=decisions,
        selection_used_training_cv_order_with_predeclared_guardrail_filtering=True,
        current_dataset_only=True, independent_final_test=False, production_promotion=False,
        pgs_calibration_changed=False))
