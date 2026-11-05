#!/usr/bin/env python3
"""Verify donor separation, out-of-fold ranking and frozen-model sweep predictions."""

import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter


def independently_balanced_tv(metadata, truth, prediction, absent):
    errors = np.abs(prediction - truth).sum(axis=1) / 2
    combinations = [tuple(np.flatnonzero(row >= absent)) for row in truth]
    scenario_means = []
    for scenario in sorted(set(metadata.scenario)):
        indices = np.flatnonzero(metadata.scenario.eq(scenario).to_numpy())
        combinations_here = sorted({combinations[index] for index in indices})
        combination_means = []
        for combination in combinations_here:
            chosen = [index for index in indices if combinations[index] == combination]
            family_means = [
                float(np.mean([errors[index] for index in chosen
                               if metadata.iloc[index].family_id == family]))
                for family in sorted(set(metadata.iloc[chosen].family_id))
            ]
            combination_means.append(float(np.mean(family_means)))
        scenario_means.append(float(np.mean(combination_means)))
    return float(np.mean(scenario_means))


def verify(out):
    from phasefree_validation import FEATURES, SOURCES, TRUTH, fold_weights
    from sweep_phasefree_calibration import checkpoint, configuration, fit_models, predict, read_table

    config = configuration(out, full=True)
    if not checkpoint(out / "cv_completion.json", out) or not checkpoint(out / "development_completion.json", out):
        raise ValueError("CV and locked-finalist development stages must complete before verification")
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    donors = {family: set(values) for family, values in json.loads((out / "donors.json").read_text()).items()}
    splits = json.loads((out / "folds.json").read_text())
    candidates = list(config["protocol"]["models"])
    cv_reports = json.loads((out / "cv_metrics.json").read_text())
    index = {identity: position for position, identity in enumerate(metadata.sample_id)}
    seen = set()
    with np.load(out / "oof_predictions.npz") as data:
        oof = data["predictions"]
        if data["sample_ids"].tolist() != metadata.sample_id.tolist() or data["candidates"].tolist() != candidates:
            raise ValueError("OOF verification candidate/sample axis mismatch")
    for split in splits:
        train_ids, test_ids = set(split["train_ids"]), set(split["test_ids"])
        held = set(split["heldout_donors"])
        train = metadata[metadata.sample_id.isin(train_ids)]
        training_donors = set(train.loc[train.track.eq("pure"), "sample_id"])
        for family in train.loc[train.track.eq("pedigree"), "family_id"]:
            training_donors.update(donors[family])
        if training_donors & held or train_ids & test_ids or seen & test_ids:
            raise ValueError(f"Independent donor/OOF overlap audit failed: {split['fold_id']}")
        seen.update(test_ids)
        path = out / "folds" / split["fold_id"]
        if not checkpoint(path / "completion.json", out):
            raise ValueError("Missing verified CV fold")
        recorded = json.loads((path / "completion.json").read_text())
        if recorded["fold"] != split:
            raise ValueError("Fold checkpoint uses different split identities")
        with np.load(path / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if set(ids) != test_ids:
                raise ValueError("Fold prediction identities differ from held-out identities")
            np.testing.assert_array_equal(data["predictions"], oof[:, [index[identity] for identity in ids]])
    if seen != set(index):
        raise ValueError("Independent OOF verification does not cover every training person exactly once")
    tolerance = config["targets"]["comparison_tolerance"]
    truth = metadata[TRUTH].to_numpy(float)
    for position, name in enumerate(candidates):
        prediction = oof[position]
        if (prediction.shape != truth.shape or not np.isfinite(prediction).all()
                or np.any(prediction < -tolerance)
                or not np.allclose(prediction.sum(axis=1), 1, atol=tolerance, rtol=0)):
            raise ValueError("Invalid out-of-fold ancestry fractions")
        expected_balanced = independently_balanced_tv(
            metadata, truth, prediction, config["gates"]["absent_ancestry_max"])
        if abs(expected_balanced - cv_reports[name]["unified"]["balanced_mean_tv"]) > 1e-12:
            raise ValueError("Independent scenario/combination/family macro TV mismatch")
        for track in ("pure", "pedigree"):
            rows = metadata.track.eq(track).to_numpy()
            error = np.abs(prediction[rows] - truth[rows])
            saved = cv_reports[name][track]
            numbers = (float(error.sum(axis=1).mean() / 2), float((error.max(axis=1) <= .10 + tolerance).mean()),
                       float((error.max(axis=1) <= .05 + tolerance).mean()))
            np.testing.assert_allclose(numbers, [
                saved["mean_tv"], saved["fraction_all_components_within_10pp"],
                saved["fraction_all_components_within_5pp"]], atol=1e-12, rtol=0)
    expected_order = sorted(candidates, key=lambda name: (
        cv_reports[name]["unified"]["balanced_mean_tv"], cv_reports[name]["pedigree"]["mean_tv"],
        -cv_reports[name]["pedigree"]["fraction_all_components_within_5pp"],
        -cv_reports[name]["pedigree"]["fraction_all_components_within_10pp"], name))
    selection = json.loads((out / "selection.json").read_text())
    if (selection["ranking"] != expected_order
            or selection["finalists"] != [name for name in expected_order if name != "refit_control"][:2]
            or selection["development_metrics_used"]):
        raise ValueError("Training-only finalist selection verification failed")
    development = json.loads((out / "development_completion.json").read_text())
    if development["selection_sha256"] != adapter.digest(out / "selection.json"):
        raise ValueError("CV selection changed after development evaluation")
    with Path(config["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    features = read_table(Path("/stack/features.tsv"))
    maximum_error, refit_checks = 0., []
    with threadpool_limits(limits=2):
        for split in (splits[0], next(row for row in splits if row["track"] == "pure")):
            train = np.flatnonzero(metadata.sample_id.isin(split["train_ids"]))
            test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
            refitted = fit_models(x[train], metadata.iloc[train][TRUTH].to_numpy(float),
                                  fold_weights(metadata.iloc[train]), config["protocol"])
            for name in ["refit_control", *selection["finalists"]]:
                expected = predict(refitted[name], x[test])
                saved = oof[candidates.index(name), test]
                error = float(np.max(np.abs(expected - saved)))
                maximum_error = max(maximum_error, error)
                if error > 1e-12:
                    raise ValueError("Separately refitted CV spot-check predictions differ")
                refit_checks.append(dict(fold=split["fold_id"], candidate=name, max_error=error))
        for name in ["refit_control", *selection["finalists"]]:
            directory = out / "candidates" / name
            with (directory / "model.pkl").open("rb") as stream:
                bundle = pickle.load(stream)
            if (bundle["sources"] != SOURCES or bundle["marker_keys"] != original["marker_keys"]
                    or bundle["experiment_sha256"] != adapter.digest(out / "experiment.json")):
                raise ValueError("Full-model bundle has incompatible source/marker/provenance contract")
            saved = read_table(directory / "predictions.tsv")
            selected_features = features[features.group.isin(["pure", "pedigree", "real"])]
            joined = saved.merge(selected_features, left_on=["dataset", "sample_id"],
                                 right_on=["group", "sample_id"], validate="one_to_one", how="outer",
                                 indicator=True)
            if not joined._merge.eq("both").all():
                raise ValueError("Development prediction/feature coverage differs")
            reproduced = predict(bundle["decoder"], joined[FEATURES].to_numpy(float))
            expected = joined[[f"pred_{source}" for source in SOURCES]].to_numpy(float)
            error = float(np.max(np.abs(expected - reproduced)))
            maximum_error = max(maximum_error, error)
            if error > 1e-12:
                raise ValueError("Saved full model does not reproduce held-out predictions")
            for track in ("pure", "pedigree"):
                scored = read_table(directory / f"{track}_predictions_scored.tsv")
                tv = np.abs(scored[[f"pred_{source}" for source in SOURCES]].to_numpy(float)
                            - scored[TRUTH].to_numpy(float)).sum(axis=1) / 2
                metrics = json.loads((directory / "metrics.json").read_text())
                if abs(float(tv.mean()) - metrics[track]["mean_tv"]) > 1e-12:
                    raise ValueError("Independent development mean TV differs")
    if adapter.digest(config["original_model_path"]) != config["original_model_sha256"]:
        raise ValueError("Original frozen model changed")
    adapter.write_json(out / "verification.json", dict(
        status="complete_verified_cv_and_development", n_donor_purged_folds=len(splits),
        n_oof_people=len(metadata), n_candidates=len(candidates),
        all_fold_donor_overlaps_zero=True, oof_each_person_exactly_once=True,
        selection_locked_without_development_metrics=True,
        balanced_macro_and_accuracy_bands_independently_recomputed=True,
        cv_refit_spot_checks=refit_checks, max_prediction_reproduction_error=maximum_error,
        experiment_sha256=adapter.digest(out / "experiment.json"),
        cv_completion_sha256=adapter.digest(out / "cv_completion.json"),
        development_completion_sha256=adapter.digest(out / "development_completion.json"),
        expanded_array_verified=False, expanded_array_status="Requires separate verified array comparison.",
        independent_final_test=False, production_promotion=False,
        verifier_sha256=adapter.digest(__file__)))
    print("Independent CV donor/metric/refit and development prediction verification passed.", flush=True)


def verify_arrays(out, arrays):
    from phasefree_validation import FEATURES, SOURCES, TRUTH
    from sweep_phasefree_calibration import checkpoint, configuration, predict, read_table

    config = configuration(out, full=True)
    if not checkpoint(out / "array_completion.json", out):
        raise ValueError("No completed standardized array comparison to verify")
    record = json.loads((out / "array_completion.json").read_text())
    if (record["source_verification_sha256"] != adapter.digest(arrays / "verification.json")
            or record["source_experiment_sha256"] != adapter.digest(arrays / "experiment.json")
            or record["selection_sha256"] != adapter.digest(out / "selection.json")):
        raise ValueError("Array source or training-only finalist selection changed")
    plan = json.loads((arrays / "experiment.json").read_text())
    selection = json.loads((out / "selection.json").read_text())
    predictions = read_table(out / "array_predictions.tsv")
    comparison = read_table(out / "array_comparison.tsv")
    candidates = ["refit_control", *selection["finalists"]]
    if (len(predictions) != len(candidates) * 28 * 2
            or predictions.duplicated(["candidate", "mask", "sample_id", "stage"]).any()
            or set(predictions.candidate) != set(candidates)):
        raise ValueError("Standardized array candidate/condition coverage mismatch")
    with np.load(out / "array_features.npz") as data:
        clean, imputed = data["clean"], data["imputed"]
        if (data["sample_ids"].tolist() != [row["sample_id"] for row in plan["samples"]]
                or data["case_ids"].tolist() != [row["case_id"] for row in plan["cases"]]
                or clean.shape != (14, len(FEATURES)) or imputed.shape != (28, len(FEATURES))):
            raise ValueError("Standardized array feature/identity axis mismatch")
    names = [row["sample_id"] for row in plan["samples"]]
    maximum_error = 0.
    with threadpool_limits(limits=2):
        for name in candidates:
            with (out / f"candidates/{name}/model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            original, recovered = predict(model["decoder"], clean), predict(model["decoder"], imputed)
            for index, case in enumerate(plan["cases"]):
                for stage, values in (
                    ("clean", original[names.index(case["sample_id"])]), ("array_imputed", recovered[index])
                ):
                    row = predictions[(predictions.candidate == name) & (predictions["mask"] == case["mask"])
                                      & (predictions.sample_id == case["sample_id"]) & (predictions.stage == stage)]
                    if len(row) != 1 or row.iloc[0]["track"] != case["track"]:
                        raise ValueError("Array prediction metadata mismatch")
                    np.testing.assert_allclose(row[TRUTH].to_numpy(float)[0], case["truth"], atol=1e-12, rtol=0)
                    error = float(np.max(np.abs(row[SOURCES].to_numpy(float)[0] - values)))
                    maximum_error = max(maximum_error, error)
                    if error > 1e-12:
                        raise ValueError("Saved-model array predictions do not reproduce")
    for row in comparison.to_dict("records"):
        group = predictions[(predictions.candidate == row["candidate"]) & (predictions["mask"] == row["mask"])
                            & (predictions.track == row["track"]) & (predictions.stage == row["stage"])]
        error = np.abs(group[SOURCES].to_numpy(float) - group[TRUTH].to_numpy(float))
        if len(group) != row["n_people"]:
            raise ValueError("Array metric support mismatch")
        np.testing.assert_allclose(
            [float(error.sum(axis=1).mean() / 2), float((error.max(axis=1) <= .10 + 1e-12).mean()),
             float((error.max(axis=1) <= .05 + 1e-12).mean())],
            [row["mean_tv"], row["fraction_all_components_within_10pp"],
             row["fraction_all_components_within_5pp"]], atol=1e-12, rtol=0)
    prior = read_table(arrays / "predictions.tsv")
    control = predictions[predictions.candidate.eq("refit_control")]
    joined = control.merge(prior, on=["mask", "sample_id", "stage", "track"],
                           suffixes=("_new", "_frozen"), how="outer", validate="one_to_one", indicator=True)
    if not joined._merge.eq("both").all():
        raise ValueError("Control does not cover the same shared array inputs")
    np.testing.assert_allclose(joined[[f"{source}_new" for source in SOURCES]].to_numpy(float),
                               joined[[f"{source}_frozen" for source in SOURCES]].to_numpy(float),
                               atol=1e-12, rtol=0)
    cv = json.loads((out / "cv_metrics.json").read_text())
    decisions = []
    for name in selection["finalists"]:
        metrics = json.loads((out / f"candidates/{name}/metrics.json").read_text())
        standard = json.loads((out / f"candidates/{name}/standard_validation.json").read_text())
        array_groups = comparison[comparison.candidate.eq(name) & comparison.stage.eq("array_imputed")]
        if len(array_groups) != 4:
            raise ValueError("Expected both pure/mixed tracks under both masks")
        checks = dict(
            improves_training_cv=cv[name]["unified"]["balanced_mean_tv"]
            < cv["refit_control"]["unified"]["balanced_mean_tv"] - 1e-12,
            unchanged_development_gates=bool(metrics["all_gates_pass"]),
            original_array5_diagnostic=standard["old_five_eur_asian_array_controls"]["mean_tv"] <= .10,
            expanded_array_diagnostics=bool((array_groups.mean_tv <= .10).all()),
            user_pedigree_mean_acceptable=standard["pedigree"]["mean_tv"] <= .10)
        decisions.append(dict(candidate=name, checks=checks, eligible=all(checks.values())))
    eligible = [row["candidate"] for row in decisions if row["eligible"]]
    adapter.write_json(out / "array_verification.json", dict(
        status="complete_verified_shared_array_comparison",
        array_completion_sha256=adapter.digest(out / "array_completion.json"),
        source_verification_sha256=adapter.digest(arrays / "verification.json"),
        shared_frozen_control_predictions_reproduced=True,
        serialized_candidate_predictions_reproduced=True,
        metrics_independently_recomputed=True, max_prediction_error=maximum_error,
        n_unique_people=14, n_conditions=28,
        independent_final_test=False, production_promotion=False))
    adapter.write_json(out / "recommendation.json", dict(
        status="eligible_for_further_validation" if eligible else "retain_frozen_baseline",
        candidate=eligible[0] if eligible else "phasefree_stack_boost_clean",
        training_cv_ranking=selection["ranking"], decisions=decisions,
        expanded_array_verification_sha256=adapter.digest(out / "array_verification.json"),
        selection_used_training_cv_order_with_predeclared_guardrail_filtering=True,
        no_minimum_coverage_fraction_invented=True,
        independent_final_test=False, production_promotion=False, pgs_calibration_changed=False))
    print("Shared-array saved-model/metric/control verification and guarded recommendation complete.", flush=True)
