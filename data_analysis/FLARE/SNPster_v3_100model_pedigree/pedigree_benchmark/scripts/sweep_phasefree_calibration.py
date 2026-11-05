#!/usr/bin/env python3
"""Training-only donor-purged calibration sweep over frozen reference-only features."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import PolynomialFeatures
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from phasefree_validation import (
    FEATURES, SOURCES, TRUTH, donor_purged_splits, fold_weights, ranking_key,
    stable_features, summarize,
)
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import score
from screen_phasefree_stack import features as genotype_features
from screen_v1_decoders import family_donors, predict_decoder, simplex_rows


PROJECT = Path(__file__).resolve().parents[2]
STACK = Path("/stack")
MODEL_SHA = "95f3e3c9abaf76fb357d4ef51b2e29bf212a033002c02a04b825ec01b0231686"
PROTOCOL = PROJECT / "config/phasefree_sweep_protocol.json"
TARGETS = PROJECT / "config/phasefree_accuracy_targets.json"
GATES = PROJECT / "config/unified_selection.json"


def read_table(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def hgb(parameters, protocol):
    return MultiOutputRegressor(HistGradientBoostingRegressor(
        loss=parameters["loss"], max_iter=protocol["boosting"]["iterations"],
        learning_rate=protocol["boosting"]["learning_rate"],
        max_leaf_nodes=parameters["leaves"], min_samples_leaf=parameters["min_leaf"],
        l2_regularization=parameters["l2"], early_stopping=False,
        random_state=protocol["seed"]), n_jobs=1)


def fit_models(x, y, weights, protocol):
    stable = stable_features(x, protocol["feature_round_decimals"])
    if (y.shape != (len(x), len(SOURCES)) or weights.shape != (len(x),)
            or np.any(weights <= 0) or not np.isclose(weights.sum(), 1)):
        raise ValueError("Invalid calibration training axes/weights")
    linear = Ridge(alpha=.001, fit_intercept=True).fit(stable, y, sample_weight=weights)
    shared_head = hgb(protocol["models"]["stable_squared"], protocol).fit(stable, y, sample_weight=weights)
    models = {}
    for name, parameters in protocol["models"].items():
        values = stable if parameters["round_features"] else x
        model = dict(parameters=parameters, round_decimals=protocol["feature_round_decimals"])
        if parameters["kind"] == "hgb":
            model["head"] = (shared_head if name == "stable_squared" else
                             hgb(parameters, protocol).fit(values, y, sample_weight=weights))
        elif parameters["kind"] == "residual":
            model["head"] = linear
            model["correction"] = hgb(parameters, protocol).fit(
                stable, y - linear.predict(stable), sample_weight=weights)
        elif parameters["kind"] == "quadratic":
            transform = PolynomialFeatures(degree=2, include_bias=False)
            model["transform"] = transform
            model["head"] = Ridge(alpha=parameters["alpha"], fit_intercept=True).fit(
                transform.fit_transform(stable), y, sample_weight=weights)
        elif parameters["kind"] == "blend":
            model["head"] = shared_head
            model["linear"] = linear
        else:
            raise ValueError(f"Unknown fixed candidate: {parameters}")
        models[name] = model
    return models


def predict(model, x):
    parameters = model["parameters"]
    values = (stable_features(x, model["round_decimals"])
              if parameters["round_features"] else np.asarray(x, float))
    kind = parameters["kind"]
    if kind == "quadratic":
        values = model["transform"].transform(values)
    result = model["head"].predict(values)
    if kind == "residual":
        result += model["correction"].predict(values)
    if kind == "blend":
        weight = parameters["linear_fraction"]
        return (1 - weight) * simplex_rows(result) + weight * simplex_rows(model["linear"].predict(values))
    return simplex_rows(result)


def configuration(out, full=False):
    config = json.loads((out / "experiment.json").read_text())
    provenance.verify_hashes(config["implementation_sha256"])
    provenance.verify_hashes(config["prepared_output_sha256"])
    if full:
        provenance.verify_hashes(config["input_sha256"])
    return config


def checkpoint(path, out):
    if not path.exists():
        return False
    record = json.loads(path.read_text())
    if (record["status"] != "complete"
            or record["experiment_sha256"] != adapter.digest(out / "experiment.json")):
        raise ValueError(f"Incompatible sweep checkpoint: {path}")
    provenance.verify_hashes(record["output_sha256"])
    return True


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Prepare a fresh calibration sweep")
    out.mkdir(parents=True, exist_ok=True)
    protocol, targets = json.loads(PROTOCOL.read_text()), json.loads(TARGETS.read_text())
    if protocol["sources"] != SOURCES or targets["sources"] != SOURCES:
        raise ValueError("Source order differs from frozen feature models")
    selected = json.loads((STACK / "selected_candidate.json").read_text())
    if adapter.digest(STACK / "experiment.json") != selected["experiment_sha256"]:
        raise ValueError("Original frozen scientific configuration changed")
    model_path = STACK / selected["model_path"]
    if adapter.digest(model_path) != MODEL_SHA or selected["model_sha256"] != MODEL_SHA:
        raise ValueError("Frozen reference candidate changed")
    verified_path = STACK / "verification/final_verification.json"
    if adapter.digest(verified_path) != selected["verification_sha256"]:
        raise ValueError("Independent frozen-candidate verification changed")
    verified = json.loads(verified_path.read_text())
    if adapter.digest(STACK / "completion.json") != verified["completion_sha256"]:
        raise ValueError("Verified original completion changed")
    completed = json.loads((STACK / "completion.json").read_text())
    paths = [STACK / name for name in ("features.tsv", "training.tsv", "feature_models.pkl")]
    for path in paths:
        if adapter.digest(path) != completed["output_sha256"][str(path)]:
            raise ValueError(f"Verified calibration input changed: {path}")
    inputs = [*paths, model_path, STACK / "selected_candidate.json", verified_path,
              STACK / "completion.json", STACK / "cohort.json", STACK / "baseline/metrics.json",
              STACK / "experiment.json"]
    cohort = json.loads((STACK / "cohort.json").read_text())
    pure = read_table(cohort["pure_train"]["labels_tsv"])
    reference = read_table(cohort["reference"]["labels_tsv"])
    inputs += [Path(cohort[group]["labels_tsv"]) for group in ("pure_train", "reference")]
    pure_annotations = pd.DataFrame(dict(
        sample_id=pure.sample_id, source_group=pure.source_group, track="pure",
        family_id="pure:" + pure.sample_id, scenario="reference_source"))
    pedigree = pd.concat([
        read_table(batch["samples_tsv"]).assign(track="pedigree", source_group="")
        for batch in cohort["pedigree_train"]["batches"]], ignore_index=True)
    inputs += [Path(batch["samples_tsv"]) for batch in cohort["pedigree_train"]["batches"]]
    plan_path = PROJECT / (
        "results_v3_pedigree_training_screen_01/cohort/pedigree/plans/plans/expected_ancestry.tsv")
    inputs.append(plan_path)
    plan = read_table(plan_path)
    donors = family_donors(plan)
    founders = plan[plan.is_founder.astype(str).str.lower().isin(["true", "1"])]
    labels = pure.set_index("sample_id").source_group.to_dict()
    for row in founders.itertuples():
        source = "EAS_SEA" if row.founder_source in ("EAS", "SEA") else row.founder_source
        if labels.get(row.founder_donor) != source:
            raise ValueError("Training founder source disagrees with pure-training labels")
    annotations = pd.concat([
        pure_annotations, pedigree[["sample_id", "source_group", "track", "family_id", "scenario"]]
    ], ignore_index=True)
    all_features = read_table(STACK / "features.tsv")
    training = read_table(STACK / "training.tsv")
    clean_features = all_features[all_features.group.eq("clean_training")]
    clean = training[training.weight_clean > 0][["sample_id", *TRUTH, "weight_clean"]]
    metadata = clean_features[["sample_id"]].merge(clean, on="sample_id", validate="one_to_one")
    metadata = metadata.merge(annotations, on="sample_id", validate="one_to_one")
    forbidden = set(reference.sample_id) | set(all_features.loc[
        all_features.group.isin(["pure", "pedigree", "real", "array_imputed"]), "sample_id"])
    if (len(metadata) != 799 or len(pure) != 307 or len(pedigree) != 492
            or set(metadata.sample_id) != set(annotations.sample_id)
            or set(metadata.sample_id) & forbidden
            or set().union(*donors.values()) & forbidden):
        raise ValueError("Training-only sample/founder/reference/development separation failed")
    x = clean_features.set_index("sample_id").loc[metadata.sample_id, FEATURES].to_numpy(float)
    stable_features(x)
    y = metadata[TRUTH].to_numpy(float)
    if not np.allclose(y.sum(axis=1), 1, atol=1e-12, rtol=0) or np.any(y < 0):
        raise ValueError("Invalid frozen clean ancestry truth")
    if not np.isclose(metadata.weight_clean.sum(), 1):
        raise ValueError("Frozen clean weights do not sum to one")
    splits = donor_purged_splits(metadata, donors, protocol["pure_cv_folds"], protocol["split_seed"])
    metadata.to_csv(out / "training.tsv", sep="\t", index=False)
    np.save(out / "training_features.npy", x)
    adapter.write_json(out / "donors.json", {family: sorted(values) for family, values in donors.items()})
    adapter.write_json(out / "folds.json", splits)
    shutil.copyfile(STACK / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(STACK / "baseline/metrics.json", out / "baseline/metrics.json")
    implementation = [
        Path(__file__), Path(__file__).with_name("phasefree_validation.py"),
        Path(__file__).with_name("verify_phasefree_sweep.py"),
        Path(adapter.__file__), Path(provenance.__file__),
        Path(__file__).with_name("report_phasefree_accuracy.py"),
        Path(__file__).with_name("screen_v1_decoders.py"),
        Path(__file__).with_name("screen_model_families.py"),
        Path(__file__).with_name("screen_phasefree_stack.py"),
        Path(__file__).with_name("screen_dosage_models.py"),
        Path(__file__).with_name("screen_fine_profiles.py"),
        Path(__file__).with_name("screen_reference_trees.py"),
        Path(__file__).with_name("screen_supervised_genotypes.py"),
        Path(__file__).with_name("train_matched_global.py"),
        Path(__file__).with_name("evaluate_finalist.py"),
        PROJECT / "scripts/fit_evaluate_candidate.py", PROTOCOL, TARGETS, GATES]
    experiment = dict(
        status="prepared", protocol=protocol, targets=targets,
        gates=json.loads(GATES.read_text()), original_model_path=str(model_path),
        original_model_sha256=MODEL_SHA, no_evaluation_in_fit=True,
        no_imputed_training=True, no_base_classifier_or_frequency_refitting=True,
        n_training=799, n_pure=307, n_pedigree=492, n_cv_folds=len(splits),
        input_sha256=provenance.fingerprint_files(inputs),
        implementation_sha256=provenance.fingerprint_files(implementation),
        prepared_output_sha256=provenance.fingerprint_files([
            out / "training.tsv", out / "training_features.npy", out / "donors.json",
            out / "folds.json", out / "cohort.json", out / "baseline/metrics.json"]),
        independent_final_test=False, production_promotion=False)
    adapter.write_json(out / "experiment.json", experiment)
    print(f"Frozen {len(protocol['models'])} candidates, {len(splits)} donor-purged folds, "
          "799 clean training people; no evaluation fitting.", flush=True)


def run_fold(out, fold_id):
    config = configuration(out)
    splits = json.loads((out / "folds.json").read_text())
    split = next(row for row in splits if row["fold_id"] == fold_id)
    directory = out / "folds" / fold_id
    if checkpoint(directory / "completion.json", out):
        return fold_id
    if directory.exists():
        raise ValueError(f"Incomplete CV fold requires inspection: {directory}")
    directory.mkdir(parents=True)
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    train = np.flatnonzero(metadata.sample_id.isin(split["train_ids"]))
    test = np.flatnonzero(metadata.sample_id.isin(split["test_ids"]))
    weights = fold_weights(metadata.iloc[train])
    started = time.monotonic()
    predictions, stability = [], {}
    with threadpool_limits(limits=2):
        models = fit_models(x[train], metadata.iloc[train][TRUTH].to_numpy(float), weights, config["protocol"])
        for name, model in models.items():
            values = predict(model, x[test])
            restored = predict(pickle.loads(pickle.dumps(model)), x[test])
            if not np.array_equal(values, restored):
                raise ValueError(f"Serialized CV prediction mismatch: {name}, {fold_id}")
            difference = max(float(np.max(np.abs(values - predict(
                model, np.nextafter(x[test], direction))))) for direction in (-np.inf, np.inf))
            if model["parameters"]["round_features"] and difference > 1e-12:
                raise ValueError(f"Stable candidate is sensitive to one-ULP features: {name}, {fold_id}")
            stability[name] = difference
            predictions.append(values)
    np.savez_compressed(directory / "predictions.npz", predictions=np.asarray(predictions),
                        sample_ids=metadata.iloc[test].sample_id.to_numpy(str),
                        candidates=np.asarray(list(config["protocol"]["models"])))
    adapter.write_json(directory / "completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        fold=split, elapsed_seconds=time.monotonic() - started,
        max_one_ulp_prediction_difference=stability,
        output_sha256=provenance.fingerprint_files([directory / "predictions.npz"])))
    return fold_id


def cross_validate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "cv_completion.json", out):
        return
    splits = json.loads((out / "folds.json").read_text())
    with ProcessPoolExecutor(max_workers=2) as pool:
        tasks = {pool.submit(run_fold, out, split["fold_id"]): split["fold_id"] for split in splits}
        for index, future in enumerate(as_completed(tasks), 1):
            print(f"CV {index}/{len(splits)}: {future.result()} verified", flush=True)
    metadata = read_table(out / "training.tsv")
    candidates = list(config["protocol"]["models"])
    values = np.empty((len(candidates), len(metadata), len(SOURCES)))
    lookup = {identity: index for index, identity in enumerate(metadata.sample_id)}
    covered = set()
    for split in splits:
        directory = out / "folds" / split["fold_id"]
        if not checkpoint(directory / "completion.json", out):
            raise ValueError("Missing out-of-fold predictions")
        with np.load(directory / "predictions.npz") as data:
            ids = data["sample_ids"].tolist()
            if data["candidates"].tolist() != candidates or covered & set(ids):
                raise ValueError("CV candidate/sample coverage mismatch")
            values[:, [lookup[identity] for identity in ids]] = data["predictions"]
            covered.update(ids)
    if covered != set(lookup):
        raise ValueError("Incomplete out-of-fold sample coverage")
    reports = {name: summarize(metadata, values[index], config["targets"], config["gates"])
               for index, name in enumerate(candidates)}
    order = sorted(candidates, key=lambda name: ranking_key(name, reports[name]))
    finalists = [name for name in order if name != "refit_control"][:config["protocol"]["n_finalists"]]
    np.savez_compressed(out / "oof_predictions.npz", predictions=values,
                        sample_ids=metadata.sample_id.to_numpy(str), candidates=np.asarray(candidates))
    adapter.write_json(out / "cv_metrics.json", reports)
    adapter.write_json(out / "selection.json", dict(
        status="locked_by_training_only_cv_before_development_evaluation",
        ranking=order, finalists=finalists, development_metrics_used=False,
        criteria=config["protocol"]["ranking"]))
    rows = [dict(candidate=name, cv_rank=order.index(name) + 1,
                 balanced_mean_tv=reports[name]["unified"]["balanced_mean_tv"],
                 pure_mean_tv=reports[name]["pure"]["mean_tv"],
                 pedigree_mean_tv=reports[name]["pedigree"]["mean_tv"],
                 pedigree_within_10pp=reports[name]["pedigree"]["fraction_all_components_within_10pp"],
                 pedigree_within_5pp=reports[name]["pedigree"]["fraction_all_components_within_5pp"])
            for name in order]
    pd.DataFrame(rows).to_csv(out / "cv_comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "cv_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        output_sha256=provenance.fingerprint_files([
            out / "oof_predictions.npz", out / "cv_metrics.json", out / "selection.json",
            out / "cv_comparison.tsv"])))
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


def evaluate(out):
    config = configuration(out, full=True)
    if checkpoint(out / "development_completion.json", out):
        return
    if not checkpoint(out / "cv_completion.json", out):
        raise ValueError("Lock training-only CV before development evaluation")
    selection = json.loads((out / "selection.json").read_text())
    metadata, x = read_table(out / "training.tsv"), np.load(out / "training_features.npy")
    with Path(config["original_model_path"]).open("rb") as stream:
        original = pickle.load(stream)
    all_features = read_table(STACK / "features.tsv")
    evaluation = all_features[all_features.group.isin(["pure", "pedigree", "real"])]
    old_array = all_features[all_features.group.eq("array_imputed")]
    evaluation_x = evaluation[FEATURES].to_numpy(float)
    output_paths, comparisons = [], []
    with threadpool_limits(limits=2):
        models = fit_models(x, metadata[TRUTH].to_numpy(float),
                            metadata.weight_clean.to_numpy(float), config["protocol"])
        for name, model in models.items():
            directory = out / "candidates" / name
            directory.mkdir(parents=True, exist_ok=True)
            bundle = dict(
                kind="phasefree_calibration", sources=SOURCES,
                marker_keys=original["marker_keys"], feature_models=original["feature_models"],
                decoder=model, protocol=config["protocol"],
                experiment_sha256=adapter.digest(out / "experiment.json"),
                input_contract="Complete hard-GT counts in marker_keys order, same frozen base features.")
            path = directory / "model.pkl"
            if path.exists():
                raise ValueError(f"Incomplete full-fit stage must be inspected: {path}")
            with path.open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            output_paths.append(path)
            if name not in ["refit_control", *selection["finalists"]]:
                continue
            predictions = predict(model, evaluation_x)
            if name == "refit_control":
                prior = predict_decoder(original["decoder"], evaluation_x)
                if float(np.max(np.abs(prior - predictions))) > 1e-12:
                    raise ValueError("Architecture-matched full refit does not reproduce the frozen control")
            with path.open("rb") as stream:
                restored = pickle.load(stream)
            if not np.array_equal(predictions, predict(restored["decoder"], evaluation_x)):
                raise ValueError("Full-model serialization changed predictions")
            frame = pd.DataFrame(predictions, columns=[f"pred_{source}" for source in SOURCES])
            frame.insert(0, "sample_id", evaluation.sample_id.to_numpy())
            frame.insert(0, "dataset", evaluation.group.to_numpy())
            frame.to_csv(directory / "predictions.tsv", sep="\t", index=False)
            score(out, name)
            scored = []
            for track in ("pure", "pedigree"):
                part = read_table(directory / f"{track}_predictions_scored.tsv")
                part["track"] = track
                if track == "pure":
                    part["source_group"] = part.source_group.replace({"EAS": "EAS_SEA", "SEA": "EAS_SEA"})
                scored.append(part)
            combined = pd.concat(scored, ignore_index=True)
            report = summarize(combined, combined[[f"pred_{source}" for source in SOURCES]].to_numpy(float),
                               config["targets"], config["gates"])
            array_predictions = predict(model, old_array[FEATURES].to_numpy(float))
            array_summary, _ = accuracy_summary(
                np.tile([.5, .5, 0, 0, 0, 0, 0], (len(old_array), 1)),
                array_predictions, config["targets"])
            report["old_five_eur_asian_array_controls"] = array_summary
            adapter.write_json(directory / "standard_validation.json", report)
            metrics = json.loads((directory / "metrics.json").read_text())
            comparisons.append(dict(
                candidate=name, cv_rank=selection["ranking"].index(name) + 1,
                pure_mean_tv=report["pure"]["mean_tv"],
                pedigree_mean_tv=report["pedigree"]["mean_tv"],
                pedigree_within_10pp=report["pedigree"]["fraction_all_components_within_10pp"],
                pedigree_within_5pp=report["pedigree"]["fraction_all_components_within_5pp"],
                array5_mean_tv=array_summary["mean_tv"], unchanged_gates_pass=metrics["all_gates_pass"],
                array5_diagnostic_pass=array_summary["mean_tv"] <= .10))
            output_paths += list(directory.glob("*scored.tsv")) + [
                directory / "predictions.tsv", directory / "metrics.json",
                directory / "standard_validation.json"]
    pd.DataFrame(comparisons).to_csv(out / "development_comparison.tsv", sep="\t", index=False)
    output_paths.append(out / "development_comparison.tsv")
    adapter.write_json(out / "development_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        selection_sha256=adapter.digest(out / "selection.json"),
        n_full_models=len(models), n_evaluated=3, development_used_to_rank=False,
        expanded_array_status="pending_independent_genotype_verification",
        independent_final_test=False, production_promotion=False,
        output_sha256=provenance.fingerprint_files(output_paths)))
    print(pd.DataFrame(comparisons).to_string(index=False), flush=True)


def array_evaluate(out, arrays):
    config = configuration(out, full=True)
    if checkpoint(out / "array_completion.json", out):
        return
    if not checkpoint(out / "development_completion.json", out):
        raise ValueError("Complete the locked-finalist development stage first")
    verification_path = arrays / "verification.json"
    if not verification_path.exists():
        adapter.write_json(out / "array_dependency.json", dict(
            status="pending", reason="The shared 28-condition genotype stress test is not independently verified.",
            source=str(arrays), independent_final_test=False, production_promotion=False))
        print("Expanded array evaluation pending independent source-genotype verification.", flush=True)
        return
    checked = json.loads(verification_path.read_text())
    if (checked["status"] != "complete_verified_seen_donor_stress_test"
            or checked["model_sha256"] != MODEL_SHA
            or checked["experiment_sha256"] != adapter.digest(arrays / "experiment.json")
            or checked["completion_sha256"] != adapter.digest(arrays / "completion.json")):
        raise ValueError("Expanded array source does not match the verified frozen experiment")
    selection = json.loads((out / "selection.json").read_text())
    array_plan = json.loads((arrays / "experiment.json").read_text())
    with Path(config["original_model_path"]).open("rb") as stream:
        frozen = pickle.load(stream)
    keys, offsets = frozen["marker_keys"], {key: i for i, key in enumerate(frozen["marker_keys"])}
    clean, imputed = np.empty((14, len(keys)), dtype=np.int8), np.empty((28, len(keys)), dtype=np.int8)
    seen = set()
    for chrom in range(1, 23):
        directory = arrays / f"chr{chrom}"
        audit = json.loads((directory / "completion.json").read_text())
        if audit["status"] != "complete" or audit["experiment_sha256"] != adapter.digest(arrays / "experiment.json"):
            raise ValueError("Invalid verified array chromosome checkpoint")
        provenance.verify_hashes(audit["output_sha256"])
        chrom_keys = [tuple(key) for key in json.loads((directory / "marker_keys.json").read_text())]
        if seen & set(chrom_keys):
            raise ValueError("Duplicate shared array marker coverage")
        seen.update(chrom_keys)
        indices = [offsets[key] for key in chrom_keys]
        with np.load(directory / "dosages.npz") as data:
            clean[:, indices], imputed[:, indices] = data["original"], data["array_imputed"]
    if (seen != set(keys) or len(array_plan["samples"]) != 14 or len(array_plan["cases"]) != 28
            or set(array_plan["masks"]) != {"lonni", "frederik"}):
        raise ValueError("Standardized array cohort/marker/mask contract failed")
    training_ids = set(read_table(out / "training.tsv").sample_id)
    if set(array_plan["donor_ids"]) & training_ids:
        raise ValueError("Shared array founders overlap calibration training")
    rows, summaries = [], []
    names = [row["sample_id"] for row in array_plan["samples"]]
    with threadpool_limits(limits=6):
        clean_features = genotype_features(frozen["feature_models"], clean)
        imputed_features = genotype_features(frozen["feature_models"], imputed)
        np.savez_compressed(out / "array_features.npz", clean=clean_features, imputed=imputed_features,
                            sample_ids=np.asarray(names),
                            case_ids=np.asarray([row["case_id"] for row in array_plan["cases"]]))
        for name in ["refit_control", *selection["finalists"]]:
            with (out / f"candidates/{name}/model.pkl").open("rb") as stream:
                model = pickle.load(stream)
            clean_predictions, imputed_predictions = predict(model["decoder"], clean_features), predict(model["decoder"], imputed_features)
            for index, case in enumerate(array_plan["cases"]):
                for stage, values in (
                    ("clean", clean_predictions[names.index(case["sample_id"])]),
                    ("array_imputed", imputed_predictions[index]),
                ):
                    rows.append(dict(
                        candidate=name, sample_id=case["sample_id"], mask=case["mask"],
                        track=case["track"], stage=stage,
                        **dict(zip(TRUTH, case["truth"])), **dict(zip(SOURCES, values))))
    frame = pd.DataFrame(rows)
    for (name, mask, track, stage), group in frame.groupby(["candidate", "mask", "track", "stage"]):
        summary, _ = accuracy_summary(group[TRUTH].to_numpy(float), group[SOURCES].to_numpy(float), config["targets"])
        summaries.append(dict(candidate=name, mask=mask, track=track, stage=stage, **summary))
    frame.to_csv(out / "array_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(out / "array_comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "array_completion.json", dict(
        status="complete", experiment_sha256=adapter.digest(out / "experiment.json"),
        source_verification_sha256=adapter.digest(verification_path),
        source_experiment_sha256=adapter.digest(arrays / "experiment.json"),
        selection_sha256=adapter.digest(out / "selection.json"),
        n_unique_people=14, n_conditions=28, no_new_imputation=True,
        independent_final_test=False, production_promotion=False,
        output_sha256=provenance.fingerprint_files([
            out / "array_predictions.tsv", out / "array_comparison.tsv",
            out / "array_features.npz"])))
    adapter.write_json(out / "array_dependency.json", dict(status="complete", source=str(arrays)))
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def run(out):
    adapter.write_json(out / "sweep_state.json", dict(status="cross_validating"))
    cross_validate(out)
    adapter.write_json(out / "sweep_state.json", dict(status="fitting_locked_finalists"))
    evaluate(out)
    from verify_phasefree_sweep import verify
    verify(out)
    arrays = Path("/validation")
    deadline = time.monotonic() + 48 * 3600
    announced = False
    while not (arrays / "verification.json").exists():
        dependency = arrays / "durable_scheduler.json"
        if dependency.exists():
            state = json.loads(dependency.read_text())
            if state["status"] == "failed":
                raise RuntimeError("Shared array experiment failed; preserve sweep results and recover that dependency")
        if time.monotonic() >= deadline:
            raise TimeoutError("Shared independently verified array results were not available within 48 hours")
        adapter.write_json(out / "sweep_state.json", dict(
            status="waiting_for_shared_array_verification", training_cv_complete=True,
            development_verified=True, production_promotion=False))
        if not announced:
            print("CV/development results verified; waiting for the existing shared array experiment.", flush=True)
            announced = True
        time.sleep(60)
    adapter.write_json(out / "sweep_state.json", dict(status="evaluating_shared_arrays"))
    array_evaluate(out, Path("/validation"))
    from verify_phasefree_sweep import verify_arrays
    verify_arrays(out, arrays)
    adapter.write_json(out / "sweep_state.json", dict(
        status="complete_verified", recommendation_sha256=adapter.digest(out / "recommendation.json"),
        production_promotion=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "smoke", "run", "arrays"), required=True)
    parser.add_argument("--output", type=Path, default=Path("/sweep"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.output)
    elif args.stage == "smoke":
        splits = json.loads((args.output / "folds.json").read_text())
        run_fold(args.output, splits[0]["fold_id"])
        completion = args.output / "folds" / splits[0]["fold_id"] / "completion.json"
        if not checkpoint(completion, args.output):
            raise ValueError("Smoke fold returned without a valid completion checkpoint")
        adapter.write_json(args.output / "smoke_passed.json", dict(
            status="complete", fold=splits[0]["fold_id"],
            experiment_sha256=adapter.digest(args.output / "experiment.json"),
            fold_completion_sha256=adapter.digest(completion)))
    elif args.stage == "arrays":
        array_evaluate(args.output, Path("/validation"))
    else:
        try:
            run(args.output)
        except (ValueError, RuntimeError, OSError) as error:
            adapter.write_json(args.output / "sweep_state.json", dict(
                status="failed", error_type=type(error).__name__, error=str(error),
                production_promotion=False))
            raise
