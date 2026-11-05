"""Frozen-scope aggregation, calibration-only interpretation, validation-only selection."""
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from core import ancestry_tv_pp, digest, probability, project_simplex, read, write
from models import dump, regression


def aggregate(folders, output):
    output.mkdir()
    metadata = [read(folder / "prediction.json") for folder in folders]
    first = metadata[0]
    expected = first["config"]["chromosomes"]
    if sorted(item["chromosome"] for item in metadata) != sorted(expected):
        raise ValueError("Incomplete or duplicate chromosome predictions")
    predictions, truths, total = [], [], 0
    for folder, item in zip(folders, metadata):
        for key in ("id", "route", "components", "samples", "sources", "config"):
            if item[key] != first[key]:
                raise ValueError(f"Chromosome identity/feature-axis mismatch: {key}")
        if digest(folder / "predictions.npz") != item["sha256"]:
            raise ValueError("Native chromosome prediction artifact changed")
        data = np.load(folder / "predictions.npz")
        size = len(item["components"])
        values = probability(data["values"], size)
        true = probability(data["truth"], len(item["sources"]))
        if values.shape != (len(item["samples"]), size) or true.shape[0] != len(item["samples"]):
            raise ValueError("Prediction person axis differs")
        length = item["span"][1] - item["span"][0]
        if item["span"] != first["config"]["spans"][str(item["chromosome"])]:
            raise ValueError("Chromosome physical denominator changed")
        predictions.append(length * values)
        truths.append(length * true)
        total += length
    np.savez_compressed(output / "aggregate.npz", values=probability(np.sum(predictions, axis=0) / total),
                        truth=probability(np.sum(truths, axis=0) / total))
    write(output / "aggregate.json", {
        **{key: first[key] for key in ("id", "family", "license_track", "components", "sources", "samples", "config", "route")},
        "chromosomes": expected, "physical_span": total,
        "sha256": digest(output / "aggregate.npz"),
    })


def interpret(folder, interpreter_id, output):
    output.mkdir()
    meta = read(folder / "aggregate.json")
    if digest(folder / "aggregate.npz") != meta["sha256"]:
        raise ValueError("Aggregated representation changed")
    config = meta["config"]
    spec = next(item for item in config["interpretations"] if item["id"] == interpreter_id)
    data = np.load(folder / "aggregate.npz")
    x, y = data["values"], data["truth"]
    calibration = [index for index, row in enumerate(meta["samples"]) if row["role"] == "calibration"]
    if len(calibration) < 2:
        raise ValueError("Require at least two disjoint calibration people")
    params, kind = spec["params"], spec["kind"]
    link, reason, trained = "identity", None, None
    if kind == "native":
        model = next(item for item in config["models"] if item["id"] == meta["id"])
        if (set(meta["components"]) == set(meta["sources"]) and len(meta["components"]) == len(meta["sources"])
                and model["params"].get("reference_label") != "fine"):
            values = x[:, [meta["components"].index(source) for source in meta["sources"]]]
        else:
            mapping = model["params"].get("component_map")
            if mapping is None:
                values = x
                reason = "Unlabelled native components retained; no ancestry percentages or native score inferred"
            else:
                mapping = probability(mapping, len(meta["sources"]))
                if mapping.shape != (len(meta["components"]), len(meta["sources"])):
                    raise ValueError("Explicit component map has the wrong named axes")
                values = probability(x @ mapping)
                trained = {"component_map": mapping, "fit_role": (
                    "reference_metadata_only" if model["params"].get("reference_label") == "fine"
                    else "user_supplied_not_evaluation")}
    elif kind == "simplex":
        features, targets = x[calibration], y[calibration]
        dimensions = (x.shape[1], y.shape[1])
        penalty = params.get("l2", .001)
        def objective(flat):
            mapping = flat.reshape(dimensions)
            return float(np.mean((features @ mapping - targets)**2) + penalty * np.mean(mapping**2))
        result = minimize(objective, np.full(dimensions, 1 / y.shape[1]).ravel(),
                          method="SLSQP", bounds=[(0, 1)] * np.prod(dimensions),
                          constraints={"type": "eq", "fun": lambda flat: flat.reshape(dimensions).sum(axis=1) - 1},
                          options={"ftol": 1e-12, "maxiter": 2000})
        if not result.success:
            raise ValueError(f"Training-only simplex mapping failed: {result.message}")
        mapping = probability(result.x.reshape(dimensions), y.shape[1])
        values = probability(x @ mapping)
        trained = {"component_map": mapping, "fit_role": "calibration_only"}
    else:
        trained = regression(kind, params, config["seed"])
        trained.fit(x[calibration], y[calibration])
        values = project_simplex(trained.predict(x))
        link = "Euclidean simplex projection of regression output"
    probability(values)
    dump(output / "interpreter.pkl", trained)
    np.savez_compressed(output / "interpreted.npz", values=values, truth=y)
    write(output / "interpretation.json", {
        **meta, "interpreter": spec, "output_link": link, "unscored_reason": reason,
        "sha256": digest(output / "interpreted.npz"), "interpreter_sha256": digest(output / "interpreter.pkl"),
        "calibration_people": [meta["samples"][index]["sample_id"] for index in calibration] if kind != "native" else [],
    })


def person_metrics(people, sources, seed):
    errors, signed, families, tracks = [], [], [], []
    for identity, items in sorted(people.items()):
        if (len({(item["family_id"], item["track"]) for item in items}) != 1
                or any(item["sources"] != sources for item in items)):
            raise ValueError("Person-level reporting axes differ across masks")
        differences = np.array([item["predictions"] for item in items]) - np.array([item["truth"] for item in items])
        errors.append(float(np.mean([item["tv_percentage_points"] for item in items])))
        signed.append(differences)
        families.append(items[0]["family_id"])
        tracks.append(items[0]["track"])
    if not errors:
        raise ValueError("Cannot report precision without evaluated people")
    values = np.asarray(errors)
    by_family = {family: values[np.array(families) == family] for family in sorted(set(families))}
    sums = np.array([group.sum() for group in by_family.values()])
    sizes = np.array([len(group) for group in by_family.values()])
    intervals = None
    if len(by_family) >= 2:
        weights = np.random.default_rng(seed).multinomial(
            len(by_family), np.full(len(by_family), 1 / len(by_family)), size=2000)
        bootstrap = (weights @ sums) / (weights @ sizes)
        intervals = np.quantile(bootstrap, [.025, .975]).tolist()
    return {
        "people": len(values), "family_components": len(by_family),
        "mean_person_tv_pp": float(values.mean()), "median_person_tv_pp": float(np.median(values)),
        "p95_person_tv_pp": float(np.quantile(values, .95)),
        "family_bootstrap_mean_tv_95pct_ci_pp": intervals,
        "uncertainty_method": "2000 family-component resamples; masks averaged within person; not IID people",
        "per_source": {source: {
            "mae_pp": float(100 * np.mean([np.abs(item[:, index]).mean() for item in signed])),
            "rmse_pp": float(100 * np.sqrt(np.mean([(item[:, index]**2).mean() for item in signed]))),
            "bias_pp": float(100 * np.mean([item[:, index].mean() for item in signed])),
        } for index, source in enumerate(sources)},
        "by_track": {track: {
            "people": tracks.count(track),
            "family_components": len({family for family, subgroup in zip(families, tracks) if subgroup == track}),
            "mean_person_tv_pp": float(values[np.array(tracks) == track].mean()),
            "median_person_tv_pp": float(np.median(values[np.array(tracks) == track])),
            "p95_person_tv_pp": float(np.quantile(values[np.array(tracks) == track], .95)),
        } for track in sorted(set(tracks))},
    }


def select(folders, output):
    output.mkdir()
    records = [read(folder / "interpretation.json") for folder in folders]
    config = records[0]["config"]
    expected = {(route["id"], model["id"], head["id"]) for route in config["routes"]
                for model in config["models"] for head in config["interpretations"]}
    identities = [(row["route"]["id"], row["id"], row["interpreter"]["id"]) for row in records]
    if len(identities) != len(set(identities)) or set(identities) != expected:
        raise ValueError("Incomplete/duplicate route x model x interpretation coverage")
    if config["version"] == 2 and any(
            sample["role"] not in ("calibration", "validation")
            for record in records for sample in record["samples"]):
        raise ValueError("Development interpretation must contain only calibration and validation people")
    predictions, validation, test = [], {}, {}
    for folder, record in zip(folders, records):
        if record["config"] != config or digest(folder / "interpreted.npz") != record["sha256"]:
            raise ValueError("Interpretation provenance mismatch")
        data = np.load(folder / "interpreted.npz")
        for index, sample in enumerate(record["samples"]):
            values = probability(data["values"][index])
            named = not record["unscored_reason"]
            tv = ancestry_tv_pp(values, data["truth"][index], record["sources"],
                combine_eas_sea=config.get("combine_eas_sea_for_metrics", False)) if named else None
            item = {
                "sample_id": sample["sample_id"], "family_id": sample["family_id"], "role": sample["role"],
                "track": sample["track"], "stage": record["route"]["stage"],
                "route": record["route"]["id"], "model": record["id"], "interpreter": record["interpreter"]["id"],
                "license_track": record["license_track"], "tv_percentage_points": tv,
                "eight_source_tv_percentage_points": ancestry_tv_pp(
                    values, data["truth"][index], record["sources"]) if named else None,
                "predictions": values.tolist(), "truth": data["truth"][index].tolist(),
                "sources": record["sources"] if named else record["components"],
            }
            predictions.append(item)
            if not named or sample["role"] not in ("validation", "test"):
                continue
            key = (item["stage"], item["model"], item["interpreter"])
            store = validation if sample["role"] == "validation" else test
            store.setdefault(key, {}).setdefault(sample["sample_id"], []).append(item)
    leaderboard = []
    for key, people in sorted(validation.items()):
        stage, model, interpreter = key
        routes = [route["id"] for route in config["routes"] if route["stage"] == stage]
        values = []
        subgroups = {}
        for sample, items in people.items():
            if sorted(item["route"] for item in items) != sorted(routes):
                raise ValueError(f"Missing/duplicate route prediction for {sample}")
            error = float(np.mean([item["tv_percentage_points"] for item in items]))
            values.append(error)
            subgroups.setdefault(items[0]["track"], []).append(error)
        if not values:
            raise ValueError("No validation people")
        leaderboard.append({
            "stage": stage, "model": model, "interpreter": interpreter, "people": len(values),
            "interpretation_track": "native" if next(item for item in config["interpretations"] if item["id"] == interpreter)["kind"] == "native" else "learned",
            "license_track": next(row["license_track"] for row in records if row["id"] == model),
            "validation_mean_person_tv_pp": float(np.mean(values)),
            "validation_by_track_tv_pp": {key: float(np.mean(value)) for key, value in subgroups.items()},
            "validation_metrics": person_metrics(people, config["sources"], config["seed"]),
        })
    winners = []
    groups = {(row["stage"], row["interpretation_track"], row["license_track"]) for row in leaderboard}
    for group in sorted(groups):
        eligible = [row for row in leaderboard if (row["stage"], row["interpretation_track"], row["license_track"]) == group]
        winner = min(eligible, key=lambda row: (row["validation_mean_person_tv_pp"], row["model"], row["interpreter"]))
        key = (winner["stage"], winner["model"], winner["interpreter"])
        if config["version"] == 2:
            if test:
                raise ValueError("Development selection must not contain test people")
            winners.append({**winner, "test_evaluated": False})
            continue
        test_people = test[key]
        routes = [route["id"] for route in config["routes"] if route["stage"] == winner["stage"]]
        for items in test_people.values():
            if sorted(item["route"] for item in items) != sorted(routes):
                raise ValueError("Incomplete final-test route coverage")
        winners.append({
            **winner, "test_people": len(test_people),
            "selected_test_mean_person_tv_pp": float(np.mean([
                np.mean([item["tv_percentage_points"] for item in items]) for items in test_people.values()])),
            "selected_test_metrics": person_metrics(test_people, config["sources"], config["seed"]),
        })
    # Test labels do not enter ranking or model/decoder fitting.
    write(output / "summary.json", {
        "scope": config["mode"], "scientific_certification": False,
        "reporting_seed": config["seed"],
        "independent_test_declared": config["independent_test"],
        "test_evaluated": config["version"] == 1,
        "metric": "100 * 0.5 * L1; physical chromosomes -> person; masks averaged within person; equal people"
                  + ("; EAS+SEA combined only during scoring" if config.get("combine_eas_sea_for_metrics", False) else ""),
        "combine_eas_sea_for_metrics": config.get("combine_eas_sea_for_metrics", False),
        "validation_leaderboard": leaderboard, "validation_selected_finalists": winners,
        "unscored_native_components": [row["unscored_reason"] for row in records if row["unscored_reason"]],
        "production_promotion": False,
    })
    # Preserve raw evidence, but publish no test leaderboard for unselected candidates.
    write(output / "prediction_evidence.json", predictions)
    with (output / "validation.tsv").open("w", newline="") as handle:
        columns = ["stage", "model", "interpreter", "interpretation_track", "license_track",
                   "people", "validation_mean_person_tv_pp"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(leaderboard)
    files = {path.name: digest(path) for path in sorted(output.iterdir()) if path.is_file()}
    write(output / "complete.json", {"files": files, "complete": True, "candidate_axes": len(expected)})


def verify(output):
    """Recompute published scores, rather than trusting Nextflow process completion."""
    complete = read(output / "complete.json")
    for filename, expected in complete["files"].items():
        if digest(output / filename) != expected:
            raise ValueError(f"Published result changed: {filename}")
    summary = read(output / "summary.json")
    evidence = read(output / "prediction_evidence.json")
    if summary["scope"].startswith("development_") and (
            summary.get("test_evaluated") is not False
            or any(row["role"] == "test" for row in evidence)
            or any("selected_test_metrics" in row for row in summary["validation_selected_finalists"])):
        raise ValueError("Development report contains test evidence or claims")
    for row in summary["validation_leaderboard"] + summary["validation_selected_finalists"]:
        for role, metric in (("validation", "validation_mean_person_tv_pp"),
                             ("test", "selected_test_mean_person_tv_pp")):
            if metric not in row:
                continue
            selected = [item for item in evidence if item["role"] == role
                        and (item["stage"], item["model"], item["interpreter"])
                        == (row["stage"], row["model"], row["interpreter"])]
            people = {}
            reporting_people = {}
            for item in selected:
                value, true = probability(item["predictions"]), probability(item["truth"])
                error = ancestry_tv_pp(value, true, item["sources"],
                    combine_eas_sea=summary.get("combine_eas_sea_for_metrics", False))
                if not np.isclose(error, item["tv_percentage_points"], rtol=0, atol=1e-9):
                    raise ValueError("Published per-person TV differs")
                if "eight_source_tv_percentage_points" in item and not np.isclose(
                        ancestry_tv_pp(value, true, item["sources"]),
                        item["eight_source_tv_percentage_points"], rtol=0, atol=1e-9):
                    raise ValueError("Published separate-source TV differs")
                people.setdefault(item["sample_id"], []).append(error)
                reporting_people.setdefault(item["sample_id"], []).append(item)
            measured = np.mean([np.mean(values) for values in people.values()])
            if not np.isclose(measured, row[metric], rtol=0, atol=1e-9):
                raise ValueError("Independent headline reconstruction failed")
            reporting_key = "validation_metrics" if role == "validation" else "selected_test_metrics"
            if reporting_key in row and row[reporting_key] != person_metrics(
                    reporting_people, selected[0]["sources"], summary["reporting_seed"]):
                raise ValueError("Published generation/source/family uncertainty metrics differ")
    return {"verified": True, "summary_sha256": digest(output / "summary.json"),
            "scope": summary["scope"], "scientific_certification": False}
