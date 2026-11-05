"""Prepared development-only fine-population comparison; never score country labels."""
import argparse
from collections import Counter, defaultdict
import copy
from pathlib import Path
import shutil

import numpy as np

import core
import guided_followup as guided
import models
import projection
from audit import verify_manifest
from evaluation import person_metrics


def population_design(samples, sources):
    reference = [row for row in samples if row["role"] == "reference"]
    groups = {}
    for row in reference:
        fine = row["fine"]
        if fine in groups and groups[fine] != row["source"]:
            raise ValueError("A fine reference population maps to multiple ancestry sources")
        groups[fine] = row["source"]
    components = sorted(groups)
    mapping = [[float(groups[name] == source) for source in sources] for name in components]
    params = {"reference_label": "fine", "components": components, "component_map": mapping}
    models.flare_reference_labels({"reference": reference, "sources": sources}, params, list(range(len(reference))))
    counts = Counter(row["fine"] for row in reference)
    return params, [{"population": name, "source": groups[name], "reference_people": counts[name],
                     "small_reference_warning": counts[name] < 10} for name in components]


def validate_specification(spec):
    if (spec["sources"] != ["EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"]
            or spec["expected_roles"] != {"reference": 1908, "train": 865, "calibration": 43, "validation": 127}
            or spec["candidate_id"] != "fine_population"):
        raise ValueError("Population comparison must retain the original eight-source development roles")
    if spec["version"] == 1:
        if spec["chromosomes"] != [20, 22] or spec["maximum_new_singleton_calls"] != 680:
            raise ValueError("Post-winner population follow-up requires chr20/22 and exactly 680 new calls")
    elif spec["version"] == 2:
        if (spec["chromosomes"] != [22] or spec["maximum_new_singleton_calls"] != 340
                or spec["baseline_id"] != "full_em100" or spec["baseline_policy"] != "predeclared_guided_EM100"):
            raise ValueError("Parallel population pilot requires only chr22, fixed guided EM100 and 340 new calls")
    else:
        raise ValueError("Unsupported population comparison protocol version")


def prepare(specification, output):
    spec = core.read(specification)
    validate_specification(spec)
    root = Path(spec["guided_run"])
    receipt = core.read(root / "screen/results/independent_verification.json")
    config_path = root / "screen/results/frozen/frozen/config.json"
    if not receipt["verified"] or core.digest(config_path) != receipt["files"]["frozen/frozen/config.json"]:
        raise ValueError("Prepared reference labels require an independently audited eight-source config")
    config = core.read(config_path)
    sample_path = Path(config["samples"])
    frozen = core.read(root / "screen/results/frozen/frozen/manifest.json")
    if core.digest(sample_path) != frozen["input_sha256"][str(sample_path)]:
        raise ValueError("Audited reference/sample labels changed")
    samples = core.rows(sample_path)
    if dict(Counter(row["role"] for row in samples)) != spec["expected_roles"]:
        raise ValueError("Population follow-up changes the original development cohort")
    core.validate_roles(samples, spec["sources"], development_only=True)
    authority_path = Path(core.read(root / "planning/design/provenance.json")["authority"])
    authority_sha = core.read(root / "planning/design/provenance.json")["input_sha256"][str(authority_path)]
    if core.digest(authority_path) != authority_sha:
        raise ValueError("Original fine-label authority changed")
    source_rows = [row for row in core.rows(authority_path) if row["role"] == "source"]
    authority = {row["sample"]: row for row in source_rows}
    if len(authority) != len(source_rows):
        raise ValueError("Ambiguous original population authority")
    for row in samples:
        if row["role"] == "reference":
            label = authority[row["sample_id"]]
            if row["fine"] != label["fine"] or row["source"] != label["source_group"]:
                raise ValueError("Reference fine/source labels disagree with the original authority")
    params, populations = population_design(samples, spec["sources"])
    output.mkdir()
    guided.write_rows(output / "reference-populations.tsv", populations, list(populations[0]))
    core.write(output / "population-parameters.json", params)
    core.write(output / "prepared.json", {
        "specification_sha256": core.digest(specification), "sources": spec["sources"],
        "reference_populations": len(populations), "reference_people": sum(row["reference_people"] for row in populations),
        "populations_by_source": dict(Counter(row["source"] for row in populations)),
        "small_populations": [row for row in populations if row["small_reference_warning"]],
        "reference_sample_manifest_sha256": core.digest(sample_path), "authority_sha256": authority_sha,
        "expected_new_singleton_calls": spec["maximum_new_singleton_calls"],
        "reference_counts_are_not_independent_accuracy_evidence": True,
        "test_accessed": False, "scientific_followup_started": False,
    })
    print(f"Prepared {len(populations)} original reference populations; fixed eight-source sums; no scientific run or test access")
    return spec, config, samples, params


def completed_predecessor(root):
    completion = root / "reporting/report/complete.json"
    cleanup = root / "reports/scoped-cache-cleanup.json"
    lock = root / "testplanning/design/winner-lock.json"
    if not completion.is_file() or not cleanup.is_file() or not lock.is_file():
        raise ValueError("Wait for the current final sweep, locked test and successful cleanup to finish; do not start a parallel follow-up")
    if not core.read(cleanup)["durable_outputs_reverified"]:
        raise ValueError("Predecessor cleanup has not reverified durable outputs")
    if (root / "java-profile").read_text().strip() != "c1":
        raise ValueError("Population follow-up requires the same pinned C1 runtime")
    results = root / "genome/results"
    receipt = verify_manifest(results)
    if not receipt["verified"]:
        raise ValueError("Predecessor development predictions lack independent verification")
    winner = core.read(lock)
    summary = core.read(results / "summary/summary/summary.json")
    if (winner["development_summary_sha256"] != receipt["summary_sha256"] or winner["test_seen_before_lock"]
            or summary["test_evaluated"] or not summary["combine_eas_sea_for_metrics"]):
        raise ValueError("Predecessor winner was not locked on development-only merged metrics")
    selected = min(summary["validation_leaderboard"], key=lambda row: (
        row["validation_mean_person_tv_pp"], row["model"]))
    config = core.read(results / "frozen/frozen/config.json")
    model = next(row for row in config["models"] if row["id"] == selected["model"])
    if model != winner["model"] or model["family"] != "flare" or set(model["params"]) - {"em", "gen", "states", "update-p"}:
        raise ValueError("Predecessor lock differs from the native full-reference validation winner")
    return config, model, receipt


def baseline_source(spec):
    validate_specification(spec)
    if spec["version"] == 1:
        predecessor = Path(spec["predecessor"])
        config, model, receipt = completed_predecessor(predecessor)
        return predecessor, predecessor / "genome/results", config, model, receipt
    source = Path(spec["guided_run"])
    if (source / "java-profile").read_text().strip() != "c1":
        raise ValueError("Parallel pilot must use the same pinned C1 runtime as its completed baseline")
    results = source / "screen/results"
    receipt = verify_manifest(results)
    config = core.read(results / "frozen/frozen/config.json")
    summary = core.read(results / "summary/summary/summary.json")
    if not receipt["verified"] or summary["test_evaluated"] or config["version"] != 2:
        raise ValueError("Parallel pilot requires an audited development-only guided baseline")
    model = next(row for row in config["models"] if row["id"] == spec["baseline_id"])
    if (model["family"] != "flare" or model["params"] != {
            "em": True, "gen": 10, "states": 100, "update-p": False}):
        raise ValueError("Parallel pilot must keep the predeclared native full-reference EM100 settings")
    return source, results, config, model, receipt


def plan(specification, root, output):
    spec = core.read(specification)
    predecessor, results, original_config, baseline, receipt = baseline_source(spec)
    spec, _, samples, params = prepare(specification, output)
    if (core.rows(original_config["samples"]) != samples or original_config["sources"] != spec["sources"]
            or len(original_config["routes"]) != 2 or any(route["stage"] != "harmonized" for route in original_config["routes"])):
        raise ValueError("Population comparison must preserve the matched observed development cohort")
    config = copy.deepcopy(original_config)
    candidate = {"id": spec["candidate_id"], "family": "flare", "params": {**baseline["params"], **params}}
    config.update(version=2, mode="development_screen", chromosomes=spec["chromosomes"],
                  models=[baseline, candidate], interpretations=[{"id": "native", "kind": "native", "params": {}}],
                  independent_test=False, combine_eas_sea_for_metrics=True,
                  samples=str(root / "planning/design/samples.tsv"), segments=str(root / "planning/design/segments.tsv"))
    for key in ("reference", "training", "maps", "spans"):
        config[key] = {str(chrom): config[key][str(chrom)] for chrom in spec["chromosomes"]}
    for route in config["routes"]:
        route["vcfs"] = {str(chrom): route["vcfs"][str(chrom)] for chrom in spec["chromosomes"]}
    core.validate_config(config)
    guided.write_rows(output / "samples.tsv", samples, list(samples[0]))
    segments_path = Path(original_config["segments"])
    manifest = core.read(results / "frozen/frozen/manifest.json")
    if core.digest(segments_path) != manifest["input_sha256"][str(segments_path)]:
        raise ValueError("Predecessor founder-grounded eight-source truth changed")
    segments = [row for row in core.rows(segments_path) if int(row["chromosome"]) in spec["chromosomes"]]
    guided.write_rows(output / "segments.tsv", segments, list(segments[0]))
    core.write(output / "config.json", config)
    phases, native = {}, {}
    queries = [row for row in samples if row["role"] in ("calibration", "validation")]
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            key = f"{route['id']}/{chrom}"
            panel = results / f"panels/{route['id']}/chr{chrom}/panel"
            phased = results / f"phased/{route['id']}/chr{chrom}/phased"
            prediction = results / f"predictions/{route['id']}/chr{chrom}/{baseline['id']}/prediction"
            meta = core.read(panel / "panel.json")
            phase_paths = (panel / "panel.json", phased / "phase.json", phased / "query.npy")
            phases[key] = {"panel": str(panel), "phased": str(phased),
                           "hashes": {str(path): receipt["files"][str(path.relative_to(results))] for path in phase_paths}}
            captured = {}
            for person in queries:
                record = projection.capture(prediction, meta, baseline["params"], person["sample_id"])
                if record is None:
                    raise ValueError("Complete predecessor baseline lacks native singleton evidence")
                for name, sha in record["files"].items():
                    if sha != receipt["files"][str((prediction / name).relative_to(results))]:
                        raise ValueError("Reusable baseline evidence differs from its independent receipt")
                captured[person["sample_id"]] = record
            native[f"{key}/{baseline['id']}"] = {"signature": projection.signature(panel, baseline["params"]),
                                                "samples": captured}
            native[f"{key}/{candidate['id']}"] = {"signature": None, "samples": {}}
    calls = len(queries) * len(config["routes"]) * len(config["chromosomes"])
    if calls != spec["maximum_new_singleton_calls"]:
        raise ValueError("Fine-population comparison differs from its frozen new-call scope")
    core.write(output / "phase-reuse.json", {"assets": phases, "seed": config["seed"],
                                            "reference_filters": config["reference_filters"]})
    core.write(output / "native-reuse.json", {
        "version": 2, "config": config, "java_profile": "c1", "source_run": str(predecessor),
        "projections": {}, "native": native,
        "reuse_policy": "Only audited baseline singletons; new native fine-panel predictions; exact phased GT/key/reference checks",
    })
    core.write(output / "scope.json", {
        "new_singleton_calls": calls, "baseline_singletons_reused": calls, "development_queries": len(queries),
        "validation_people": sum(row["role"] == "validation" for row in queries),
        "baseline": baseline, "candidate": candidate["id"], "chromosomes": config["chromosomes"],
        "baseline_policy": spec.get("baseline_policy", "final_validation_lock"),
        "populations": len(params["components"]), "test_evaluated": False,
        "predecessor_test_results_read": False, "automatic_expansion": False,
        "predecessor_development_summary_sha256": receipt["summary_sha256"],
    })
    print(f"Population comparison ready: {calls} new candidate calls, {calls} reused baseline calls; no test or automatic expansion")


def report(root, output):
    results = root / "screen/results"
    receipt = verify_manifest(results)
    if not receipt["verified"]:
        raise ValueError("Population report requires an independent native prediction audit")
    config = core.read(results / "frozen/frozen/config.json")
    evidence = core.read(results / "summary/summary/prediction_evidence.json")
    if config["version"] != 2 or any(row["role"] == "test" for row in evidence):
        raise ValueError("Population comparison must never contain test evidence")
    design = root / "planning/design"
    scope = core.read(design / "scope.json")
    samples = core.rows(config["samples"])
    metadata = {row["sample_id"]: row for row in samples}
    params, populations = population_design(samples, config["sources"])
    candidate = next(row for row in config["models"] if row["id"] == scope["candidate"])
    for name in params:
        if candidate["params"][name] != params[name]:
            raise ValueError("Candidate mapping differs from reference-only population metadata")
    reference = [row for row in samples if row["role"] == "reference"]
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            base_id = next(row["base_id"] for row in models.plan(config)["models"] if row["id"] == scope["candidate"])
            base = results / f"bases/{route['id']}/chr{chrom}/{base_id}/base"
            labels = [line.split() for line in (base / "reference.panel").read_text().splitlines()]
            if labels != [[row["sample_id"], row["fine"]] for row in reference]:
                raise ValueError("Native candidate did not fit exactly the declared fine-panel reference")
    models_report = []
    people_errors = defaultdict(dict)
    for model in config["models"]:
        selected = [row for row in evidence if row["model"] == model["id"] and row["role"] == "validation"]
        people = defaultdict(list)
        for row in selected:
            people[row["sample_id"]].append(row)
        if (len(people) != scope["validation_people"] or any(
                sorted(row["route"] for row in rows) != sorted(route["id"] for route in config["routes"])
                for rows in people.values())):
            raise ValueError("Population report lacks complete validation person/mask coverage")
        metrics = {}
        for source_index, source in enumerate(config["sources"]):
            pure = [row for row in selected if row["track"] == "pure" and metadata[row["sample_id"]]["source"] == source]
            mixed = [row for row in selected if row["track"] != "pure"]
            absent = [row for row in selected if row["truth"][source_index] < 1e-10]
            metrics[source] = {
                "pure_people": len({row["sample_id"] for row in pure}),
                "pure_mean_own_source_pp": float(100 * np.mean([row["predictions"][source_index] for row in pure])) if pure else None,
                "mixed_component_MAE_pp": float(100 * np.mean([
                    abs(row["predictions"][source_index] - row["truth"][source_index]) for row in mixed])) if mixed else None,
                "absent_source_people": len({row["sample_id"] for row in absent}),
                "absent_source_calls_at_least_2pct_per_mask_fraction": float(np.mean([
                    row["predictions"][source_index] >= .02 for row in absent])) if absent else None,
            }
        eas, sea = config["sources"].index("EAS"), config["sources"].index("SEA")
        discrimination = {}
        for source, own, other in (("EAS", eas, sea), ("SEA", sea, eas)):
            rows = [row for row in selected if row["track"] == "pure" and metadata[row["sample_id"]]["source"] == source]
            discrimination[source] = float(np.mean([
                row["predictions"][own] > row["predictions"][other] for row in rows])) if rows else None
        route_metrics = {route["id"]: float(np.mean([row["tv_percentage_points"] for row in selected
                                                   if row["route"] == route["id"]])) for route in config["routes"]}
        for identity, rows in people.items():
            people_errors[identity][model["id"]] = float(np.mean([row["tv_percentage_points"] for row in rows]))
        models_report.append({
            "model": model["id"], "merged_metrics": person_metrics(people, config["sources"], config["seed"]),
            "separate_source_mean_TV_pp": float(np.mean([row["eight_source_tv_percentage_points"] for row in selected])),
            "source_precision": metrics, "EAS_SEA_per_mask_discrimination": discrimination,
            "merged_TV_pp_by_route": route_metrics,
        })
    families = defaultdict(list)
    for identity, errors in people_errors.items():
        families[metadata[identity]["family_id"]].append(errors[scope["baseline"]["id"]] - errors[scope["candidate"]])
    groups = list(families.values())
    sums, sizes = np.asarray([sum(g) for g in groups]), np.asarray([len(g) for g in groups])
    draws = np.random.default_rng(config["seed"]).integers(0, len(groups), size=(20000, len(groups)))
    bootstrap = sums[draws].sum(axis=1) / sizes[draws].sum(axis=1)
    output.mkdir()
    core.write(output / "population-comparison.json", {
        "models": models_report, "reference_populations": populations,
        "paired_baseline_minus_candidate_merged_TV_pp": float(sums.sum() / sizes.sum()),
        "family_bootstrap_95pct_CI_pp": np.quantile(bootstrap, [.025, .975]).tolist(),
        "validation_people": scope["validation_people"], "family_components": len(groups), "development_only": True,
        "chromosomes": config["chromosomes"],
        "baseline_policy": scope.get("baseline_policy", "final_validation_lock"),
        "test_evaluated": False, "predecessor_test_results_read": False, "production_promotion": False,
        "fine_population_or_country_accuracy_claimed": False, "automatic_expansion": False,
        "caveat": "Repeated development cohort; partial-chromosome development scores are not genome-wide accuracy. "
                  "The final sweep's locked test belongs to its own model only, not a fresh test for this comparison.",
    })
    core.write(output / "complete.json", {"summary_sha256": receipt["summary_sha256"],
               "report_sha256": core.digest(output / "population-comparison.json")})
    print("Independent fine-native integration and eight-source summation verified; development-only comparison complete")


def cleanup(root, output):
    if (root / "population-followup").read_text().strip() != "true" or not (root / "source.sha256").is_file():
        raise ValueError("Cleanup requires the exact dedicated population follow-up")
    root = root.resolve()
    receipt = core.read(root / "reporting/report/complete.json")
    report_path = root / "reporting/report/population-comparison.json"
    if (verify_manifest(root / "screen/results")["summary_sha256"] != receipt["summary_sha256"]
            or core.digest(report_path) != receipt["report_sha256"]):
        raise ValueError("Population report/audit changed before cache cleanup")
    paths = [root / stage / "work" for stage in ("planning", "screen", "reporting")]
    if any(path.is_symlink() or path.resolve() != path for path in paths):
        raise ValueError("Refuse redirected population work caches")
    removed = []
    for path in paths:
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(str(path))
    verify_manifest(root / "screen/results")
    if core.digest(report_path) != receipt["report_sha256"]:
        raise ValueError("Population report changed after cleanup")
    core.write(output, {"removed_exact_successful_caches": removed, "durable_outputs_reverified": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "check-ready", "plan", "report", "cleanup"))
    parser.add_argument("--specification", type=Path)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.specification, args.output)
    elif args.stage == "check-ready":
        spec = core.read(args.specification)
        _, _, _, model, _ = baseline_source(spec)
        print(f"Audited baseline {model['id']} ready; launch is allowed, but not started by this check")
    elif args.stage == "plan":
        plan(args.specification, args.root, args.output)
    elif args.stage == "report":
        report(args.root, args.output)
    else:
        cleanup(args.root, args.output)
