"""Three genome-wide native finalists, one locked winner and one test evaluation."""
import argparse
from collections import Counter, defaultdict
import copy
from pathlib import Path
import shutil

import numpy as np

import core
import guided_followup as guided
import projection
from audit import verify_manifest


def source_inputs(spec):
    guided_run, original_run = Path(spec["guided_run"]), Path(spec["original_run"])
    verify_manifest(guided_run / "screen/results")
    original = core.read(original_run / "experiment.json")
    release = Path(original["samples"]).parents[3]
    receipt = core.read(release / "verification/independent.json")
    if not receipt["verified"]:
        raise ValueError("Final comparison requires the independently verified standardized release")
    provenance = core.read(guided_run / "planning/design/provenance.json")
    authority_path = Path(provenance["authority"])
    inventory = core.read(release / "inventory/inventory.json")
    for path in (Path(original["samples"]), Path(original["segments"]), release / "inventory/inventory.json"):
        if core.digest(path) != receipt["files_sha256"][str(path.relative_to(release))]:
            raise ValueError("Original roles, truth or inventory changed")
    if core.digest(authority_path) != inventory["input_sha256"][str(authority_path)]:
        raise ValueError("Original source authority changed")
    authority_rows = [row for row in core.rows(authority_path) if row["role"] == "source"]
    authority = {row["sample"]: row["source_group"] for row in authority_rows}
    if len(authority) != len(authority_rows):
        raise ValueError("Ambiguous source authority")
    return guided_run, original_run, original, release, receipt, authority


def restore_segments(samples, original, release, receipt, authority, config):
    wanted = {row["sample_id"] for row in samples}
    expected = defaultdict(list)
    for row in core.rows(original["segments"]):
        if row["sample_id"] in wanted:
            expected[(row["sample_id"], row["chromosome"], row["haplotype"])].append(row)
    all_segments, hashes = [], {}
    original_samples = core.rows(original["samples"])
    for chrom in config["chromosomes"]:
        path = release / f"panels/observed/genotypes/all/chr{chrom}/master/inheritance.tsv"
        sha = receipt["files_sha256"][str(path.relative_to(release))]
        if core.digest(path) != sha:
            raise ValueError("Original founder-haplotype transmissions changed")
        segments = guided.expand_inheritance(core.rows(path), samples, original_samples, authority,
                                              chrom, config["spans"][str(chrom)])
        count = sum(len(v) for key, v in expected.items() if key[1] == str(chrom))
        if len(segments) != count or any(
                {**row, "source": guided.merged(row["source"])} not in expected[
                    (row["sample_id"], row["chromosome"], row["haplotype"])] for row in segments):
            raise ValueError("Final eight-source truth differs from original founder transmissions")
        all_segments.extend(segments)
        hashes[str(path)] = sha
    return all_segments, hashes


def phase_assets(config, sources):
    assets = {}
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            root = sources[chrom]
            receipt = core.read(root / "independent_verification.json")
            panel = root / f"panels/{route['id']}/chr{chrom}/panel"
            phased = root / f"phased/{route['id']}/chr{chrom}/phased"
            hashes = {}
            for path in (panel / "panel.json", phased / "phase.json", phased / "query.npy"):
                expected = receipt["files"][str(path.relative_to(root))]
                if core.digest(path) != expected:
                    raise ValueError("Reusable phase differs from independent audit")
                hashes[str(path)] = expected
            assets[f"{route['id']}/{chrom}"] = {"panel": str(panel), "phased": str(phased), "hashes": hashes}
    return {"assets": assets, "seed": config["seed"], "reference_filters": config["reference_filters"]}


def write_design(config, samples, segments, output, durable, reuse_sources, native_source):
    output.mkdir()
    config = copy.deepcopy(config)
    config.update(samples=str(durable / "samples.tsv"), segments=str(durable / "segments.tsv"),
                  combine_eas_sea_for_metrics=True)
    core.validate_config(config)
    core.validate_roles(samples, config["sources"], development_only=config["version"] == 2)
    guided.write_rows(output / "samples.tsv", samples, list(samples[0]))
    guided.write_rows(output / "segments.tsv", segments,
                      ["sample_id", "chromosome", "haplotype", "start", "end", "source"])
    core.write(output / "config.json", config)
    core.write(output / "phase-reuse.json", phase_assets(config, reuse_sources))
    projection.build_reuse(config, samples, native_source, output / "native-reuse.json", stopped=False)
    queries = [row for row in samples if row["role"] not in ("reference", "train")]
    jobs = len(config["chromosomes"]) * len(config["routes"]) * len(config["models"])
    reused = sum(len(job["samples"]) for job in core.read(output / "native-reuse.json")["native"].values())
    core.write(output / "scope.json", {
        "roles": dict(Counter(row["role"] for row in samples)), "sources": config["sources"],
        "chromosomes": config["chromosomes"], "models": config["models"],
        "total_singleton_predictions": len(queries) * jobs,
        "verified_native_predictions_reused": reused,
        "new_singleton_calls": len(queries) * jobs - reused,
        "selection_metric": "EAS/SEA merged only during scoring; outputs and truth remain eight-source",
    })
    print(core.read(output / "scope.json"))


def plan(specification, root, output):
    spec = core.read(specification)
    if (spec["maximum_finalists"] != 3 or spec["chromosomes"] != list(range(1, 23))
            or spec["primary_selection_metric"] != "validation_mean_person_EAS_SEA_merged_TV_percentage_points"):
        raise ValueError("Final protocol must lock three genome-wide finalists and merged scoring")
    probe = Path(spec["copying_probe"]) / "results/report"
    result = core.read(probe / "copying-probe.json")
    if (core.digest(probe / "copying-probe.json") != core.read(probe / "complete.json")["report_sha256"]
            or result["accepted_strength"] is not None or result["validation_used"] or result["test_used"]
            or result["native_predictions_independently_reconstructed"] != 258):
        raise ValueError("This final protocol requires the audited rejected copying hypothesis")
    guided_run, original_run, original, release, receipt, authority = source_inputs(spec)
    verify_manifest(original_run / "genome/results")
    config = copy.deepcopy(original)
    previous = core.read(guided_run / "screen/results/frozen/frozen/config.json")
    config.update(version=2, mode="development_full", sources=previous["sources"], independent_test=False,
                  interpretations=previous["interpretations"],
                  models=[next(row for row in previous["models"] if row["id"] == name)
                          for name in ("full_em100", "full_fixed2", "full_em280")])
    samples = core.rows(previous["samples"])
    if Counter(row["role"] for row in samples) != Counter(reference=1908, train=865, calibration=43, validation=127):
        raise ValueError("Final development cohort changed")
    segments, hashes = restore_segments(samples, original, release, receipt, authority, config)
    sources = {chrom: guided_run / "screen/results" if chrom in (20, 22)
               else original_run / "genome/results" for chrom in config["chromosomes"]}
    write_design(config, samples, segments, output, root / "planning/design", sources, guided_run / "screen")
    core.write(output / "truth-provenance.json", {
        "input_sha256": hashes, "original_founder_truth_reconstructed": True,
        "test_evaluated": False, "copying_probe_report_sha256": core.digest(probe / "copying-probe.json"),
    })


def lock_and_test_plan(specification, root, output):
    spec = core.read(specification)
    receipt = verify_manifest(root / "genome/results")
    config = core.read(root / "genome/results/frozen/frozen/config.json")
    summary = core.read(root / "genome/results/summary/summary/summary.json")
    evidence = core.read(root / "genome/results/summary/summary/prediction_evidence.json")
    if (not receipt["verified"] or not summary["combine_eas_sea_for_metrics"] or summary["test_evaluated"]
            or any(row["role"] == "test" for row in evidence) or len(config["models"]) != 3):
        raise ValueError("Cannot lock a winner without the exact independently audited development screen")
    winner = min(summary["validation_leaderboard"],
                 key=lambda row: (row["validation_mean_person_tv_pp"], row["model"]))
    selected = next(row for row in config["models"] if row["id"] == winner["model"])
    lock = {"model": selected, "selection": winner, "development_summary_sha256": receipt["summary_sha256"],
            "selection_metric": "mean person EAS/SEA-merged validation TV", "test_seen_before_lock": False}
    lock_path = root / "testplanning/design/winner-lock.json"
    if lock_path.exists() and core.read(lock_path) != lock:
        raise ValueError("Existing winner lock differs; never change the winner after test planning")
    _, _, original, release, receipt, authority = source_inputs(spec)
    development = core.rows(config["samples"])
    test = [dict(row, source=authority[row["sample_id"]] if row["track"] == "pure" else row["source"])
            for row in core.rows(original["samples"]) if row["role"] == "test"]
    if len(test) != 310 or {row["sample_id"] for row in test} & {row["sample_id"] for row in development}:
        raise ValueError("Original untouched test role changed or overlaps development")
    samples = development + test
    config.update(version=1, mode="full", models=[selected], independent_test=False)
    test_segments, _ = restore_segments(test, original, release, receipt, authority, config)
    segments = core.rows(root / "planning/design/segments.tsv") + test_segments
    sources = {chrom: root / "genome/results" for chrom in config["chromosomes"]}
    write_design(config, samples, segments, output, root / "testplanning/design", sources, root / "genome")
    core.write(output / "winner-lock.json", lock)


def report(root, output):
    development = verify_manifest(root / "genome/results")
    test = verify_manifest(root / "test/results")
    lock = core.read(root / "testplanning/design/winner-lock.json")
    summary = core.read(root / "test/results/summary/summary/summary.json")
    config = core.read(root / "test/results/frozen/frozen/config.json")
    if (config["models"] != [lock["model"]] or development["summary_sha256"] != lock["development_summary_sha256"]
            or len(summary["validation_selected_finalists"]) != 1 or not summary["combine_eas_sea_for_metrics"]):
        raise ValueError("Test result differs from the locked single model or metric")
    winner = summary["validation_selected_finalists"][0]
    if winner["model"] != lock["model"]["id"] or winner["test_people"] != 310:
        raise ValueError("Incomplete or wrong locked test evaluation")
    evidence = core.read(root / "test/results/summary/summary/prediction_evidence.json")
    test_rows = [row for row in evidence if row["role"] == "test"]
    samples = {row["sample_id"]: row for row in core.rows(config["samples"])}
    source_metrics = {}
    for index, source in enumerate(config["sources"]):
        pure = [row for row in test_rows if row["track"] == "pure" and samples[row["sample_id"]]["source"] == source]
        mixed = [row for row in test_rows if row["track"] != "pure"]
        absent = [row for row in test_rows if row["truth"][index] < 1e-10]
        source_metrics[source] = {
            "pure_people": len({row["sample_id"] for row in pure}),
            "pure_mean_predicted_own_source_pp": float(100 * np.mean([row["predictions"][index] for row in pure])) if pure else None,
            "mixed_people": len({row["sample_id"] for row in mixed}),
            "mixed_component_MAE_pp": float(100 * np.mean([
                abs(row["predictions"][index] - row["truth"][index]) for row in mixed])) if mixed else None,
            "absent_source_people": len({row["sample_id"] for row in absent}),
            "absent_source_predicted_at_least_2pct_per_mask_fraction": float(np.mean([
                row["predictions"][index] >= .02 for row in absent])) if absent else None,
        }
    eas, sea = config["sources"].index("EAS"), config["sources"].index("SEA")
    discrimination = {}
    for source, index, other in (("EAS", eas, sea), ("SEA", sea, eas)):
        rows = [row for row in test_rows if row["track"] == "pure" and samples[row["sample_id"]]["source"] == source]
        discrimination[source] = {
            "people": len({row["sample_id"] for row in rows}),
            "correct_per_mask_fraction": float(np.mean([
                row["predictions"][index] > row["predictions"][other] for row in rows])) if rows else None,
        }
    output.mkdir()
    core.write(output / "final.json", {
        "locked_model": lock["model"], "development_selection": lock["selection"],
        "held_out_test": winner, "test_eight_source_mean_TV_pp": float(np.mean([
            row["eight_source_tv_percentage_points"] for row in test_rows])),
        "test_models_evaluated": 1, "test_people": 310,
        "test_used_for_selection_or_tuning": False, "sources": config["sources"],
        "primary_metric_merges_EAS_SEA_only_at_scoring": True,
        "source_precision": source_metrics, "EAS_SEA_discrimination": discrimination,
        "production_promotion": False, "no_further_tuning": True,
        "caveat": "Population-labelled truth; related pedigree components and masks are not IID; not an external clinical validation",
    })
    core.write(output / "complete.json", {"development_summary_sha256": development["summary_sha256"],
        "test_summary_sha256": test["summary_sha256"], "final_report_sha256": core.digest(output / "final.json")})
    print("Final three-model comparison and locked single-winner test independently verified")


def cleanup(root, output):
    if (root / "final-comparison").read_text().strip() != "true" or not (root / "source.sha256").is_file():
        raise ValueError("Cleanup requires the exact dedicated final comparison")
    record = core.read(root / "reporting/report/complete.json")
    for stage, key in (("genome", "development_summary_sha256"), ("test", "test_summary_sha256")):
        if verify_manifest(root / stage / "results")["summary_sha256"] != record[key]:
            raise ValueError("Final audit receipt differs before cleanup")
    final = root / "reporting/report/final.json"
    if core.digest(final) != record["final_report_sha256"]:
        raise ValueError("Final report changed before cleanup")
    paths = [root / stage / "work" for stage in ("planning", "genome", "testplanning", "test", "reporting")]
    if any(path.is_symlink() or path.resolve() != path for path in paths):
        raise ValueError("Refuse redirected final caches")
    removed = []
    for path in paths:
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(str(path))
    verify_manifest(root / "genome/results")
    verify_manifest(root / "test/results")
    if core.digest(final) != record["final_report_sha256"]:
        raise ValueError("Final report changed after cleanup")
    core.write(output, {"removed_exact_successful_caches": removed, "durable_outputs_reverified": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("plan", "test-plan", "report", "cleanup"))
    parser.add_argument("--specification", type=Path)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "plan":
        plan(args.specification, args.root, args.output)
    elif args.stage == "test-plan":
        lock_and_test_plan(args.specification, args.root, args.output)
    elif args.stage == "report":
        report(args.root, args.output)
    else:
        cleanup(args.root, args.output)
