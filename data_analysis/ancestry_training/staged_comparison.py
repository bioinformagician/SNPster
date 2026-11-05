"""Deterministic development scope and validation-only FLARE advancement."""
import argparse
from collections import Counter
import copy
import csv
import hashlib
from pathlib import Path
import shutil

import core


def development_samples(samples, sources, seed, sizes=(32, 64)):
    core.validate_roles(samples, sources)
    chosen = set()
    report = {}
    for role, size in zip(("calibration", "validation"), sizes):
        available = [row for row in samples if row["role"] == role]
        pure_size = size // 2
        if pure_size < len(sources):
            raise ValueError("Development scope cannot represent every pure source")
        quotas = {"pure": pure_size}
        for index, track in enumerate(("mixed_g1", "mixed_g2", "mixed_g3")):
            quotas[track] = (size - pure_size) // 3 + (index < (size - pure_size) % 3)
        families, labels, fine = Counter(), Counter(), Counter()
        selected = []
        for track, count in quotas.items():
            pool = [row for row in available if row["track"] == track]
            if len(pool) < count:
                raise ValueError(f"Insufficient {role}/{track} people for the declared development scope")
            for _ in range(count):
                def priority(row):
                    identity = row["sample_id"]
                    tie = hashlib.sha256(f"{seed}/{role}/{identity}".encode()).hexdigest()
                    return (labels[row["source"]] if track == "pure" else 0,
                            families[row["family_id"]],
                            fine[(row["source"], row["fine"])], tie, identity)
                row = min(pool, key=priority)
                pool.remove(row)
                selected.append(row)
                chosen.add(row["sample_id"])
                families[row["family_id"]] += 1
                labels[row["source"]] += 1
                fine[(row["source"], row["fine"])] += 1
        if ({row["source"] for row in selected if row["track"] == "pure"} != set(sources)
                or set(families) != {row["family_id"] for row in available}):
            raise ValueError(f"Development subset loses a source or biological-family component: {role}")
        report[role] = {
            "people": len(selected), "families": len(families), "by_track": dict(Counter(
                row["track"] for row in selected)),
            "pure_sources": dict(Counter(row["source"] for row in selected if row["track"] == "pure")),
            "family_counts": dict(sorted(families.items())),
        }
    result = [row for row in samples if row["role"] in ("reference", "train") or row["sample_id"] in chosen]
    core.validate_roles(result, sources, development_only=True)
    return result, report


def call_count(config, samples):
    queries = sum(row["role"] in ("calibration", "validation", "test") for row in samples)
    return queries * len(config["routes"]) * len(config["chromosomes"]) * sum(
        model["family"] == "flare" for model in config["models"])


def plan(config_path, source_run, output, durable):
    from projection import build_reuse
    config = core.read(config_path)
    core.validate_config(config)
    if config["mode"] != "full":
        raise ValueError("Staged development requires the original full-autosome comparison config")
    samples = core.rows(config["samples"])
    selected, counts = development_samples(samples, config["sources"], config["seed"])
    output.mkdir()
    with (output / "samples.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)
    full = copy.deepcopy(config)
    full.update(version=2, mode="development_full", independent_test=False,
                samples=str(durable / "samples.tsv"))
    screen = copy.deepcopy(full)
    screen.update(mode="development_screen", chromosomes=[22])
    for key in ("reference", "training", "maps", "spans"):
        screen[key] = {"22": screen[key]["22"]}
    for route in screen["routes"]:
        route["vcfs"] = {"22": route["vcfs"]["22"]}
    for name, value in (("screen.json", screen), ("genome-template.json", full)):
        core.validate_config(value)
        core.write(output / name, value)
    if sum(model["family"] == "flare" for model in full["models"]) < 5:
        raise ValueError("The five-finalist screen requires at least five FLARE candidates")
    core.write(output / "scope.json", {
        "source_config_sha256": core.digest(config_path), "source_samples_sha256": core.digest(config["samples"]),
        "selection": "seeded metadata-only source/fine/family/generation balancing; original role ownership",
        "counts": counts, "sources": config["sources"], "test_people": 0,
        "reference_people": sum(row["role"] == "reference" for row in selected),
        "training_people": sum(row["role"] == "train" for row in selected),
        "screen_flare_calls": call_count(screen, selected),
        "superseded_full_flare_calls": call_count(config, samples),
        "planned_genome_flare_calls": (
            sum(row["role"] in ("calibration", "validation") for row in selected)
            * 5 * len(full["routes"]) * len(full["chromosomes"])),
        "automatic_flare_finalists": 5,
        "test_policy": "no test genotypes or scores in model execution; no automatic test expansion",
        "source_resolution": "frozen source axis preserved; merged truth cannot certify separate EAS/SEA",
    })
    build_reuse(screen, selected, source_run, output / "screen.reuse.json", stopped=True)


def finalists(config, summary, count=5):
    core.validate_config(config)
    if config["version"] != 2 or summary["scope"] != "development_screen" or summary["test_evaluated"]:
        raise ValueError("Only an audited development-screen validation leaderboard may advance models")
    candidates = {model["id"] for model in config["models"] if model["family"] == "flare"}
    heads = {head["id"] for head in config["interpretations"]}
    ranked = []
    for identity in candidates:
        entries = [row for row in summary["validation_leaderboard"] if row["model"] == identity]
        if not entries:
            raise ValueError(f"Missing FLARE validation candidate: {identity}")
        if {row["interpreter"] for row in entries} != heads or len(entries) != len(heads):
            raise ValueError("FLARE advancement requires complete, nonduplicated interpretation coverage")
        if len({row["stage"] for row in entries}) != 1:
            raise ValueError("Do not combine different assay stages to select FLARE finalists")
        ranked.append(min(entries, key=lambda row: (
            row["validation_mean_person_tv_pp"], row["interpreter"])))
    ranked.sort(key=lambda row: (row["validation_mean_person_tv_pp"], row["model"]))
    if count != 5 or len(ranked) < count:
        raise ValueError("Automatic advancement is exactly five FLARE models; expansions require review")
    result = copy.deepcopy(config)
    selected = {row["model"] for row in ranked[:count]}
    result["models"] = [model for model in config["models"]
                        if model["family"] != "flare" or model["id"] in selected]
    core.validate_config(result)
    return result, ranked


def advance(plan_folder, screen_results, output):
    from audit import verify_manifest
    from projection import build_reuse
    verified = verify_manifest(screen_results)
    summary_path = screen_results / "summary/summary/summary.json"
    if not verified["verified"] or verified["summary_sha256"] != core.digest(summary_path):
        raise ValueError("Screen summary differs from independently audited evidence")
    if core.read(screen_results / "frozen/frozen/config.json") != core.read(plan_folder / "screen.json"):
        raise ValueError("Audited screen belongs to a different frozen development plan")
    config, ranked = finalists(core.read(plan_folder / "genome-template.json"), core.read(summary_path))
    samples = core.rows(config["samples"])
    if call_count(config, samples) > 25000:
        raise ValueError("Genome stage exceeds the explicit 25,000-singleton automatic development budget")
    output.mkdir()
    core.write(output / "genome.json", config)
    core.write(output / "finalists.json", {
        "selection_role": "validation_only", "screen_summary_sha256": core.digest(summary_path),
        "selected": ranked[:5], "reserves": ranked[5:],
        "rule": "best validation-only interpretation per distinct FLARE config; then mean-person TV, model ID",
        "flare_singleton_calls": call_count(config, samples), "test_evaluated": False,
        "additional_finalists": "not automatic; require a documented result-based reason and revised cost approval",
        "partial_chromosome_screen": True, "genome_wide_scores_require_all_22_autosomes": True,
    })
    build_reuse(config, samples, screen_results.parent, output / "genome.reuse.json", stopped=False)


def cleanup(root, output):
    from audit import verify_manifest
    root = root.resolve()
    if ((root / "staged-development").read_text().strip() != "true"
            or not (root / "source.sha256").is_file() or not (root / "experiment.json").is_file()):
        raise ValueError("Cleanup requires this exact successful, dedicated staged run")
    verified = {stage: verify_manifest(root / stage / "results") for stage in ("screen", "genome")}
    if any(not row["verified"] for row in verified.values()):
        raise ValueError("Both scientific stages must pass independent verification before cleanup")
    removed = []
    for stage in ("planning", "screen", "advancement", "genome"):
        work = root / stage / "work"
        if work.is_symlink() or work.resolve() != root / stage / "work":
            raise ValueError("Cleanup refuses redirected stage work directories")
        if work.is_dir():
            shutil.rmtree(work)
            removed.append(str(work))
    for stage in ("screen", "genome"):
        verify_manifest(root / stage / "results")
    core.write(output, {"removed_exact_successful_caches": removed, "durable_outputs_reverified": True,
                        "failed_superseded_runs_untouched": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("plan", "advance", "cleanup"))
    parser.add_argument("--config", type=Path)
    parser.add_argument("--source-run", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--screen-results", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--durable", type=Path)
    args = parser.parse_args()
    if args.stage == "plan":
        plan(args.config, args.source_run, args.output, args.durable)
    elif args.stage == "advance":
        advance(args.plan, args.screen_results, args.output)
    else:
        cleanup(args.source_run, args.output)
