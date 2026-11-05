#!/usr/bin/env python3
"""Independent CSV/VCF checks of the completed observed-only ancestry pilot."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


SOURCES = ("EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA")


def read(path):
    return json.loads(Path(path).read_text())


def records(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_hashes(values):
    for path, expected in values.items():
        if digest(path) != expected:
            raise ValueError(f"Saved scientific artifact changed: {path}")


def calls(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                if len(fields) != 10:
                    raise ValueError("Expected singleton observed-only input")
                yield ("sample", fields[9]), None
                continue
            gt = fields[9].split(":")[fields[8].split(":").index("GT")]
            if len(gt) != 3 or gt[0] not in "01" or gt[2] not in "01" or gt[1] not in "|/":
                raise ValueError("Incomplete or invalid observed diploid GT")
            yield (fields[0], int(fields[1]), fields[3], fields[4]), (int(gt[0]) + int(gt[2]), gt)


def check_case(arguments):
    case, markers = arguments
    known = list(calls(case["inherited_phase"]))
    typed = list(calls(Path(case["inherited_phase"]).with_name("typed.unphased.vcf.gz")))
    phased = list(calls(case["array_phase_only"]))
    if not len(known) == len(typed) == len(phased) == len(markers) + 1:
        raise ValueError("Phase-only marker count changed")
    for a, b, c, expected in zip(known[1:], typed[1:], phased[1:], markers):
        if (a[0] != b[0] or a[0] != c[0] or a[0] != expected
                or a[1][0] != b[1][0] or a[1][0] != c[1][0] or c[1][1][1] != "|"):
            raise ValueError("Supplied observed GT/key changed or missing phase")
    if not known[0][0] == typed[0][0] == phased[0][0] == ("sample", case["sample_id"]):
        raise ValueError("Observed sample axis changed")
    heterozygotes = [(a[1][1][0], b[1][1][0]) for a, b in zip(known[1:], phased[1:]) if a[1][0] == 1]
    orientation = [a != b for a, b in heterozygotes]
    switches = sum(a != b for a, b in zip(orientation, orientation[1:]))
    if switches != case["phase_diagnostic"]["switches"]:
        raise ValueError("Saved phase diagnostic differs")
    return len(markers), switches, max(0, len(heterozygotes) - 1)


def local_physical(path, ancestry_names, span, markers):
    positions, probabilities = [], []
    order = [SOURCES.index(name) for name in ancestry_names]
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key != markers[len(positions)]:
                raise ValueError("Local posterior exact-allele axis mismatch")
            names, values = fields[8].split(":"), fields[9].split(":")
            p = np.array([[float(v) for v in values[names.index(name)].split(",")]
                          for name in ("ANP1", "ANP2")])
            p /= p.sum(axis=1)[:, None]
            vector = np.zeros(7)
            vector[order] = p.mean(axis=0)
            probabilities.append(vector)
            positions.append(key[1])
    if len(positions) != len(markers) or np.any(np.diff(positions) <= 0):
        raise ValueError("Incomplete local posterior coverage")
    boundaries = np.r_[span[0], (np.array(positions[:-1]) + positions[1:]) / 2, span[1]]
    return np.average(probabilities, axis=0, weights=np.diff(boundaries))


def verify(out):
    experiment, completion = read(out / "experiment.json"), read(out / "completion.json")
    verify_hashes(experiment["implementation_sha256"])
    verify_hashes(experiment["input_sha256"])
    verify_hashes(read(out / "recovery.json")["implementation_sha256"])
    verify_hashes(completion["output_sha256"])
    protocol = experiment["protocol"]
    legacy = read("/legacy/experiment.json")
    spans = read("/results/cohort.json")["chromosome_spans"]
    truth = {sample["sample_id"]: np.asarray(sample["truth"], float) for sample in legacy["samples"]}
    samples = {sample["sample_id"]: sample for sample in legacy["samples"]}
    selection = read(out / "selection.json")
    if (set(selection["reference_ids"]) & set(selection["pure_ids"])
            or set(legacy["donor_ids"]) & set(selection["reference_ids"] + selection["pure_ids"])):
        raise ValueError("Evaluation founder leakage")
    donors = read(Path(__file__).resolve().parents[2] / "results_v3_phasefree_rich_01/donors.json")
    for batch in selection["pedigree_train"]:
        for row in records(batch["samples_tsv"]):
            if not set(donors[row["family_id"]]) <= set(selection["pure_ids"]):
                raise ValueError("Unavailable training founder retained")
    totals, weight, phase_jobs, phase_tracks = {}, {}, [], []
    spot_checks, calls_reused = 0, 0
    coverage = {}
    for mask in protocol["masks"]:
        observed = {(c.removeprefix("chr"), int(p), r, a)
                    for line in Path(legacy["masks"][mask]).read_text().splitlines()
                    for c, p, r, a in [line.split("\t")]}
        coverage[mask] = dict(typed=len(observed), retained=0)
        for chrom in protocol["chromosomes"]:
            directory = out / "tasks" / f"{mask}_chr{chrom}"
            for name in ("complete.json", "inputs_complete.json"):
                verify_hashes(read(directory / name)["output_sha256"])
            markers = [(r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                       for r in records(directory / "markers.tsv")]
            if not set(markers) <= observed or len(set(markers)) != len(markers):
                raise ValueError("Added/untyped or duplicate input marker")
            coverage[mask]["retained"] += len(markers)
            cases = read(directory / "cases.json")
            if {c["sample_id"] for c in cases} != set(truth) or len(cases) != 14:
                raise ValueError("Case cohort changed")
            for case in cases:
                phase_jobs.append((case, markers))
                phase_tracks.append((mask, case["track"]))
            task_rows = records(directory / "predictions.tsv")
            expected = {(candidate, stage, person) for candidate in protocol["candidates"]
                        for stage in ("inherited_phase", "array_phase_only") for person in truth}
            found = {(r["candidate"], r["stage"], r["sample_id"]) for r in task_rows}
            if found != expected or len(task_rows) != len(expected):
                raise ValueError("Incomplete or duplicate native-method task coverage")
            for row in task_rows:
                p = np.array([float(row[f"pred_{s}"]) for s in SOURCES])
                if (not np.isfinite(p).all() or np.any(p < 0)
                        or not np.isclose(p.sum(), 1, rtol=0, atol=1e-10)):
                    raise ValueError("Invalid native posterior fraction")
                bp = spans[str(chrom)][1] - spans[str(chrom)][0]
                if int(row["span_bp"]) != bp:
                    raise ValueError("Incorrect physical-span weight")
                key = (row["candidate"], mask, row["stage"], row["sample_id"])
                totals[key] = totals.get(key, np.zeros(7)) + p * bp
                weight[key] = weight.get(key, 0) + bp
            if chrom == 22:
                for case in [next(c for c in cases if c["track"] == track) for track in ("pure", "pedigree")]:
                    for stage in ("inherited_phase", "array_phase_only"):
                        prefix = directory / "flare_array_em" / f"{case['sample_id']}_{stage}"
                        names = [line.split() for line in Path(str(prefix) + ".model").read_text().splitlines()
                                 if line.strip() and not line.startswith("#")][0]
                        p = local_physical(str(prefix) + ".anc.vcf.gz", names, spans["22"], markers)
                        saved = next(r for r in task_rows if r["candidate"] == "flare_array_em"
                                     and r["stage"] == stage and r["sample_id"] == case["sample_id"])
                        np.testing.assert_allclose(p, [float(saved[f"pred_{s}"]) for s in SOURCES], rtol=0, atol=1e-10)
                        spot_checks += 1
    with ProcessPoolExecutor(max_workers=4) as pool:
        phase_values = list(pool.map(check_case, phase_jobs, chunksize=4))
    n_gt = sum(value[0] for value in phase_values)
    phase_summary = []
    for mask in protocol["masks"]:
        for track in ("pure", "pedigree"):
            subset = [v for key, v in zip(phase_tracks, phase_values) if key == (mask, track)]
            pairs = sum(v[2] for v in subset)
            phase_summary.append(dict(mask=mask, track=track, switches=sum(v[1] for v in subset),
                                      pairs=pairs, pooled_switch_error=sum(v[1] for v in subset) / pairs))
    global_rows = records(out / "predictions.tsv")
    if len(global_rows) != len(totals):
        raise ValueError("Incomplete final global rows")
    by_group, details = {}, []
    for row in global_rows:
        key = (row["candidate"], row["mask"], row["stage"], row["sample_id"])
        p = totals[key] / weight[key]
        np.testing.assert_allclose(p, [float(row[f"pred_{s}"]) for s in SOURCES], rtol=0, atol=1e-12)
        absolute = np.abs(p - truth[row["sample_id"]])
        tv, worst = absolute.sum() / 2, absolute.max()
        if abs(tv - float(row["tv"])) > 1e-12 or abs(worst - float(row["worst_component_error"])) > 1e-12:
            raise ValueError("Per-person error differs")
        track = "pure" if samples[row["sample_id"]]["track"] == "pure" else "pedigree"
        if row["track"] != track:
            raise ValueError("Changed pure/pedigree track")
        group = key[:3] + (track,)
        by_group.setdefault(group, []).append((tv, worst, absolute, float((p * (truth[row["sample_id"]] < .01)).sum())))
        dominant = sorted(zip(SOURCES, p), key=lambda item: -item[1])[:3]
        details.append(dict(candidate=key[0], mask=key[1], stage=key[2], sample_id=key[3],
                            track=track, tv=tv, worst_component_error=worst,
                            true_sources=";".join(f"{s}:{v:.3f}" for s, v in zip(SOURCES, truth[key[3]]) if v > 0),
                            leading_predictions=";".join(f"{s}:{v:.3f}" for s, v in dominant)))
    comparisons = records(out / "comparison.tsv")
    if len(comparisons) != len(by_group):
        raise ValueError("Missing group summaries")
    for row in comparisons:
        group = (row["candidate"], row["mask"], row["stage"], row["track"])
        values = by_group[group]
        tv, worst = np.array([v[0] for v in values]), np.array([v[1] for v in values])
        expected = dict(mean_tv=tv.mean(), p95_tv=np.quantile(tv, .95),
                        p95_worst_component_error=np.quantile(worst, .95),
                        n_all_components_within_10pp=int((worst <= .10 + 1e-12).sum()),
                        n_all_components_within_5pp=int((worst <= .05 + 1e-12).sum()),
                        off_target_mass=np.mean([v[3] for v in values]))
        expected.update({f"mae_{s}": np.mean([v[2][i] for v in values]) for i, s in enumerate(SOURCES)})
        for name, value in expected.items():
            if abs(float(row[name]) - value) > 1e-12:
                raise ValueError(f"Summary differs: {group}/{name}")
    with (out / "per_person_diagnostic.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(details[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(details)
    result = dict(status="complete_independently_verified_observed_only_development_pilot",
                  n_tasks=44, n_autosomes=22, n_unique_people=14, n_masks=2,
                  n_chromosome_method_predictions=44 * 5 * 2 * 14,
                  n_global_method_predictions=len(global_rows),
                  n_supplied_typed_diploid_gt_independently_verified=n_gt,
                  n_local_flare_physical_aggregation_spot_checks=spot_checks,
                  no_added_untyped_snps=True, no_changed_diploid_gt=True,
                  original_and_augmented_metrics_independently_recomputed=True,
                  coverage=coverage, phase_summary=phase_summary,
                  independent_final_test=False, imputation_performed=False, production_promotion=False,
                  completion_sha256=digest(out / "completion.json"),
                  recovery_sha256=digest(out / "recovery.json"),
                  verifier_sha256=digest(Path(__file__)))
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    print("Worst primary FLARE EM pure controls:", flush=True)
    for row in sorted([r for r in details if r["candidate"] == "flare_array_em"
                       and r["stage"] == "array_phase_only" and r["track"] == "pure"],
                      key=lambda r: -r["tv"])[:8]:
        print(row, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("/benchmark"))
    verify(parser.parse_args().outdir)
