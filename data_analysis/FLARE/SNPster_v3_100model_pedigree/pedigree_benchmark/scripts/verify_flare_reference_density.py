#!/usr/bin/env python3
"""Independent persisted-table, raw Ped-sim BP and singleton GT verification."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import tempfile

import numpy as np

from realized_ancestry import ancestry_for_batch
from verify_dense_observed import calls, check_case, local_physical


SOURCES = ("EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA")
PROJECT = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def check_hashes(values):
    for path, expected in values.items():
        if digest(path) != expected:
            raise ValueError(f"Changed verification input/output: {path}")


def parent_calls(path, identities, keys):
    values = {identity: bytearray() for identity in identities}
    wanted, found, columns = set(keys), [], None
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("##"):
                continue
            fields = line.rstrip().split("\t")
            if line.startswith("#CHROM"):
                if len(set(fields[9:])) != len(fields[9:]):
                    raise ValueError("Duplicate parent sample header")
                columns = [fields.index(identity, 9) for identity in identities]
                continue
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key not in wanted:
                continue
            if columns is None:
                raise ValueError("Parent VCF lacks sample header")
            found.append(key)
            gt_index = fields[8].split(":").index("GT")
            for identity, column in zip(identities, columns):
                gt = fields[column].split(":")[gt_index]
                if len(gt) != 3 or gt[1] != "|" or gt[0] not in "01" or gt[2] not in "01":
                    raise ValueError("Invalid source haplotype")
                values[identity].append(2 * int(gt[0]) + int(gt[2]))
    if found != keys:
        raise ValueError("Original inherited genotype marker axis differs")
    return values


def verify_truth(selection, protocol):
    root = PROJECT / "results_v3_pedigree_training_screen_01/cohort/pedigree"
    inputs, verified = {}, 0
    with tempfile.TemporaryDirectory() as temporary:
        for chrom in protocol["chromosomes"]:
            for batch in sorted({r["batch"] for r in selection["people"] if r["track"] == "pedigree"}):
                plan = root / f"plans/plans/{batch}/pedigrees.json"
                bp = root / f"chromosomes/pdb_{batch}_chr{chrom}.bp"
                target = Path(temporary) / "truth.tsv"
                ancestry_for_batch(plan, bp, chrom, target)
                actual = {r["sample_id"]: r for r in rows(target)}
                inputs.update({str(path): digest(path) for path in (plan, bp)})
                for person in selection["people"]:
                    if person.get("batch") != batch:
                        continue
                    row = actual[person["sample_id"]]
                    value = [float(row[f"realized_{s}"]) if s != "EAS_SEA"
                             else float(row["realized_EAS"]) + float(row["realized_SEA"]) for s in SOURCES]
                    np.testing.assert_allclose(
                        value, selection["truth"][str(chrom)][person["sample_id"]], rtol=0, atol=1e-12)
                    verified += 1
    return verified, inputs


def verify(out, smoke):
    experiment, selection = read(out / "experiment.json"), read(out / "selection.json")
    check_hashes(experiment["input_sha256"])
    check_hashes(experiment["implementation_sha256"])
    protocol = experiment["protocol"]
    reference_rows = rows("/benchmark/inputs/reference_labels.tsv")
    by_id = {r["sample_id"]: r for r in reference_rows}
    references = selection["references"]
    counts = {name: Counter(by_id[x]["source_group"] for x in identities)
              for name, identities in references.items()}
    if counts["finecap16"] != counts["size_control"]:
        raise ValueError("Reference-size control source counts differ")
    fine = Counter(by_id[x]["fine"] for x in references["finecap16"])
    if max(fine.values()) > protocol["reference_fine_cap"]:
        raise ValueError("Reference fine-population cap exceeded")
    truth_checks, truth_inputs = verify_truth(selection, protocol)
    destination = out / "smoke" if smoke else out
    completion = read(destination / "completion.json")
    check_hashes(completion["output_sha256"])
    people = selection["people"]
    metadata = {r["sample_id"]: r for r in people}
    if len(people) != 316 or len(metadata) != 316:
        raise ValueError("Training people changed")
    seed = protocol["selection_seed"]
    if smoke:
        ranked = lambda identity: hashlib.sha256(f"{seed}:smoke:{identity}".encode()).digest()
        pure = min((r for r in people if r["source_group"] == "SAS"), key=lambda r: ranked(r["sample_id"]))
        child = min((r for r in people if r["track"] == "pedigree"), key=lambda r: (
            -sum(x > 0 for x in selection["truth"]["22"][r["sample_id"]]), ranked(r["sample_id"])))
        identities = {pure["sample_id"], child["sample_id"]}
    else:
        identities = set(metadata)
    candidates = [f"{reference}_{mode}" for reference in protocol["reference_representations"]
                  for mode in ("fixed", "em")]
    groups, total, preserved, posterior_checks = defaultdict(list), 0, 0, 0
    output_hashes = {}
    for chrom, panel in itertools.product([22] if smoke else protocol["chromosomes"], protocol["panels"]):
        directory = destination / f"{panel}_chr{chrom}" if smoke else out / f"tasks/{panel}_chr{chrom}"
        checkpoint = read(directory / "complete.json")
        check_hashes(checkpoint["output_sha256"])
        output_hashes[str(directory / "complete.json")] = digest(directory / "complete.json")
        keys = [(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in rows(directory / "markers.tsv")]
        parent_task = Path("/benchmark/tasks") / f"{'frederik' if panel == 'frederik_thin' else panel}_chr{chrom}"
        cohort = read(parent_task / "cohort.json")
        if panel == "frederik_thin":
            dense_keys = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                          for r in rows(parent_task / "markers.tsv")}
            sparse_keys = rows(Path("/benchmark/tasks") / f"lonni_chr{chrom}/markers.tsv")
            if not set(keys) <= dense_keys or len(keys) != len(sparse_keys):
                raise ValueError("Nested density contrast changed keys/count")
        cases = read(directory / "cases.json")
        if {case["sample_id"] for case in cases} != identities or len(cases) != len(identities):
            raise ValueError("Case identity coverage differs")
        paths = defaultdict(list)
        batches = {r["batch"]: r["vcf"] for r in cohort["pedigree_train"]["batches"]}
        for identity in identities:
            person = metadata[identity]
            path = cohort["pure_train"]["vcf"] if person["track"] == "pure" else batches[person["batch"]]
            paths[path].append(identity)
        source_calls = {}
        for path, ids in paths.items():
            source_calls.update(parent_calls(path, sorted(ids), keys))
        for case in cases:
            n, _, _ = check_case((case, keys))
            preserved += n
            inherited = list(calls(case["inherited_phase"]))[1:]
            expected = source_calls[case["sample_id"]]
            if any(expected[i] != 2 * int(gt[1][1][0]) + int(gt[1][1][2])
                   for i, gt in enumerate(inherited)):
                raise ValueError("Inherited control differs from original person/GT phase")
        predictions = rows(directory / "predictions.tsv")
        expected = set(itertools.product(candidates, protocol["stages"], identities))
        found = {(r["candidate"], r["stage"], r["sample_id"]) for r in predictions}
        if found != expected or len(predictions) != len(expected):
            raise ValueError("Incomplete/duplicate method predictions")
        spot = {min(ids) for ids in (
            [x for x in identities if metadata[x]["track"] == "pure"],
            [x for x in identities if metadata[x]["track"] == "pedigree"])}
        for row in predictions:
            truth = np.array(selection["truth"][str(chrom)][row["sample_id"]])
            p = np.array([float(row[f"pred_{s}"]) for s in SOURCES])
            np.testing.assert_allclose(truth, [float(row[f"truth_{s}"]) for s in SOURCES], rtol=0, atol=1e-12)
            if not np.isfinite(p).all() or (p < 0).any() or abs(p.sum() - 1) > 1e-10:
                raise ValueError("Invalid posterior fractions")
            error = np.abs(p - truth)
            for name, value in (("tv", error.sum() / 2), ("max_component_error", error.max()),
                                ("off_target_mass", p[truth == 0].sum())):
                if abs(value - float(row[name])) > 1e-12:
                    raise ValueError(f"Incorrect saved metric: {name}")
            key = (panel, str(chrom), row["candidate"], row["stage"], row["track"])
            groups[key].append(row)
            if smoke or (row["sample_id"] in spot and row["candidate"] in ("full_fixed", "full_em")):
                prefix = directory / row["candidate"] / f"{row['sample_id']}_{row['stage']}"
                names = [line.split() for line in Path(str(prefix) + ".model").read_text().splitlines()
                         if line.strip() and not line.lstrip().startswith("#")][0]
                value = local_physical(str(prefix) + ".anc.vcf.gz", names,
                                       cohort["chromosome_spans"][str(chrom)], keys)
                np.testing.assert_allclose(p, value, rtol=0, atol=1e-10)
                posterior_checks += 1
        total += len(predictions)
    summaries = rows(destination / "comparison.tsv")
    if len(summaries) != len(groups):
        raise ValueError("Incomplete summary coverage")
    for row in summaries:
        group = groups[tuple(row[k] for k in ("panel", "chromosome", "candidate", "stage", "track"))]
        tv = np.array([float(r["tv"]) for r in group])
        maximum = np.array([float(r["max_component_error"]) for r in group])
        for name, value in (
                ("mean_tv", tv.mean()), ("max_tv", tv.max()), ("p95_tv", np.quantile(tv, .95)),
                ("all_components_within_10pp", (maximum <= .1 + 1e-12).mean()),
                ("all_components_within_5pp", (maximum <= .05 + 1e-12).mean())):
            if abs(float(row[name]) - value) > 1e-12 or int(row["n_people"]) != len(group):
                raise ValueError(f"Saved summary differs: {name}")
    result = dict(status="independently_verified_smoke" if smoke else "independently_verified_training_screen",
                  n_unique_people=len(identities), n_predictions=total, n_summary_groups=len(groups),
                  raw_bp_descendant_truth_vectors=truth_checks, raw_bp_plan_sha256=truth_inputs,
                  unchanged_original_typed_diploid_gts=preserved, independent_local_aggregations=posterior_checks,
                  reference_source_counts={k: dict(v) for k, v in counts.items()}, whole_genome_score=False,
                  verifier_sha256=digest(__file__), block_manifest_sha256=output_hashes)
    (out / ("independent_smoke_verification.json" if smoke else "independent_verification.json")).write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if not k.endswith("sha256")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("/screen"))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    verify(args.outdir, args.smoke)
