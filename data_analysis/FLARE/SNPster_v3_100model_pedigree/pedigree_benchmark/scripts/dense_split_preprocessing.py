#!/usr/bin/env python3
"""Synthetic template genomes and array uploads for the approved dense split."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np

from audit_production_ancestry import SOURCES, digest, write_json
import model_family_gnomix as native
import production_matched_preprocessing as base
from sgdp_genetic_qc import query_keys
from verify_sgdp_full_manifest import require_hashes


def mask_keys(cohort, mask, chromosome):
    span = cohort["chromosome_spans"][str(chromosome)]
    points = native.map_points(base.ROOT / "plink.GRCh38.map" /
                               f"plink.chr{chromosome}.GRCh38.map", chromosome)
    return {(fields[0], int(fields[1]), fields[2], fields[3])
            for line in Path(cohort["masks"][mask]).read_text().splitlines()
            if (fields := line.split()) and fields[0] == str(chromosome)
            and span[0] <= int(fields[1]) < span[1]
            and points[0, 0] <= int(fields[1]) <= points[-1, 0]}


def rsids(cohort_path, chromosomes, output):
    cohort = base.read(cohort_path)
    require_hashes(cohort["input_sha256"])
    canonical = base.PROJECT / "results_v3_model_family_screen_05/inputs/pure_pure.canonical.vcf.gz"
    keys = set()
    for chromosome in chromosomes:
        span = cohort["chromosome_spans"][str(chromosome)]
        keys.update(key for key in query_keys(canonical, chromosome)
                    if span[0] <= key[1] < span[1])
        for mask in cohort["protocol"]["masks"]:
            keys.update(mask_keys(cohort, mask, chromosome))
    base.bind_rsids(cohort_path, chromosomes, sorted(keys, key=lambda key: (int(key[0]), *key[1:])), output)


def synthesise(cohort, chromosome, keys, identities, founder_calls, people):
    if founder_calls.shape != (len(identities), 2, len(keys)) or not np.isin(founder_calls, [0, 1]).all():
        raise ValueError("Synthetic genomes require complete, explicitly paired template homologs")
    offsets = {identity: index for index, identity in enumerate(identities)}
    positions = np.asarray([key[1] for key in keys])
    calls = np.full((len(people), 2, len(keys)), -1, dtype=np.int8)
    labels = np.full(calls.shape, -1, dtype=np.int8)
    for person_index, person in enumerate(people):
        covered = np.zeros((2, len(keys)), dtype=np.int8)
        for segment in cohort["local_ancestry"][str(chromosome)][person["sample_id"]]:
            first, end = np.searchsorted(positions, [segment["start"], segment["end"]])
            hap = segment["haplotype"]
            donor_hap = segment["donor_haplotype"]
            if hap not in (0, 1) or donor_hap not in (0, 1) or segment["donor"] not in offsets:
                raise ValueError("Synthetic inheritance references an unknown donor or homolog")
            calls[person_index, hap, first:end] = founder_calls[offsets[segment["donor"]], donor_hap, first:end]
            labels[person_index, hap, first:end] = SOURCES.index(segment["source"])
            covered[hap, first:end] += 1
        if not np.all(covered == 1):
            raise ValueError("Synthetic genotype/truth inheritance does not cover every key exactly once")
    return calls, labels


def founders(cohort_path, chromosome, markers, output):
    output.mkdir()
    cohort, rows = base.read(cohort_path), base.marker_rows(markers)
    source = base.reference(chromosome)
    require_hashes({str(source): cohort["template_sha256"][str(source)]})
    region = base.positions(output / "template_positions.tsv", [row[:4] for row in rows])
    available = query_keys(source, chromosome, region)
    counts = Counter(key[1] for key in available)
    available = {key for key in available if counts[key[1]] == 1}
    rows = [row for row in rows if row[:4] in available]
    keys = [row[:4] for row in rows]
    if len(keys) < 100:
        raise ValueError("Insufficient exact synthetic template SNPs")
    people = cohort["training_founders"] + cohort["evaluation_founders"]
    identities = [row["sample_id"] for row in people]
    calls = base.read_call_matrix(source, keys, identities, output, phased=True)
    np.savez_compressed(output / "template.npz", calls=calls)
    write_json(output / "founders.json", dict(
        chromosome=chromosome, people=people, founder_ids=identities, keys=keys,
        rsids=[row[4] for row in rows], template_sha256=digest(source),
        truth_design=cohort["truth_design"], cohort_sha256=digest(cohort_path)))


def make_cases(cohort_path, chromosome, folder, masks, output):
    output.mkdir()
    cohort, data = base.read(cohort_path), base.read(folder / "founders.json")
    if data["cohort_sha256"] != digest(cohort_path):
        raise ValueError("Synthetic template extraction belongs to a different frozen cohort")
    keys = [tuple(key) for key in data["keys"]]
    calls, labels = synthesise(cohort, chromosome, keys, data["founder_ids"],
                              np.load(folder / "template.npz")["calls"], cohort["people"])
    base.write_vcf(output / f"truth_chr{chromosome}.vcf", keys,
                   [row["sample_id"] for row in cohort["people"]], calls, phased=True)
    np.savez_compressed(output / "local_truth.npz", labels=labels)
    (output / "pieces").mkdir()
    jobs, coverage = [], {}
    for mask in masks:
        typed = mask_keys(cohort, mask, chromosome)
        use = [i for i, key in enumerate(keys) if key in typed]
        fraction = len(use) / len(typed) if typed else 0
        if len(use) < 100 or fraction < .9:
            raise ValueError(f"{mask}/chr{chromosome}: genuine mask/rsID/template support "
                             f"{len(use)}/{len(typed)} ({fraction:.4%}) is below the unchanged90% gate")
        coverage[mask] = dict(expected=len(typed), retained=len(use), fraction=fraction)
        for person_index, person in enumerate(cohort["people"]):
            identity = person["production_ids"][mask]
            with (output / "pieces" / f"{identity}.txt").open("w") as handle:
                for i in use:
                    chrom, position, ref, alt = keys[i]
                    genotype = "".join((ref, alt)[int(value)] for value in calls[person_index, :, i])
                    handle.write(f"{data['rsids'][i]}\t{chrom}\t{position}\t{genotype}\n")
            jobs.append(dict(mask=mask, production_id=identity,
                             sample_id=person["sample_id"], markers=len(use)))
    write_json(output / "cases.json", dict(chromosome=chromosome, jobs=jobs, keys=keys,
                                           mask_coverage=coverage, truth_design=cohort["truth_design"]))


def aligned_training_labels(truth, inferred, labels, positions):
    if truth.shape != inferred.shape or labels.shape != truth.shape or truth.ndim != 3 or truth.shape[1] != 2:
        raise ValueError("Training phase-alignment axes disagree")
    result = labels.copy()
    evidence = []
    positions = np.asarray(positions)
    for index, (expected, actual) in enumerate(zip(truth, inferred)):
        informative = np.flatnonzero((expected.sum(axis=0) == actual.sum(axis=0)) &
                                     (expected[0] != expected[1]))
        if not len(informative):
            if not np.array_equal(labels[index, 0], labels[index, 1]):
                raise ValueError("Mixed training homologs cannot be aligned without informative genotypes")
            evidence.append(0)
            continue
        insertion = np.searchsorted(positions[informative], positions)
        left = informative[np.maximum(insertion - 1, 0)]
        right = informative[np.minimum(insertion, len(informative) - 1)]
        nearest = np.where(positions - positions[left] <= positions[right] - positions, left, right)
        swapped = expected[0, nearest] != actual[0, nearest]
        result[index, 0, swapped] = labels[index, 1, swapped]
        result[index, 1, swapped] = labels[index, 0, swapped]
        evidence.append(int(len(informative)))
    return result, evidence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("rsids", "founders", "cases", "assemble", "verify"))
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--chromosomes", type=int, nargs="+")
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--markers", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--folders", type=Path, nargs="+")
    parser.add_argument("--masks", nargs="+")
    parser.add_argument("--mask")
    parser.add_argument("--merged", type=Path)
    parser.add_argument("--imputed", type=Path)
    parser.add_argument("--qc", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "rsids":
        rsids(args.cohort, args.chromosomes, args.output)
    elif args.stage == "founders":
        founders(args.cohort, args.chromosome, args.markers, args.output)
    elif args.stage == "cases":
        make_cases(args.cohort, args.chromosome, args.folder, args.masks, args.output)
    elif args.stage == "assemble":
        base.assemble(args.cohort, args.folders, args.masks, args.output)
    else:
        base.verify(args.cohort, args.chromosome, args.mask, args.folder, args.merged,
                    args.imputed, args.qc, args.output)
