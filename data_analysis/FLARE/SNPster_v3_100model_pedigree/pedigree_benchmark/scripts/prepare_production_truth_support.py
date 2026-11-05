#!/usr/bin/env python3
"""Recover observed assay support without the old four-dataset marker intersection."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from audit_production_ancestry import digest, write_json
from production_matched_preprocessing import read, reference, positions, read_call_matrix, write_vcf
from sgdp_genetic_qc import query_keys, sample_ids
from verify_sgdp_full_manifest import require_hashes

BASE = Path("/srv/dependencies/ancestry")
RAW = BASE / "GSE77508_ProcessedMatrix.txt.gz"


def assay(output):
    output.mkdir()
    counts = Counter()
    with gzip.open(RAW, "rt") as source, gzip.open(output / "assay.tsv.gz", "wt") as destination:
        fields = source.readline().rstrip("\n").split("\t")
        names = [field.removesuffix(".GType") for field in fields if field.endswith(".GType")]
        if len(names) != 196 or len(set(names)) != len(names) or "Bajo25KN" not in names:
            raise ValueError("Unexpected original GSE77508 assay sample axis")
        indices = [(fields.index(name + ".GType"), fields.index(name + ".Score")) for name in names]
        destination.write("rsid\t" + "\t".join(names) + "\n")
        for line in source:
            row = line.rstrip("\n").split("\t")
            if len(row) != len(fields):
                raise ValueError("Malformed original assay row")
            if not re.fullmatch(r"rs[0-9]+", row[0]):
                counts["non_rsid_rows"] += 1
                continue
            calls = []
            for genotype, score in indices:
                call = row[genotype]
                if call in ("AA", "AB", "BA", "BB"):
                    quality = float(row[score])
                    if not np.isfinite(quality) or not 0 <= quality <= 1:
                        raise ValueError("Invalid original assay genotype score")
                    if quality < .8:
                        call = "--"
                        counts["low_score_genotypes"] += 1
                elif call not in ("--", "NC", "NoCall", "NOCALL", ""):
                    raise ValueError(f"Unknown raw assay genotype: {call}")
                else:
                    call = "--"
                    counts["missing_genotypes"] += 1
                calls.append(call)
            destination.write(row[0] + "\t" + "\t".join(calls) + "\n")
            counts["rsid_rows"] += 1
    write_json(output / "assay.json", dict(
        samples=names, minimum_score=.8, counts=dict(counts),
        input_sha256={str(RAW): digest(RAW)}))


def raw_dosage(path, identities):
    result = {}
    with gzip.open(path, "rt") as handle:
        samples = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            identity, _, rest = line.rstrip("\n").partition("\t")
            if identity in identities:
                if identity in result:
                    raise ValueError("Duplicate raw assay rsID")
                result[identity] = np.asarray([
                    {"AA": 0, "AB": 1, "BA": 1, "BB": 2, "--": -1}[value]
                    for value in rest.split("\t")], dtype=np.int8)
    return samples, result


def reconstruct_assay(raw, phased, candidate):
    """Learn A/B orientation only from other people, never the evaluation donor."""
    valid = (raw >= 0) & (phased >= 0)
    valid[candidate] = False
    if np.count_nonzero(valid & (raw != 1)) < 4:
        if raw[candidate] != 1:
            return None
        observed = 1
    else:
        forward = np.array_equal(raw[valid], phased[valid])
        reverse = np.array_equal(2 - raw[valid], phased[valid])
        if forward == reverse:
            raise ValueError("Original assay and published genotypes disagree on A/B allele encoding")
        observed = int(raw[candidate] if forward else 2 - raw[candidate])
    if phased[candidate] >= 0 and observed != phased[candidate]:
        raise ValueError("Published Bajo genotype differs from independently reconstructed original assay")
    return observed


def observed_completeness(calls, assay_eligible):
    eligible = np.asarray(assay_eligible)
    if (calls.ndim != 3 or calls.shape[1] != 2 or eligible.shape != (calls.shape[2],)
            or eligible.dtype != np.dtype(bool) or not np.isin(calls, [-1, 0, 1]).all()):
        raise ValueError("Invalid observed genotypes or assay eligibility axis")
    denominator = int(eligible.sum())
    if denominator == 0:
        raise ValueError("No loci pass original assay quality/orientation checks")
    called = (calls >= 0).all(axis=1)
    retained = called.all(axis=0)
    if (retained & ~eligible).any():
        raise ValueError("Complete genotype outside original assay eligibility")
    return retained, dict(
        source_compatible_loci=calls.shape[2], assay_eligible_loci=denominator,
        assay_excluded_loci=int((~eligible).sum()), complete_loci=int(retained.sum()),
        complete_fraction_of_eligible=float(retained.sum() / denominator),
        complete_fraction_of_source_compatible=float(retained.mean()),
        per_person_called_loci=called[:, eligible].sum(axis=1).tolist())


def chromosome(cohort_path, chromosome, assay_folder, run_dir, output):
    output.mkdir()
    cohort = read(cohort_path)
    people = cohort["evaluation_founders"] + cohort["training_founders"]
    gse = BASE / f"harmonized/gse73996/vcf/gse73996.chr{chromosome}.grch38.vcf.gz"
    bajo = BASE / f"gse77508/GSE77508.chr{chromosome}.phased.vcf.gz"
    public = Path(cohort["normalized_vcfs"][str(chromosome)])
    shared = set(query_keys(gse, chromosome)) & set(query_keys(bajo, chromosome))
    region = positions(output / "positions.tsv", sorted(shared))
    shared &= set(query_keys(public, chromosome, region))
    shared &= set(query_keys(reference(chromosome), chromosome, region))
    keys = sorted(shared, key=lambda key: key[1])
    if len(keys) < 100:
        raise ValueError("Insufficient expanded exact observed truth support")
    ids = sample_ids(bajo)
    if ids != read(assay_folder / "assay.json")["samples"]:
        if set(ids) != set(read(assay_folder / "assay.json")["samples"]):
            raise ValueError("Published and original assay sample identities differ")
    stream = subprocess.check_output(
        ["bcftools", "query", "-R", str(region), "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n", str(bajo)],
        text=True)
    rsids = {}
    for line in stream.splitlines():
        c, p, r, a, identity = line.split()
        key = (c, int(p), r, a)
        if key in shared:
            if key in rsids:
                raise ValueError("Duplicated published assay allele key")
            rsids[key] = identity
    raw_names, raw = raw_dosage(assay_folder / "assay.tsv.gz", set(rsids.values()))
    reorder = [raw_names.index(identity) for identity in ids]
    bajo_calls = read_call_matrix(bajo, keys, ids, output)
    candidate = ids.index("Bajo25KN")
    dosage = np.full(len(keys), -1, dtype=np.int8)
    excluded = Counter()
    for index, key in enumerate(keys):
        original = raw.get(rsids[key])
        if original is None or original[reorder[candidate]] < 0:
            excluded["missing_or_low_quality_original_assay"] += 1
            continue
        values = original[reorder]
        published = bajo_calls[:, :, index]
        target = published.sum(axis=1)
        target[(published < 0).any(axis=1)] = -1
        reconstructed = reconstruct_assay(values, target, candidate)
        if reconstructed is None:
            excluded["ambiguous_raw_assay_allele_orientation"] += 1
        else:
            dosage[index] = reconstructed
    calls = np.full((len(people), 2, len(keys)), -1, dtype=np.int8)
    gse_ids = [i for i, row in enumerate(people[:8])
               if not row["public"] and row["sample_id"] != "Bajo25KN"]
    public_ids = [i for i, row in enumerate(people[:8]) if row["public"]]
    for path, indices in ((gse, gse_ids), (public, public_ids),
                          (reference(chromosome), list(range(8, len(people))))):
        calls[indices] = read_call_matrix(path, keys, [people[i]["sample_id"] for i in indices], output)
    bajo_index = next(i for i, row in enumerate(people) if row["sample_id"] == "Bajo25KN")
    calls[bajo_index, 0] = np.where(dosage < 0, -1, dosage == 2)
    calls[bajo_index, 1] = np.where(dosage < 0, -1, dosage > 0)
    retained, completeness = observed_completeness(calls, dosage >= 0)
    selected = [key for key, keep in zip(keys, retained) if keep]
    counts = {}
    for mask, path in cohort["masks"].items():
        typed = {(c, int(p), r, a) for c, p, r, a in (line.split() for line in Path(path).read_text().splitlines())}
        counts[mask] = len(set(selected) & typed)
    write_json(output / "eligibility.json", dict(
        chromosome=chromosome, completeness=completeness, typed_counts=counts,
        people=[row["sample_id"] for row in people], assay_exclusions=dict(excluded),
        minimum_assay_score=.8, minimum_post_assay_completeness=.9,
        minimum_typed_markers_per_mask=100, cohort_sha256=digest(cohort_path)))
    if completeness["complete_fraction_of_eligible"] < .9:
        raise ValueError(
            f"Only {len(selected)}/{completeness['assay_eligible_loci']} eligible loci have complete "
            f"observed genotypes after assay QC on chr{chromosome}; require90%")
    for mask in cohort["masks"]:
        if counts[mask] < 100:
            raise ValueError(f"Expanded observed truth still has fewer than100 {mask}/chr{chromosome} markers")
    write_vcf(output / f"observed_chr{chromosome}.vcf", selected,
              [row["sample_id"] for row in people], calls[:, :, retained])
    sources = {str(path): digest(path) for path in (gse, bajo, public, reference(chromosome))}
    sources.update({str(run_dir / "assay/assay" / name): digest(assay_folder / name)
                    for name in ("assay.json", "assay.tsv.gz")})
    write_json(output / "support.json", dict(
        chromosome=chromosome, original_common_loci=len(keys), retained_observed_loci=len(selected),
        typed_counts=counts, completeness=completeness, excluded_assay_loci=dict(excluded),
        source_sha256=sources,
        original_assay_orientation="at_least_four_other_informative_people; never_evaluation_donor",
        training_truth="fixed_dense_training_role_genotypes_not_evaluation_truth",
        cohort_sha256=digest(cohort_path)))


def collect(cohort_path, run_dir, folders, output):
    output.mkdir()
    cohort = read(cohort_path)
    require_hashes(cohort["input_sha256"])
    support, sources = {}, {}
    if len(folders) != 22:
        raise ValueError("All22 chromosome support outputs required")
    for folder in folders:
        record = read(folder / "support.json")
        chromosome = record["chromosome"]
        if chromosome not in range(1, 23) or str(chromosome) in support or record["cohort_sha256"] != digest(cohort_path):
            raise ValueError("Incorrect expanded support chromosome/cohort binding")
        require_hashes(record["source_sha256"])
        vcf = folder / f"observed_chr{chromosome}.vcf.gz"
        destination = run_dir / f"chromosomes/chr{chromosome}/support"
        for path in (vcf, Path(str(vcf) + ".csi"), folder / "support.json", folder / "eligibility.json"):
            sources[str(destination / path.name)] = digest(path)
        sources.update(record["source_sha256"])
        support[str(chromosome)] = str(destination / vcf.name)
    cohort["observed_truth_vcfs"] = support
    cohort["input_sha256"].update(sources)
    cohort["input_sha256"][str(cohort_path)] = digest(cohort_path)
    cohort["truth_support_recovery"] = dict(
        original_cohort=str(cohort_path), reason="old_four_dataset_backbone_has_only75_Lonni_chr22_markers",
        original_protocol_and_donor_selection_unchanged=True, no_imputed_evaluation_truth=True,
        observed_completeness_denominator="source_compatible_loci_after_original_assay_quality_and_orientation_qc")
    write_json(output / "cohort.json", cohort)
    write_json(output / "complete.json", dict(
        status="all22_observed_truth_support_verified", published_sha256=sources,
        cohort_sha256=digest(output / "cohort.json")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("assay", "chromosome", "collect"))
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--assay", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--folders", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.stage == "assay":
        assay(args.output)
    elif args.stage == "chromosome":
        chromosome(args.cohort, args.chromosome, args.assay, args.run_dir, args.output)
    else:
        collect(args.cohort, args.run_dir, args.folders, args.output)
