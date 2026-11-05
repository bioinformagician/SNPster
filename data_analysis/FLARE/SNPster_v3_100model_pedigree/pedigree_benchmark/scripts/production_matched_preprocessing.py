#!/usr/bin/env python3
"""Task implementations for the unchanged production preprocessing benchmark."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from audit_production_ancestry import digest, write_json
from sgdp_genetic_qc import genotype_rows, query_keys, sample_ids, run
from verify_sgdp_full_manifest import require_hashes

PROJECT = Path(__file__).resolve().parents[2]
PVAR = Path("/srv/dependencies/imputation_runner/harmonizer/all_hg38.pvar.zst")
ROOT = Path("/srv/dependencies/imputation_runner/imputer")


def read(path):
    return json.loads(Path(path).read_text())


def reference(chromosome):
    return ROOT / "beagle_references" / (
        f"hgdp1kgp_chr{chromosome}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz")


def compress(path):
    subprocess.run(["bgzip", str(path)], check=True)
    result = Path(str(path) + ".gz")
    subprocess.run(["bcftools", "index", str(result)], check=True)
    return result


def positions(path, keys):
    path.write_text("".join(f"{key[0]}\t{key[1]}\t{key[1]}\n" for key in keys))
    return path


def write_vcf(path, keys, identities, calls, phased=False):
    if calls.shape != (len(identities), 2, len(keys)) or not np.isin(calls, [0, 1]).all():
        raise ValueError("VCF writing requires complete paired biallelic genotypes")
    separator = "|" if phased else "/"
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chromosome in sorted({key[0] for key in keys}, key=int):
            handle.write(f"##contig=<ID={chromosome}>\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                     "\t".join(identities) + "\n")
        for index, (chromosome, position, ref, alt) in enumerate(keys):
            genotypes = [f"{pair[0]}{separator}{pair[1]}" for pair in calls[:, :, index]]
            handle.write(f"{chromosome}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" +
                         "\t".join(genotypes) + "\n")
    return compress(path)


def read_call_matrix(path, keys, identities, directory, phased=False):
    wanted = {key: index for index, key in enumerate(keys)}
    region = positions(directory / "query_positions.tsv", keys)
    result = np.full((len(identities), 2, len(keys)), -1, dtype=np.int8)
    seen = set()
    for key, calls in genotype_rows(path, region, identities):
        if key not in wanted:
            continue
        if key in seen:
            raise ValueError("Duplicate exact genotype key")
        seen.add(key)
        for index, call in enumerate(calls):
            if call in ("./.", ".|."):
                if phased:
                    raise ValueError("Missing production haplotype; filling is forbidden")
                continue
            if not re.fullmatch(r"[01][|/][01]", call) or (phased and "|" not in call):
                raise ValueError(f"Invalid or unphased required genotype: {call}")
            result[index, :, wanted[key]] = [int(call[0]), int(call[2])]
    if seen != set(keys):
        raise ValueError(f"Missing exact genotype loci: {len(set(keys) - seen)}")
    return result


def rsids(cohort_path, chromosomes, output):
    cohort = read(cohort_path)
    require_hashes(cohort["input_sha256"])
    canonical = PROJECT / "results_v3_model_family_screen_05/inputs/pure_pure.canonical.vcf.gz"
    keys = [key for chromosome in chromosomes for key in query_keys(
        Path(cohort["observed_truth_vcfs"][str(chromosome)]) if "observed_truth_vcfs" in cohort else canonical,
        chromosome)]
    bind_rsids(cohort_path, chromosomes, keys, output)


def bind_rsids(cohort_path, chromosomes, keys, output):
    output.mkdir()
    regions = output / "loci.range"
    regions.write_text("".join(f"{c}\t{p}\t{p}\tlocus{i}\n" for i, (c, p, _, _) in enumerate(keys)))
    run(["plink2", "--pvar", PVAR, "--extract", "range", regions, "--make-just-pvar",
         "--threads", "1", "--memory", "8000", "--out", output / "rsids"], output / "extract_rsids")
    wanted, mapping = set(keys), {}
    with (output / "rsids.pvar").open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chromosome, position, identity, ref, alts = line.split()[:5]
            if not re.fullmatch(r"rs[0-9]+", identity):
                continue
            for alt in alts.split(","):
                key = (chromosome, int(position), ref, alt)
                if key in wanted and (key not in mapping or int(identity[2:]) < int(mapping[key][2:])):
                    mapping[key] = identity
    id_counts = Counter(mapping.values())
    mapping = {key: identity for key, identity in mapping.items() if id_counts[identity] == 1}
    for chromosome in chromosomes:
        selected = [(key, mapping[key]) for key in keys if key[0] == str(chromosome) and key in mapping]
        if len(selected) < 100:
            raise ValueError(f"Insufficient real rsID-supported truth SNPs on chr{chromosome}")
        (output / f"chr{chromosome}.tsv").write_text(
            "".join("\t".join(map(str, (*key, identity))) + "\n" for key, identity in selected))
    write_json(output / "rsid_binding.json", dict(
        original_loci=len(keys), real_rsid_loci=len(mapping), pvar_sha256=digest(PVAR),
        cohort_sha256=digest(cohort_path), chromosomes=chromosomes))


def marker_rows(path):
    with Path(path).open() as handle:
        return [(fields[0], int(fields[1]), fields[2], fields[3], fields[4])
                for line in handle if (fields := line.split())]


def founders(cohort_path, chromosome, markers, output):
    output.mkdir()
    cohort, rows = read(cohort_path), marker_rows(markers)
    keys = [row[:4] for row in rows]
    people = cohort["evaluation_founders"] + cohort["training_founders"]
    inputs = PROJECT / "results_v3_model_family_screen_05/inputs"
    specifications = [
        (inputs / "pure_pure.canonical.vcf.gz", [i for i, row in enumerate(people[:8]) if not row["public"]]),
        (inputs / "pure_train.canonical.vcf.gz", list(range(8, len(people)))),
        (Path(cohort["normalized_vcfs"][str(chromosome)]), [i for i, row in enumerate(people[:8]) if row["public"]]),
    ]
    if "observed_truth_vcfs" in cohort:
        specifications = [(Path(cohort["observed_truth_vcfs"][str(chromosome)]), list(range(len(people))))]
    region = positions(output / "candidate_positions.tsv", keys)
    shared = set(keys) & set(query_keys(reference(chromosome), chromosome, region))
    unavailable = {}
    for path, indices in specifications:
        if indices:
            available = set(query_keys(path, chromosome, region))
            unavailable[str(path)] = len(set(keys) - available)
            shared &= available
    keys = [key for key in keys if key in shared]
    rows = [row for row in rows if row[:4] in shared]
    values = np.full((len(people), 2, len(keys)), -1, dtype=np.int8)
    for path, indices in specifications:
        if indices:
            values[indices] = read_call_matrix(path, keys, [people[i]["sample_id"] for i in indices], output)
    retained = (values >= 0).all(axis=(0, 1))
    original_count = len(marker_rows(markers))
    if retained.sum() / original_count < .9:
        raise ValueError("Observed truth/reference intersection retains fewer than90% of known loci")
    selected_keys = [key for key, keep in zip(keys, retained) if keep]
    values = values[:, :, retained]
    ids = [f"FOUNDER_{index}" for index in range(len(people))]
    write_vcf(output / f"founders_chr{chromosome}.vcf", selected_keys, ids, values)
    np.savez_compressed(output / "observed.npz", calls=values)
    write_json(output / "founders.json", dict(
        chromosome=chromosome, people=people, founder_ids=ids, keys=selected_keys,
        rsids=[row[4] for row, keep in zip(rows, retained) if keep],
        source_loci=original_count, retained_loci=len(selected_keys),
        unavailable_exact_loci_by_input=unavailable,
        excluded_uncalled_loci=int((~retained).sum()),
        cohort_sha256=digest(cohort_path)))


def transmission(seed, chromosome, sample, donor):
    return hashlib.sha256(f"{seed}:{chromosome}:{sample}:{donor}".encode("ascii")).digest()[0] % 2


def make_cases(cohort_path, chromosome, folder, phased, masks, output):
    output.mkdir()
    cohort, data = read(cohort_path), read(folder / "founders.json")
    keys = [tuple(row) for row in data["keys"]]
    calls = read_call_matrix(phased, keys, data["founder_ids"], output, phased=True)
    observed = np.load(folder / "observed.npz")["calls"]
    if not np.array_equal(observed.sum(axis=1), calls.sum(axis=1)):
        raise ValueError("Production founder phasing changed observed diploid sequencing/array truth")
    index = {row["sample_id"]: i for i, row in enumerate(data["people"])}
    cases, arrays = [], []
    for person in cohort["people"]:
        donors = person["donors"]
        if person["kind"] == "pure":
            value = calls[index[donors[0]]].copy()
        elif person["kind"] == "F1" and len(donors) == 2 and donors[0] != donors[1]:
            value = np.stack([calls[index[donor], transmission(
                cohort["protocol"]["seed"], chromosome, person["sample_id"], donor)] for donor in donors])
        else:
            raise ValueError("Unsupported or biologically invalid founder transmission")
        arrays.append(value)
        cases.append(person)
    arrays = np.stack(arrays)
    write_vcf(output / f"truth_chr{chromosome}.vcf", keys, [row["sample_id"] for row in cases], arrays, phased=True)
    (output / "pieces").mkdir()
    jobs = []
    for mask in masks:
        typed = {tuple([fields[0], int(fields[1]), fields[2], fields[3]])
                 for line in Path(cohort["masks"][mask]).read_text().splitlines() if (fields := line.split())}
        use = [i for i, key in enumerate(keys) if key in typed]
        if len(use) < 100:
            raise ValueError(f"Too few truly observed {mask}/chr{chromosome} calls")
        for person_index, person in enumerate(cases):
            identity = person["production_ids"][mask]
            path = output / "pieces" / f"{identity}.txt"
            with path.open("w") as handle:
                for i in use:
                    chrom, pos, ref, alt = keys[i]
                    genotype = "".join((ref, alt)[int(allele)] for allele in arrays[person_index, :, i])
                    handle.write(f"{data['rsids'][i]}\t{chrom}\t{pos}\t{genotype}\n")
            jobs.append(dict(mask=mask, production_id=identity, sample_id=person["sample_id"], markers=len(use)))
    write_json(output / "cases.json", dict(chromosome=chromosome, jobs=jobs, keys=keys,
                                           founder_genotypes_verified=int(observed.shape[0] * observed.shape[2])))


def assemble(cohort_path, folders, masks, output):
    output.mkdir()
    cohort = read(cohort_path)
    ordered = sorted(folders, key=lambda folder: read(folder / "cases.json")["chromosome"])
    jobs = []
    for mask in masks:
        for person in cohort["people"]:
            identity = person["production_ids"][mask]
            target = output / f"{identity}.txt"
            with target.open("w") as handle:
                handle.write("# 23andMe synthetic benchmark; GRCh38 forward strand\n"
                             "# reference human assembly build 38\n"
                             "# rsid\tchromosome\tposition\tgenotype\n")
                for folder in ordered:
                    handle.write((folder / "pieces" / f"{identity}.txt").read_text())
            jobs.append(dict(mask=mask, production_id=identity, raw=target.name))
    write_json(output / "jobs.json", jobs)


def nucleotide_calls(path, regions, identity):
    result = {}
    for key, fields in genotype_rows(path, regions, [identity]):
        if key[1] in result:
            raise ValueError("Ambiguous harmonized position")
        genotype = fields[0]
        if not re.fullmatch(r"[01][|/][01]", genotype):
            raise ValueError("Missing or malformed observed genotype")
        alleles = (key[2], key[3])
        bases = "".join(sorted(alleles[int(char)] for char in (genotype[0], genotype[2])))
        if any(base not in "ACGT" for base in bases) or len(bases) != 2:
            raise ValueError("Invalid harmonized nucleotide genotype")
        result[key[1]] = bases
    return result


def validate_concordance(expected, actual, observed, policy):
    if expected.shape != actual.shape or observed.shape != expected.shape[1:]:
        raise ValueError("Genotype concordance axes differ")
    if not np.isin(expected, [0, 1, 2]).all() or not np.isin(actual, [0, 1, 2]).all():
        raise ValueError("Invalid or missing genotype dosage in concordance check")
    results = {}
    for label, selected in (("observed", observed), ("untyped", ~observed)):
        count = int(selected.sum())
        if label == "untyped" and count < policy["minimum_untyped_validation_markers"]:
            raise ValueError("Insufficient independent untyped genotype validation markers")
        if not count:
            raise ValueError("Empty observed genotype validation stratum")
        person_concordance = np.mean(expected[:, selected] == actual[:, selected], axis=1)
        concordance = float(person_concordance.mean())
        threshold = policy[f"minimum_{label}_genotype_concordance"]
        if person_concordance.min() < threshold:
            raise ValueError(f"{label} minimum-person genotype concordance {person_concordance.min():.6f} < {threshold}")
        results[label] = dict(markers=count, concordance=concordance,
                              minimum_person_concordance=float(person_concordance.min()))
    return results


def verify(cohort_path, chromosome, mask, truth_folder, merged, imputed, qc_files, output):
    output.mkdir()
    cohort, data = read(cohort_path), read(truth_folder / "cases.json")
    keys = [tuple(row) for row in data["keys"]]
    people = cohort["people"]
    expected_ids = [f"IMPID_{row['production_ids'][mask]}" for row in people]
    if set(sample_ids(merged)) != set(expected_ids) or set(sample_ids(imputed)) != set(expected_ids):
        raise ValueError("Production merge/imputation sample identity mismatch")
    truth = read_call_matrix(truth_folder / f"truth_chr{chromosome}.vcf.gz", keys,
                             [row["sample_id"] for row in people], output, phased=True)
    inferred = read_call_matrix(imputed, keys, expected_ids, output, phased=True)
    typed = {tuple([fields[0], int(fields[1]), fields[2], fields[3]])
             for line in Path(cohort["masks"][mask]).read_text().splitlines() if (fields := line.split())}
    observed = np.asarray([key in typed for key in keys])
    concordance = validate_concordance(truth.sum(axis=1), inferred.sum(axis=1), observed,
                                      cohort["protocol"]["production"])
    regions = positions(output / "observed_positions.tsv", [key for key, keep in zip(keys, observed) if keep])
    harmonized = []
    for index, identity in enumerate(expected_ids):
        original = {key[1]: "".join(sorted((key[2], key[3])[int(x)] for x in truth[index, :, i]))
                    for i, key in enumerate(keys) if observed[i]}
        actual = nucleotide_calls(merged, regions, identity)
        if set(actual) - set(original) or any(actual[pos] != original[pos] for pos in actual):
            raise ValueError("Production standardization/harmonization changed observed nucleotide genotypes")
        fraction = len(actual) / len(original)
        if fraction < cohort["protocol"]["production"]["minimum_harmonized_marker_retention"]:
            raise ValueError("Too many typed markers lost by production standardization/harmonization")
        harmonized.append(dict(sample_id=people[index]["sample_id"], retention=fraction,
                               retained=len(actual), excluded=len(original) - len(actual)))
    if len(qc_files) != len(people):
        raise ValueError("Require exactly one production split/QC result per person")
    qc_ids, qc_counts = [], {}
    for path in qc_files:
        ids = sample_ids(path)
        if len(ids) != 1 or ids[0] in qc_ids:
            raise ValueError("Split/QC sample identity duplicated or malformed")
        qc_ids += ids
        count = int(subprocess.check_output(["bcftools", "index", "-n", str(path)], text=True))
        if count < 100:
            raise ValueError("Production QC retained fewer than100 variants")
        expected_person = expected_ids.index(ids[0])
        key_indices = {key: index for index, key in enumerate(keys)}
        checked = 0
        for key, calls in genotype_rows(path, positions(output / "qc_positions.tsv", keys), ids):
            if key not in key_indices:
                continue
            call = calls[0]
            if not re.fullmatch(r"[01]\|[01]", call):
                raise ValueError("Production QC changed or removed valid haplotype phase")
            if [int(call[0]), int(call[2])] != inferred[expected_person, :, key_indices[key]].tolist():
                raise ValueError("Production split/QC changed imputed haplotypes")
            checked += 1
        if checked < 100:
            raise ValueError("Too few known genotype anchors survive production split/QC")
        qc_counts[ids[0]] = count
    if set(qc_ids) != set(expected_ids):
        raise ValueError("Split/QC output sample coverage incomplete")
    hashes = {str(path.resolve()): digest(path) for path in [merged, imputed, *qc_files]}
    write_json(output / "preprocessing_verified.json", dict(
        status="production_preprocessing_genotypes_and_coverage_verified", chromosome=chromosome,
        mask=mask, people=len(people), concordance=concordance, harmonization=harmonized,
        qc_retained=qc_counts, input_sha256=hashes, production_changed=False))
    print(f"{mask}/chr{chromosome}: real production merge/Beagle/split/QC passed; {concordance}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("rsids", "founders", "cases", "assemble", "verify"))
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--chromosomes", type=int, nargs="+")
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--markers", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--folders", type=Path, nargs="+")
    parser.add_argument("--phased", type=Path)
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
        make_cases(args.cohort, args.chromosome, args.folder, args.phased, args.masks, args.output)
    elif args.stage == "assemble":
        assemble(args.cohort, args.folders, args.masks, args.output)
    else:
        verify(args.cohort, args.chromosome, args.mask, args.folder, args.merged, args.imputed, args.qc, args.output)
