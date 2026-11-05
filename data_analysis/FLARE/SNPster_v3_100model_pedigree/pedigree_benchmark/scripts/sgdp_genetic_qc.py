#!/usr/bin/env python3
"""Construct exact-allele identity panels and audit every public candidate pair."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import subprocess

from audit_production_ancestry import digest, write_json


def run(command, label):
    write_json(str(label) + ".command.json", list(map(str, command)))
    with Path(str(label) + ".log").open("w") as handle:
        subprocess.run(list(map(str, command)), check=True, stdout=handle, stderr=subprocess.STDOUT)


def sample_ids(path):
    ids = subprocess.check_output(["bcftools", "query", "-l", str(path)], text=True).splitlines()
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Empty or duplicated genotype sample IDs")
    return ids


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def query_keys(path, chromosome, regions=None):
    command = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n"]
    command += ["-R", str(regions)] if regions else ["-r", str(chromosome)]
    text = subprocess.check_output(command + [str(path)], text=True)
    result = []
    seen = set()
    for line in text.splitlines():
        chrom, pos, ref, alt = line.split("\t")
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            continue
        key = (chrom, int(pos), ref, alt)
        if key in seen:
            raise ValueError("Duplicate exact allele in identity-panel input")
        seen.add(key)
        result.append(key)
    return result


def genotype_rows(path, regions, ids):
    command = ["bcftools", "query", "-R", str(regions), "-s", ",".join(ids),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4 + len(ids):
                raise ValueError("Identity-panel sample axis mismatch")
            key = (fields[0], int(fields[1]), fields[2], fields[3])
            yield key, fields[4:]
        if process.wait():
            raise ValueError("Identity-panel genotype query failed")
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait()


def panel(project, chromosome, normalized, verification, production_root, output):
    output.mkdir()
    verified = json.loads(verification.read_text())
    if (verified["status"] != "coordinate_allele_and_genotype_preservation_verified"
            or verified["chromosome"] != chromosome
            or digest(normalized) != verified["normalized_sha256"]):
        raise ValueError("Independently verified normalized genotypes required")
    preflight = json.loads((project / "results_v3_production_preflight_20261021_193442/summary/preflight.json").read_text())
    reference_report = next(row for row in preflight["references"] if row["chromosome"] == chromosome)
    production = production_root / "beagle_references" / Path(reference_report["reference"]).name
    for filename, item in reference_report["files"].items():
        if digest(filename) != item["sha256"]:
            raise ValueError(f"Fixed production reference/map changed: {filename}")
    fixed = sample_ids(production)
    if fixed != reference_report["sample_ids"]:
        raise ValueError("Fixed imputation reference sample order changed")
    inputs = project / "results_v3_model_family_screen_05/inputs"
    specifications = [
        ("ancestry_reference", "reference.retained.vcf.gz"),
        ("pure_training", "pure_train.canonical.vcf.gz"),
        ("pure_development", "pure_pure.canonical.vcf.gz"),
    ]
    groups, union = [], set(fixed)
    for role, filename in specifications:
        path = inputs / filename
        ids = sample_ids(path)
        expected = {row["sample_id"] for row in
                    preflight["cohorts"][role]["reference_overlap"] +
                    preflight["cohorts"][role]["direct_id_disjoint_people"]}
        if set(ids) != expected:
            raise ValueError(f"Ancestry fitting/development genotype identities changed: {role}")
        additional = [identity for identity in ids if identity not in union]
        union.update(ids)
        groups.append((path, additional))
    candidates = sample_ids(normalized)
    if candidates != verified["sample_ids"] or union & set(candidates):
        raise ValueError("Candidate genotype axis or reference exclusion mismatch")
    keys = query_keys(inputs / specifications[0][1], chromosome)
    positions = output / "positions.tsv"
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({key[1] for key in keys})))
    shared = set(keys) & set(query_keys(production, chromosome, positions))
    shared &= set(query_keys(normalized, chromosome, positions))
    for path, _ in groups:
        shared &= set(query_keys(path, chromosome, positions))
    if not shared or len(shared) != len({key[1] for key in shared}):
        raise ValueError("Empty or ambiguous shared identity-marker panel")
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({key[1] for key in shared})))
    columns = []
    for path, identities in groups + [(normalized, candidates)]:
        if identities:
            mapping = {key: values for key, values in genotype_rows(path, positions, identities) if key in shared}
            if set(mapping) != shared:
                raise ValueError("Incomplete exact-allele genotype block")
            columns.append((identities, mapping))
    ids = fixed + [identity for identities, _ in columns for identity in identities]
    if len(ids) != len(set(ids)) or set(ids) != union | set(candidates):
        raise ValueError("Incomplete or duplicated identity-panel comparison universe")
    called, retained, seen = Counter(), 0, set()
    path = output / "panel.vcf"
    with path.open("w") as handle:
        handle.write(f'##fileformat=VCFv4.2\n##contig=<ID={chromosome}>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased observed genotype">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + "\t".join(ids) + "\n")
        for key, fixed_values in genotype_rows(production, positions, fixed):
            if key not in shared:
                continue
            if key in seen:
                raise ValueError("Duplicate fixed-reference genotype row")
            seen.add(key)
            values = fixed_values + [value for _, mapping in columns for value in mapping[key]]
            values = [value.replace("|", "/") for value in values]
            if any(value not in ("0/0", "0/1", "1/0", "1/1", "./.") for value in values):
                raise ValueError("Invalid biallelic diploid identity genotype")
            chrom, pos, ref, alt = key
            handle.write(f"{chrom}\t{pos}\t{chrom}:{pos}:{ref}:{alt}\t{ref}\t{alt}\t.\t.\t.\tGT\t" +
                         "\t".join(values) + "\n")
            for identity, value in zip(ids, values):
                called[identity] += value != "./."
            retained += 1
    if seen != shared:
        raise ValueError("Shared identity-panel coverage differs")
    subprocess.run(["bgzip", str(path)], check=True)
    subprocess.run(["bcftools", "index", str(path) + ".gz"], check=True)
    write_json(output / "panel.json", dict(
        chromosome=chromosome, sites=retained, sample_ids=ids, candidate_ids=candidates,
        called=dict(called), normalized_sha256=digest(normalized),
        production_reference_sha256=digest(production),
        original_genotype_sha256={str(path): digest(path) for path, _ in groups},
        panel_sha256=digest(str(path) + ".gz"), scientific_ancestry_launch_authorized=False,
    ))
    print(f"chr{chromosome}: exact identity panel has {retained:,} markers and {len(ids):,} people.", flush=True)


def audit_pairs(records, ids, candidates, policy):
    candidate_set = set(candidates)
    expected = len(candidates) * (len(ids) - len(candidates)) + len(candidates) * (len(candidates) - 1) // 2
    known, seen, related, excluded, minimum = set(ids), set(), [], set(), None
    for row in records:
        first, second = row["IID1"], row["IID2"]
        key = tuple(sorted((first, second)))
        if first == second or not {first, second} <= known or not {first, second} & candidate_set or key in seen:
            raise ValueError("Invalid, irrelevant or duplicate KING comparison pair")
        seen.add(key)
        nsnp = int(row["NSNP"])
        coefficient = float(row["KINSHIP"])
        ibs0 = float(row["IBS0"])
        discordance = ibs0 + float(row["HET1_HOM2"]) + float(row["HET2_HOM1"])
        if not math.isfinite(coefficient) or not math.isfinite(discordance) or not 0 <= discordance <= 1:
            raise ValueError("Nonfinite or invalid KING estimates")
        if nsnp < policy["minimum_pair_called_markers"]:
            raise ValueError(f"Insufficient independently called markers for {first}/{second}: {nsnp}")
        minimum = nsnp if minimum is None else min(minimum, nsnp)
        if coefficient >= policy["king_kinship_exclusion_threshold"] or 1 - discordance >= policy["duplicate_concordance_threshold"]:
            related.append(dict(first=first, second=second, kinship=coefficient,
                                genotype_concordance=1 - discordance, called_markers=nsnp,
                                both_public_candidates=first in candidate_set and second in candidate_set))
            if not ({first, second} <= candidate_set):
                excluded.update({first, second} & candidate_set)
    if len(seen) != expected:
        raise ValueError(f"Incomplete candidate-relatedness pair coverage: {len(seen)}/{expected}")
    return dict(
        pairs_checked=len(seen), minimum_pair_called_markers=minimum,
        excluded_candidate_ids=sorted(excluded), candidate_involving_related_pairs=related,
        eligible_candidate_ids=sorted(candidate_set - excluded),
    )


def kinship(project, panels, output):
    output.mkdir()
    reports = [json.loads((folder / "panel.json").read_text()) for folder in panels]
    if sorted(row["chromosome"] for row in reports) != list(range(1, 23)):
        raise ValueError("All22 identity-marker panels required")
    ids, candidates = reports[0]["sample_ids"], reports[0]["candidate_ids"]
    if any(row["sample_ids"] != ids or row["candidate_ids"] != candidates for row in reports):
        raise ValueError("Identity sample axis varies between chromosomes")
    policy = json.loads((project / "config/sgdp_normalization.json").read_text())["relatedness"]
    sites = sum(row["sites"] for row in reports)
    totals = {identity: sum(row["called"][identity] for row in reports) for identity in ids}
    if any(totals[identity] / sites < policy["minimum_sample_called_fraction"] for identity in ids):
        raise ValueError("Identity-panel sample callability below predeclared threshold")
    ordered = sorted(panels, key=lambda folder: json.loads((folder / "panel.json").read_text())["chromosome"])
    concat = ["bcftools", "concat", "-Oz", "-o", str(output / "all.vcf.gz")]
    run(concat + [folder / "panel.vcf.gz" for folder in ordered], output / "concat")
    run(["plink2", "--vcf", output / "all.vcf.gz", "--const-fid", "0", "--threads", "1",
         "--memory", "4096", "--make-pgen", "--out", output / "all"], output / "import")
    run(["plink2", "--pfile", output / "all", "--maf", str(policy["maf_min"]),
         "--geno", str(policy["marker_missing_max"]), "--indep-pairwise",
         f'{policy["ld_window_kb"]}kb', str(policy["ld_r2_max"]), "--threads", "1",
         "--memory", "4096", "--out", output / "prune"], output / "pruning")
    n_markers = len((output / "prune.prune.in").read_text().splitlines())
    if n_markers < policy["minimum_pruned_markers"]:
        raise ValueError(f"Insufficient LD-pruned identity markers: {n_markers}")
    candidate_file = output / "candidates.id"
    candidate_file.write_text("".join(f"0\t{identity}\n" for identity in candidates))
    run(["plink2", "--pfile", output / "all", "--extract", output / "prune.prune.in",
         "--make-king-table", "cols=id,nsnp,hethet,ibs0,ibs1,kinship",
         "--king-table-require", candidate_file, "--threads", "1", "--memory", "4096",
         "--out", output / "relatedness"], output / "king")
    table = output / "relatedness.kin0"
    with table.open() as handle:
        names = handle.readline().lstrip("#").split()
        pairs = [dict(zip(names, line.split())) for line in handle if line.strip()]
    result = audit_pairs(pairs, ids, candidates, policy)
    result.update(
        status="candidate_identity_and_king_screen_complete_not_source_purity_certification",
        all_people=len(ids), all22_panel_sites=sites, pruned_markers=n_markers,
        candidate_ids=candidates, sample_called_fractions={identity: totals[identity] / sites for identity in ids},
        genotype_relatedness_screen_complete=True, scientific_ancestry_launch_authorized=False,
    )
    write_json(output / "genetic_eligibility.json", result)
    print(json.dumps({key: value for key, value in result.items()
                      if key != "sample_called_fractions"}, indent=2), flush=True)


def collect(project, run_dir, folders, reports, panels, kinship_path, output):
    records = [json.loads(path.read_text()) for path in reports]
    if sorted(row["chromosome"] for row in records) != list(range(1, 23)) or len(records) != 22:
        raise ValueError("Exactly one independently normalized report per autosome required")
    eligible = json.loads((kinship_path / "genetic_eligibility.json").read_text())
    roots, hashes = {}, {}
    for folder in folders:
        summary = json.loads((folder / "normalization.json").read_text())
        chromosome = summary["chromosome"]
        record = next(row for row in records if row["chromosome"] == chromosome)
        if digest(folder / "normalized.vcf.gz") != record["normalized_sha256"]:
            raise ValueError("Normalized published evidence differs")
        destination = run_dir / f"normalized/chr{chromosome}"
        roots[str(chromosome)] = str(destination / "normalized.vcf.gz")
        for path in folder.iterdir():
            if path.is_file():
                hashes[str(destination / path.name)] = digest(path)
    for path in reports:
        hashes[str(run_dir / "verification" / path.name)] = digest(path)
    for folder in panels:
        chromosome = json.loads((folder / "panel.json").read_text())["chromosome"]
        hashes[str(run_dir / f"identity_panels/chr{chromosome}/panel.json")] = digest(folder / "panel.json")
    for path in kinship_path.iterdir():
        if path.is_file():
            hashes[str(run_dir / "kinship" / path.name)] = digest(path)
    hashes[str(run_dir / "source_sha256.txt")] = digest(run_dir / "source_sha256.txt")
    result = dict(
        status="normalization_and_genetic_identity_screen_complete_benchmark_gates_remain",
        chromosomes=list(range(1, 23)), sample_ids=records[0]["sample_ids"],
        normalized_vcfs=roots, sites=sum(row["sites"] for row in records),
        genotype_cells_verified=sum(row["genotype_cells_verified"] for row in records),
        coordinate_harmonization_complete=True, genetic_identity_screen=eligible,
        genotype_relatedness_screen_complete=True, source_purity_certified=False,
        production_preprocessing_changed=False, scientific_ancestry_launch_authorized=False,
        required_ancestry_resolutions=["eight_category", "seven_category"],
        remaining_gates=json.loads((project / "config/sgdp_normalization.json").read_text())["remaining_gates"],
        published_sha256=hashes,
    )
    write_json(output, result)
    print(f"All22 normalized chromosomes and {eligible['pairs_checked']:,} identity pairs verified.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("panel", "kinship", "collect"))
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--normalized", type=Path)
    parser.add_argument("--verification", type=Path)
    parser.add_argument("--production-root", type=Path)
    parser.add_argument("--panels", type=Path, nargs="+")
    parser.add_argument("--folders", type=Path, nargs="+")
    parser.add_argument("--reports", type=Path, nargs="+")
    parser.add_argument("--kinship", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "panel":
        panel(args.project, args.chromosome, args.normalized, args.verification, args.production_root, args.output)
    elif args.stage == "kinship":
        kinship(args.project, args.panels, args.output)
    else:
        collect(args.project, args.run_dir, args.folders, args.reports, args.panels, args.kinship, args.output)
