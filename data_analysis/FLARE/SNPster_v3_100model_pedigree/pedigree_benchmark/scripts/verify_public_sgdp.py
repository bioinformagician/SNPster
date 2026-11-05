#!/usr/bin/env python3
"""Independently rederive every published SGDP sequencing genotype before cleanup."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re
import subprocess

from audit_production_ancestry import digest, write_json


def expected_call(pl_text, ad_text, protocol):
    if "." in pl_text or "." in ad_text:
        return "./.", ".", "."
    pl = tuple(map(int, pl_text.split(",")))
    ad = tuple(map(int, ad_text.split(",")))
    if len(pl) != 3 or len(ad) != 2 or any(value < 0 for value in pl + ad):
        raise ValueError("Invalid published sequencing likelihoods/depths")
    best = pl.index(min(pl))
    gap = min(value for i, value in enumerate(pl) if i != best) - pl[best]
    depth = sum(ad)
    gt = "./."
    if (protocol["min_depth"] <= depth <= protocol["max_depth"]
            and gap >= protocol["min_likelihood_gap_phred"]):
        gt = ("0/0", "0/1", "1/1")[best]
    return gt, str(gap), str(depth)


def verify(run, project, chromosome=22, metadata_run=None, output=None):
    if chromosome not in range(1, 23):
        raise ValueError("Expected an autosome from 1 to 22")
    metadata_run = run if metadata_run is None else metadata_run
    protocol_path = project / "config/sgdp_public_import.json"
    protocol = json.loads(protocol_path.read_text())
    selection = json.loads((metadata_run / "metadata/selected.json").read_text())
    call_summary = json.loads((run / "basecalls/basecalls.json").read_text())
    control = json.loads((run / "verification/verification.json").read_text())
    source = json.loads((run / "source/download.json").read_text())
    public_metadata = metadata_run / "metadata/metadata.original.txt"
    if digest(public_metadata) != selection["source_metadata_sha256"]:
        raise ValueError("Published public metadata hash mismatch")
    with public_metadata.open(encoding=protocol["metadata_encoding"], newline="") as handle:
        records = {row["SGDP_ID"] + ".DG": row for row in csv.DictReader(handle, delimiter="\t")}
    preflight_path = project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
    if digest(preflight_path) != selection["production_preflight_sha256"]:
        raise ValueError("Production reference exclusion evidence changed")
    preflight = json.loads(preflight_path.read_text())
    forbidden = set(preflight["references"][0]["sample_ids"])
    for group in preflight["cohorts"].values():
        forbidden.update(row["sample_id"] for row in group["reference_overlap"] + group["direct_id_disjoint_people"])
    ids = [row["sample_id"] for row in selection["people"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate published donor IDs")
    for person in selection["people"]:
        row = records[person["sample_id"]]
        names = {
            item for field in ("Illumina_ID", "Sample_ID", "Sample_ID(Aliases)", "SGDP_ID")
            for item in re.split(r"[,;\s]+", row[field]) if item
        }
        if (row["Embargo"] != "FullyPublic" or row[list(row)[-1]] != "X"
                or names & forbidden or sorted(names) != person["aliases"]
                or protocol["candidate_population_sources"][row["Population_ID"]] != person["source"]):
            raise ValueError("Public access, ancestry-label or alias exclusion mismatch")
    if dict(Counter(row["source"] for row in selection["people"])) != protocol["expected_metadata_disjoint_counts"]:
        raise ValueError("Unexpected public donor source coverage")
    calls = run / "basecalls/basecalls.vcf.gz"
    expected_hash = digest(calls)
    if expected_hash != call_summary["output_sha256"] or expected_hash != control["basecalls_sha256"]:
        raise ValueError("Published basecall hash mismatch")
    if source["sha256"] != call_summary["source_sha256"]:
        raise ValueError("Source likelihood identity mismatch")
    raw = metadata_run / "raw_control/raw_control.vcf.gz"
    if digest(raw) != control["raw_control_sha256"]:
        raise ValueError("Original-call control hash mismatch")
    if (control["concordance"] < protocol["raw_control_min_concordance"]
            or control["n_compared"] < protocol["raw_control_min_shared_calls"]
            or control["scientific_ancestry_launch_authorized"] is not False
            or control["chromosome"] != chromosome or call_summary["chromosome"] != chromosome
            or control["selected_sha256"] != digest(metadata_run / "metadata/selected.json")
            or call_summary["genotypes_from"] != "mpileup_PL_and_AD_only"
            or call_summary["input_imputed_GT_DS_GP_used"] is not False):
        raise ValueError("Original-call validation or readiness declaration invalid")
    actual_ids = subprocess.check_output(["bcftools", "query", "-l", str(calls)], text=True).splitlines()
    if actual_ids != ids:
        raise ValueError("Published genotype sample axis differs")
    command = json.loads((run / "basecalls/likelihood_query.command.json").read_text())
    if command[command.index("-f") + 1] != "%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%AD]\n":
        raise ValueError("Reference-imputed fields used to construct observed genotypes")
    sites, cells, counts, previous_position, keys = 0, 0, Counter(), -1, set()
    process = subprocess.Popen([
        "bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%DP\t%PL\t%AD]\n",
        str(calls),
    ], stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            chrom, position, ref, alt = fields[:4]
            position = int(position)
            if (len(fields) != 4 + 5 * len(ids) or chrom != str(chromosome) or position < 1
                    or position < previous_position or len(ref) != 1 or len(alt) != 1
                    or ref not in "ACGT" or alt not in "ACGT" or ref == alt):
                raise ValueError("Invalid genotype locus/sample axis")
            if position != previous_position:
                keys.clear()
            if (ref, alt) in keys:
                raise ValueError("Duplicated exact allele")
            keys.add((ref, alt))
            previous_position = position
            for i, sample in enumerate(ids):
                gt, gq, dp, pl, ad = fields[4 + 5 * i:9 + 5 * i]
                if (gt, gq, dp) != expected_call(pl, ad, protocol):
                    raise ValueError(f"Genotype not reproducible from raw sequencing evidence: {sample}")
                counts[sample] += gt != "./."
                cells += 1
            sites += 1
        if process.wait():
            raise ValueError("Independent genotype query failed")
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait()
    if sites != call_summary["sites"] or dict(counts) != call_summary["called"]:
        raise ValueError("Published callability counts differ")
    if int(subprocess.check_output(["bcftools", "index", "-n", str(calls)], text=True)) != sites:
        raise ValueError("Published genotype index disagrees with streamed site count")
    if not sites or any(counts[sample] / sites < protocol["minimum_called_fraction"] for sample in ids):
        raise ValueError("Published callability fails declared threshold")
    paths = sorted(path for folder in ("source", "basecalls", "verification")
                   for path in (run / folder).iterdir() if path.is_file())
    shared = sorted(path for folder in ("metadata", "raw_control")
                    for path in (metadata_run / folder).iterdir() if path.is_file())
    if metadata_run == run:
        paths += shared
    result = dict(
        status="published_genotypes_verified_not_ancestry_execution_ready",
        chromosome=chromosome, sites=sites, genotype_cells_independently_rederived=cells,
        sample_ids=ids,
        source_counts=protocol["expected_metadata_disjoint_counts"],
        original_raw_control_concordance=control["concordance"],
        original_raw_control_compared=control["n_compared"],
        coordinate_harmonization_complete=False, genotype_relatedness_verified=False,
        scientific_ancestry_launch_authorized=False,
        sha256={str(path.relative_to(run)): digest(path) for path in paths},
        shared_sha256={str(path.relative_to(metadata_run)): digest(path) for path in shared},
    )
    write_json(run / "published_verification.json" if output is None else output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "sha256"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--metadata-run", type=Path)
    parser.add_argument("--chromosome", type=int, default=22)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    verify(args.run, Path(__file__).resolve().parents[2], args.chromosome, args.metadata_run, args.output)
