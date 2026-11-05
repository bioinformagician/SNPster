#!/usr/bin/env python3
"""Bind independent chromosome checks to durable all-autosome import outputs."""
import argparse
import json
from pathlib import Path

from audit_production_ancestry import digest, write_json


def require_hashes(hashes):
    for filename, expected in hashes.items():
        if digest(filename) != expected:
            raise ValueError(f"Durable evidence changed: {filename}")


def source_hashes(path):
    hashes = {}
    for line in path.read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        if not filename.endswith("nextflow-26.04.4-one.jar"):
            hashes[filename] = expected
    return hashes


def collect(import_run, smoke, verification_run, reports, output):
    producer = json.loads((import_run / "summary/full_import.json").read_text())
    if (producer["chromosomes"] != list(range(1, 23))
            or producer["reused_chr22"] != str(smoke)
            or producer["new_chromosome_root"] != str(import_run / "chromosomes")
            or producer["scientific_ancestry_launch_authorized"] is not False):
        raise ValueError("Unexpected full-import coverage, paths or authorization")
    records = [json.loads(path.read_text()) for path in reports]
    if (len(records) != 22
            or sorted(item["chromosome"] for item in records) != list(range(1, 23))):
        raise ValueError("Exactly one independent report per autosome is required")
    expected_ids = [row["sample_id"] for row in
                    json.loads((smoke / "metadata/selected.json").read_text())["people"]]
    producer_qc = {row["chromosome"]: row for row in producer["raw_call_qc"]}
    hashes = source_hashes(import_run / "source_sha256.txt")
    hashes.update(source_hashes(verification_run / "source_sha256.txt"))
    hashes[str(import_run / "source_sha256.txt")] = digest(import_run / "source_sha256.txt")
    hashes[str(import_run / "summary/full_import.json")] = digest(import_run / "summary/full_import.json")
    for record in records:
        chromosome = record["chromosome"]
        if (record["status"] != "published_genotypes_verified_not_ancestry_execution_ready"
                or record["scientific_ancestry_launch_authorized"] is not False
                or record["sample_ids"] != expected_ids
                or record["sites"] < 1
                or record["genotype_cells_independently_rederived"] != record["sites"] * len(expected_ids)
                or record["original_raw_control_concordance"] != producer_qc[chromosome]["concordance"]
                or record["original_raw_control_compared"] != producer_qc[chromosome]["n_compared"]):
            raise ValueError(f"Independent report is incomplete or disagrees: chr{chromosome}")
        root = smoke if chromosome == 22 else import_run / f"chromosomes/chr{chromosome}"
        for base, manifest in ((root, record["sha256"]), (smoke, record["shared_sha256"])):
            for relative, expected in manifest.items():
                filename = (base / relative).resolve()
                if not filename.is_relative_to(base.resolve()):
                    raise ValueError("Evidence path escapes published directory")
                if str(filename) in hashes and hashes[str(filename)] != expected:
                    raise ValueError(f"Conflicting evidence hashes: {filename}")
                hashes[str(filename)] = expected
    require_hashes(hashes)
    result = dict(
        status="all_autosome_public_genotypes_independently_verified_not_ancestry_ready",
        import_run=str(import_run), smoke_run=str(smoke),
        chromosomes=list(range(1, 23)), sample_ids=expected_ids,
        sites=sum(item["sites"] for item in records),
        genotype_cells_independently_rederived=sum(
            item["genotype_cells_independently_rederived"] for item in records),
        original_raw_control_compared=sum(item["original_raw_control_compared"] for item in records),
        minimum_chromosome_raw_concordance=min(item["original_raw_control_concordance"] for item in records),
        chromosome_reports=sorted(records, key=lambda record: record["chromosome"]),
        coordinate_harmonization_complete=False, genotype_relatedness_verified=False,
        scientific_ancestry_launch_authorized=False, inputs_sha256=hashes,
    )
    write_json(output, result)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("inputs_sha256", "chromosome_reports")}, indent=2))


def recheck(manifest, output):
    record = json.loads(manifest.read_text())
    if record["status"] != "all_autosome_public_genotypes_independently_verified_not_ancestry_ready":
        raise ValueError("An independent all-autosome verification is required")
    require_hashes(record["inputs_sha256"])
    write_json(output, dict(
        status="durable_evidence_rechecked_after_exact_import_cache_cleanup",
        verification_sha256=digest(manifest), verified_files=len(record["inputs_sha256"]),
        scientific_ancestry_launch_authorized=False,
    ))
    print("Durable all-autosome evidence rechecked after exact-session cleanup.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("collect", "recheck"))
    parser.add_argument("--import-run", type=Path)
    parser.add_argument("--smoke", type=Path)
    parser.add_argument("--verification-run", type=Path)
    parser.add_argument("--reports", type=Path, nargs="+")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "collect":
        collect(args.import_run, args.smoke, args.verification_run, args.reports, args.output)
    else:
        recheck(args.manifest, args.output)
