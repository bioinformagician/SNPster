#!/usr/bin/env python3
"""Check durable smoke provenance or collect all-autosome public-data QC."""
import argparse
import json
from pathlib import Path

from audit_production_ancestry import digest, write_json


def require_verified(smoke, project):
    manifest = json.loads((smoke / "published_verification.json").read_text())
    if (manifest["status"] != "published_genotypes_verified_not_ancestry_execution_ready"
            or manifest["chromosome"] != 22
            or manifest["scientific_ancestry_launch_authorized"] is not False):
        raise ValueError("A verified public chr22 sequencing smoke is required")
    for relative, expected in manifest["sha256"].items():
        if digest(smoke / relative) != expected:
            raise ValueError(f"Published smoke changed: {relative}")
    for line in (smoke / "source_sha256.txt").read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        if filename.endswith("nextflow-26.04.4-one.jar"):
            continue
        if digest(filename) != expected:
            raise ValueError(f"Frozen import producer/protocol changed: {filename}")
    protocol = json.loads((project / "config/sgdp_public_import.json").read_text())
    if manifest["source_counts"] != protocol["expected_metadata_disjoint_counts"]:
        raise ValueError("Public source coverage changed")
    write_json("smoke_gate.json", dict(
        status="verified_smoke_reusable_for_data_import_only",
        smoke=str(smoke), manifest_sha256=digest(smoke / "published_verification.json"),
        ancestry_execution_authorized=False,
    ))


def collect(smoke, reports, run):
    expected = set(range(1, 22))
    records = [json.loads(path.read_text()) for path in reports]
    if {item["chromosome"] for item in records} != expected or len(records) != len(expected):
        raise ValueError("Require exactly one public data QC result for chromosomes 1-21")
    reused = json.loads((smoke / "verification/verification.json").read_text())
    records.append(reused)
    if any(item["status"] != "raw_control_concordance_passed_not_ancestry_execution_ready"
           for item in records):
        raise ValueError("Incomplete or failed raw-call concordance")
    write_json("full_import.json", dict(
        status="producer_qc_complete_independent_all_autosome_verification_required",
        chromosomes=sorted(item["chromosome"] for item in records),
        reused_chr22=str(smoke),
        new_chromosome_root=str(run / "chromosomes"),
        raw_call_qc=sorted(records, key=lambda item: item["chromosome"]),
        public_metadata=str(smoke / "metadata/selected.json"),
        required_ancestry_resolutions=["eight_category", "seven_category"],
        independent_all_autosome_genotype_verification_complete=False,
        coordinate_harmonization_complete=False, genotype_relatedness_verified=False,
        scientific_ancestry_launch_authorized=False,
    ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("require", "collect"))
    parser.add_argument("--smoke", type=Path, required=True)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--reports", type=Path, nargs="+")
    args = parser.parse_args()
    if args.stage == "require":
        require_verified(args.smoke, args.project)
    else:
        collect(args.smoke, args.reports, args.run)
