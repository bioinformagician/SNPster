#!/usr/bin/env python3
"""Check durable normalized outputs and gate exact-session successful-run cleanup."""
import argparse
import json
from pathlib import Path

from audit_production_ancestry import digest, write_json
from verify_sgdp_full_manifest import require_hashes, source_hashes


def smoke(run, output):
    record = json.loads((run / "verification/chr22.verified.json").read_text())
    panel = json.loads((run / "identity_panels/chr22/panel.json").read_text())
    calls = run / "normalized/chr22/normalized.vcf.gz"
    if (record["status"] != "coordinate_allele_and_genotype_preservation_verified"
            or record["chromosome"] != 22 or panel["chromosome"] != 22
            or record["sites"] < 1000 or panel["sites"] < 100
            or digest(calls) != record["normalized_sha256"]
            or digest(str(calls) + ".csi") != record["normalized_index_sha256"]
            or digest(calls) != panel["normalized_sha256"]
            or record["sample_ids"] != panel["candidate_ids"]
            or digest(run / "identity_panels/chr22/panel.vcf.gz") != panel["panel_sha256"]):
        raise ValueError("Real normalization/verification/identity-panel smoke failed")
    producers = source_hashes(run / "source_sha256.txt")
    require_hashes(producers)
    evidence = {str(path): digest(path) for folder in ("inputs", "normalized", "verification", "identity_panels")
                for path in (run / folder).rglob("*") if path.is_file()}
    write_json(output, dict(
        status="normalization_and_identity_panel_smoke_passed",
        chromosome=22, normalized_sites=record["sites"], identity_markers=panel["sites"],
        comparison_people=len(panel["sample_ids"]), source_sha256=producers,
        published_sha256=evidence, scientific_ancestry_launch_authorized=False,
    ))
    print("Real chr22 normalization, independent genotype verification and identity-panel smoke passed.", flush=True)


def require_smoke(path):
    record = json.loads(path.read_text())
    if record["status"] != "normalization_and_identity_panel_smoke_passed":
        raise ValueError("Successful real normalization/identity-panel smoke required")
    require_hashes(record["source_sha256"])
    require_hashes(record["published_sha256"])
    print("Frozen smoke sources and reusable evidence verified.", flush=True)


def verify_full(run, output):
    path = run / "summary/normalization_complete.json"
    record = json.loads(path.read_text())
    if (record["chromosomes"] != list(range(1, 23))
            or record["coordinate_harmonization_complete"] is not True
            or record["genotype_relatedness_screen_complete"] is not True
            or record["scientific_ancestry_launch_authorized"] is not False):
        raise ValueError("Complete normalized evidence and genetic identity screen required")
    require_hashes(source_hashes(run / "source_sha256.txt"))
    require_hashes(record["published_sha256"])
    write_json(output, dict(
        status="durable_normalized_data_and_genetic_audit_hashes_verified",
        normalization_manifest_sha256=digest(path), verified_files=len(record["published_sha256"]),
        scientific_ancestry_launch_authorized=False,
    ))
    print("All published normalized genotypes and genetic identity evidence verified.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("smoke", "require-smoke", "full"))
    parser.add_argument("--run", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.stage == "smoke":
        smoke(args.run, args.output)
    elif args.stage == "require-smoke":
        require_smoke(args.manifest)
    else:
        verify_full(args.run, args.output)
