#!/usr/bin/env python3
"""Bind the staging-only recovery sources without changing frozen producers."""
import argparse
import json
from pathlib import Path

from audit_production_ancestry import digest, write_json
from finalize_sgdp_normalization import smoke
from verify_sgdp_full_manifest import require_hashes, source_hashes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--attempt-sources", type=Path, required=True)
    args = parser.parse_args()
    output = args.run / "smoke_passed.json"
    smoke(args.run, output)
    evidence = source_hashes(args.attempt_sources)
    require_hashes(evidence)
    record = json.loads(output.read_text())
    record["source_sha256"].update(evidence)
    record["published_sha256"][str(args.attempt_sources)] = digest(args.attempt_sources)
    write_json(output, record)
