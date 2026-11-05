#!/usr/bin/env python3
"""Audit unchanged production inputs before a matched ancestry benchmark."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

SOURCES = ("EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA")


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(1024 * 1024):
            checksum.update(block)
    return checksum.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def reference_audit(chromosome, production_root):
    reference = production_root / "beagle_references" / (
        f"hgdp1kgp_chr{chromosome}.filtered.SNV_INDEL.phased.shapeit5."
        "SNV_biallelic.numericCHR.vcf.gz"
    )
    index = Path(str(reference) + ".tbi")
    genetic_map = production_root / "plink.GRCh38.map" / f"plink.chr{chromosome}.GRCh38.map"
    identities = subprocess.check_output(
        ["bcftools", "query", "-l", str(reference)], text=True
    ).splitlines()
    if not identities or len(set(identities)) != len(identities):
        raise ValueError(f"Empty/duplicate production reference identities: chr{chromosome}")
    return {
        "chromosome": chromosome,
        "reference": str(reference),
        "sample_ids": identities,
        "files": {
            str(path): {"sha256": digest(path), "bytes": path.stat().st_size}
            for path in (reference, index, genetic_map)
        },
    }


def rows(path):
    with path.open(newline="") as handle:
        result = list(csv.DictReader(handle, delimiter="\t"))
    if not result:
        raise ValueError(f"Empty metadata: {path}")
    return result


def cohort_audit(records, reference_ids, authority):
    identities = [row["sample_id"] for row in records]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate cohort identities")
    eligible, overlap = [], []
    for row in records:
        sample = row["sample_id"]
        source = row.get("original_source_group", row["source_group"])
        if source not in SOURCES:
            raise ValueError(f"Invalid eight-source label: {sample}: {source}")
        if authority.get(sample) != (source, row["fine"]):
            raise ValueError(f"Authoritative source/fine label mismatch: {sample}")
        item = dict(sample_id=sample, source=source, fine=row["fine"])
        (overlap if sample in reference_ids else eligible).append(item)
    return dict(
        total=len(records), reference_overlap_count=len(overlap),
        reference_overlap=overlap, direct_id_disjoint_count=len(eligible),
        direct_id_disjoint_counts=dict(Counter(row["source"] for row in eligible)),
        direct_id_disjoint_people=eligible,
        alias_and_relatedness_audit_complete=False,
    )


def assemble(reference_paths, runtime_path, project, authority_path):
    audits = [json.loads(path.read_text()) for path in reference_paths]
    if sorted(item["chromosome"] for item in audits) != list(range(1, 23)):
        raise ValueError("Require exactly one production-reference audit per autosome")
    identities = audits[0]["sample_ids"]
    if any(item["sample_ids"] != identities for item in audits):
        raise ValueError("Production reference sample identity/order differs across chromosomes")
    reference_ids = set(identities)
    authority = {}
    for row in rows(authority_path):
        sample, label = row["sample"], (row["source_group"], row["fine"])
        if not label[0]:
            continue
        if sample in authority and authority[sample] != label:
            raise ValueError(f"Conflicting authority entries: {sample}")
        authority[sample] = label
    inputs = project / "results_v3_model_family_screen_05/inputs"
    cohorts = {
        "pure_training": inputs / "pure_train_labels.tsv",
        "pure_development": inputs / "pure.samples.tsv",
        "ancestry_reference": inputs / "reference_labels.tsv",
    }
    reports = {
        name: cohort_audit(rows(path), reference_ids, authority)
        for name, path in cohorts.items()
    }
    role_ids = {
        name: {row["sample_id"] for row in rows(path)}
        for name, path in cohorts.items()
    }
    for first, second in (
        ("pure_training", "pure_development"),
        ("pure_training", "ancestry_reference"),
        ("pure_development", "ancestry_reference"),
    ):
        if role_ids[first] & role_ids[second]:
            raise ValueError(f"Original ancestry role overlap: {first}/{second}")
    legacy_path = project / "results_v3_phasefree_array_validation_01/experiment.json"
    legacy = json.loads(legacy_path.read_text())
    donor_ids = legacy["donor_ids"]
    if not donor_ids or len(donor_ids) != len(set(donor_ids)):
        raise ValueError("Empty/duplicate pilot founder identities")
    counts = reports["pure_development"]["direct_id_disjoint_counts"]
    missing = [source for source in SOURCES if not counts.get(source)]
    if missing:
        status = "blocked_missing_independent_source_coverage"
    else:
        status = "not_ready_alias_relatedness_and_remaining_launch_gates_required"
    contract_path = project / "config/ancestry_benchmark_contract_v2.json"
    contract = json.loads(contract_path.read_text())
    if contract["version"] != 2 or contract["sources"]["prediction"] != list(SOURCES):
        raise ValueError("Unexpected production benchmark contract/source axis")
    for relative, expected in contract["production_preprocessing"]["source_sha256"].items():
        source = project.parents[2] / "data pipeline" / relative
        if digest(source) != expected:
            raise ValueError(f"Production source changed since contract freeze: {relative}")
    evidence = [*cohorts.values(), authority_path, legacy_path, contract_path, runtime_path]
    return dict(
        status=status, scientific_launch_authorized=False,
        production_reference_sample_count=len(identities),
        reference_chromosomes_verified=22,
        reported_resolutions=["eight_category", "seven_category"],
        metric="total_variation_error_pct",
        current_pilot_founders={
            "total": len(donor_ids),
            "reference_overlapping_ids": sorted(set(donor_ids) & reference_ids),
            "direct_id_disjoint_ids": sorted(set(donor_ids) - reference_ids),
        },
        cohorts=reports,
        missing_direct_id_disjoint_development_sources=missing,
        runtime=json.loads(runtime_path.read_text()),
        references=sorted(audits, key=lambda item: item["chromosome"]),
        input_sha256={str(path): digest(path) for path in evidence},
        interpretation=(
            "Direct-ID exclusion is necessary but not sufficient: aliases/relatedness and "
            "all remaining launch gates require review. No accuracy results were generated. "
            "Keep production preprocessing unchanged; do not silently replace this gate "
            "with an overlapping or reduced-source benchmark."
        ),
    )


def runtime_snapshot():
    import config

    jar, java = Path(config.BEAGLE_JAR), Path(config.JAVA_EXE)
    version = subprocess.run([str(java), "-version"], capture_output=True, text=True, check=True)
    help_result = subprocess.run(
        [str(java), "-jar", str(jar)], capture_output=True, text=True, check=True
    )
    return {
        "image_id": "sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1",
        "beagle_jar": str(jar), "beagle_sha256": digest(jar),
        "beagle_version": help_result.stdout.splitlines()[0],
        "java": str(java), "java_sha256": digest(java),
        "java_version": version.stderr.strip(),
        "production_source_sha256": {
            name: digest(Path("/app") / name)
            for name in ("main.py", "config.py", "imputation_classes.py")
        },
        "resolved_environment": {
            name: os.environ[name] for name in (
                "BEAGLE_JAR", "JAVA_EXE", "HEAP_GB", "THREADS",
                "PLINK_MAP_DIR", "BEAGLE_REFERENCE_DIR", "OUTPUT_DIR",
            )
        },
        "remaining_runtime_freeze": (
            "Cohort-specific batching/merge/split, harmonization and QC identities "
            "still require freezing before scientific execution."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("reference", "runtime", "report", "gate"))
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--production-root", type=Path)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--authority", type=Path)
    parser.add_argument("--runtime", type=Path)
    parser.add_argument("--references", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "reference":
        result = reference_audit(args.chromosome, args.production_root)
    elif args.stage == "runtime":
        result = runtime_snapshot()
    elif args.stage == "report":
        result = assemble(args.references, args.runtime, args.project, args.authority)
    else:
        report = json.loads(args.output.read_text())
        raise SystemExit(
            f"Scientific launch blocked: {report['status']}; missing direct-ID-disjoint "
            f"development sources: {report['missing_direct_id_disjoint_development_sources']}. "
            "Review published preflight.json; no ancestry model has been run."
        )
    write_json(args.output, result)
    print(f"{args.stage}: published {args.output}", flush=True)


if __name__ == "__main__":
    main()
