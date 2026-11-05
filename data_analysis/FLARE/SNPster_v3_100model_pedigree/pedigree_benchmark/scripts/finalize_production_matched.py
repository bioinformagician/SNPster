#!/usr/bin/env python3
"""Bind independently checked task results to durable copies before exact cleanup."""
import argparse
import json
from pathlib import Path

from audit_production_ancestry import digest, write_json
from verify_sgdp_full_manifest import require_hashes, source_hashes


def verify(run, post_cleanup=False):
    manifest_path = run / "published_verification.json"
    if post_cleanup:
        manifest = json.loads(manifest_path.read_text())
        require_hashes(manifest["published_sha256"])
        require_hashes(manifest["source_sha256"])
        write_json(run / "post_cleanup_verification.json", dict(
            status="all_published_benchmark_evidence_verified_after_cleanup",
            files_verified=len(manifest["published_sha256"]),
            manifest_sha256=digest(manifest_path)))
        print("Post-cleanup benchmark evidence reverified.", flush=True)
        return
    result = run / "summary/summary/independent_verification.json"
    record = json.loads(result.read_text())
    if record["status"] not in (
            "complete_verified_development_comparison", "smoke_passed_not_a_whole_genome_leaderboard"):
        raise ValueError("Independently verified complete scope required")
    require_hashes(record["independently_verified_model_sha256"])
    sources = source_hashes(run.parent / "source_sha256.txt")
    require_hashes(sources)
    directories = ("preparation", "founders", "founder_phase", "truth", "raw", "standardized", "combined",
                   "harmonized", "merged", "imputed", "split", "qc", "indexed_qc", "preprocessing",
                   "preprocessing_gate", "panels", "models", "summary")
    if record.get("truth_design") == "synthetic_reference_haplotype_development":
        directories = tuple(name for name in directories if name != "founder_phase")
    hashes = {}
    for name in directories:
        folder = run / name
        if not folder.is_dir():
            raise ValueError(f"Missing durable benchmark output: {name}")
        for path in folder.rglob("*"):
            if path.is_file():
                hashes[str(path)] = digest(path)
    write_json(manifest_path, dict(
        status=record["status"], source_sha256=sources, published_sha256=hashes,
        independent_verification_sha256=digest(result),
        production_changed=False, production_promotion=False, independent_final_test=False))
    require_hashes(hashes)
    print(f"All {len(hashes)} published benchmark evidence files independently bound and verified.", flush=True)


def verify_references(run, root, cohort_path, chromosomes):
    cohort = json.loads(cohort_path.read_text())
    hashes = {}
    for chromosome in chromosomes:
        folder = root / "beagle_references"
        metadata_path = folder / f"chr{chromosome}.reference.json"
        metadata = json.loads(metadata_path.read_text())
        filename = f"hgdp1kgp_chr{chromosome}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"
        original = Path("/srv/dependencies/imputation_runner/imputer/beagle_references") / filename
        if (metadata["status"] != "all_retained_reference_alleles_and_phased_genotypes_preserved"
                or metadata["cohort_sha256"] != digest(cohort_path)
                or metadata["retained_ids"] != cohort["imputation_reference_ids"]
                or metadata["excluded_ids"] != cohort["imputation_excluded"]
                or metadata["source_sha256"] != cohort["template_sha256"][str(original)]):
            raise ValueError("Published benchmark reference differs from the frozen purge/genotype audit")
        hashes[str(folder / filename)] = metadata["copy_sha256"]
        hashes[str(folder / (filename + ".csi"))] = metadata["index_sha256"]
        hashes[str(metadata_path)] = digest(metadata_path)
    require_hashes(hashes)
    sources = source_hashes(run.parent / "source_sha256.txt")
    require_hashes(sources)
    write_json(run / "published_verification.json", dict(
        status="donor_purged_benchmark_references_and_all_retained_genotypes_verified",
        chromosomes=chromosomes, published_sha256=hashes, source_sha256=sources,
        production_resources_modified=False))
    print(f"Durable benchmark reference copies verified for {chromosomes}.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--post-cleanup", action="store_true")
    parser.add_argument("--references", type=Path)
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--chromosomes", nargs="+", type=int)
    args = parser.parse_args()
    if args.references is not None and not args.post_cleanup:
        verify_references(args.run, args.references, args.cohort, args.chromosomes)
    else:
        verify(args.run, args.post_cleanup)
