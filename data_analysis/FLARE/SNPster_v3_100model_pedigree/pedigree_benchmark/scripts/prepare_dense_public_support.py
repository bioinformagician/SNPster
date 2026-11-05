#!/usr/bin/env python3
"""Prepare dense observed donor candidates without choosing from ancestry results."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import re
import subprocess

from audit_production_ancestry import SOURCES, digest, rows, write_json
from import_sgdp_public import select_public
from recover_sgdp_duplicates import select_record
from sgdp_genetic_qc import genotype_rows, sample_ids
from verify_public_sgdp import expected_call
from verify_sgdp_full_manifest import require_hashes


def select_candidates(records, forbidden, protocol):
    people = select_public(records, forbidden, protocol)
    counts = Counter(person["source"] for person in people)
    if set(counts) != set(SOURCES) or min(counts.values()) < protocol["minimum_metadata_candidates_per_source"]:
        raise ValueError(f"Insufficient public metadata-disjoint eight-source coverage: {dict(counts)}")
    return people, dict(counts)


def metadata(project, protocol_path):
    protocol = json.loads(protocol_path.read_text())
    metadata_path = project / protocol["metadata"]
    preflight_path = project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
    preflight = json.loads(preflight_path.read_text())
    forbidden = set(preflight["references"][0]["sample_ids"])
    fitting = [project / "results_v3_model_family_screen_05/inputs" / name
               for name in ("reference_labels.tsv", "pure_train_labels.tsv")]
    for path in fitting:
        forbidden.update(row["sample_id"] for row in rows(path))
    previous_path = project / protocol["reusable_normalization"] / "summary/normalization_complete.json"
    cohort_path = project / "results_v3_production_cohort_20261022_110441/cohort/frozen/cohort.json"
    frozen = json.loads(cohort_path.read_text())
    require_hashes({str(previous_path): frozen["input_sha256"][str(previous_path)]})
    previous = json.loads(previous_path.read_text())
    known = set(previous["genetic_identity_screen"]["candidate_ids"])
    excluded = set()
    for pair in previous["genetic_identity_screen"]["candidate_involving_related_pairs"]:
        for candidate, partner in ((pair["first"], pair["second"]), (pair["second"], pair["first"])):
            if candidate in known and partner in forbidden:
                excluded.add(candidate)
    with metadata_path.open(encoding="latin-1") as handle:
        records = [row for row in csv.DictReader(handle, delimiter="\t")
                   if row["SGDP_ID"] + ".DG" not in excluded]
    people, counts = select_candidates(records, forbidden, protocol)
    reused = [person for person in people if person["sample_id"] in known]
    recalled = [person for person in people if person["sample_id"] not in known]
    control = next(person for person in people if person["sample_id"] == "S_Chane-1.DG")
    inputs = [metadata_path, preflight_path, protocol_path, previous_path, cohort_path, *fitting]
    write_json("recall.json", dict(
        people=sorted([*recalled, control], key=lambda person: person["sample_id"]),
        genotypes_from="mpileup_PL_AD_only", reused_control_not_a_new_individual=True))
    write_json("selected.json", dict(
        people=people, counts=counts, genotype_relatedness_verified=False,
        reused_ids=[person["sample_id"] for person in reused],
        new_ids=[person["sample_id"] for person in recalled],
        previously_excluded_fitting_relatives=sorted(excluded),
        donor_roles_frozen=False, ancestry_models_authorized=False,
        production_preprocessing_changed=False,
        legacy_development_identity_bridges_allowed=True,
        input_sha256={str(path): digest(path) for path in inputs}))
    print(f"Public dense candidate inventory: {len(people)} people; {counts}", flush=True)
    print(f"Reusing {len(reused)} verified people; recalling {len(recalled)} new people "
          "plus the existing original-call control.", flush=True)


def verify_raw(calls, summary_path, selection_path, protocol_path, output):
    output.mkdir()
    protocol = json.loads(protocol_path.read_text())
    summary = json.loads(summary_path.read_text())
    identities = [person["sample_id"] for person in json.loads(selection_path.read_text())["people"]]
    if digest(calls) != summary["output_sha256"] or sample_ids(calls) != identities:
        raise ValueError("Raw sequencing evidence hash or sample order changed")
    counts, original_sites, retained, removed, cells = Counter(), 0, 0, 0, 0
    evidence, pending, position = [], [], None

    def emit(destination):
        nonlocal retained, removed
        groups = {}
        for line in pending:
            fields = line.rstrip("\n").split("\t")
            key = tuple(fields[index] for index in (0, 1, 3, 4))
            groups.setdefault(key, []).append(line)
        for key, records in groups.items():
            chosen = select_record(records)
            if len(records) > 1:
                evidence.append(dict(key=key, input_records=records, retained_record=chosen))
            removed += len(records) - 1
            destination.write(chosen)
            retained += 1
            for identity, value in zip(identities, chosen.rstrip("\n").split("\t")[9:]):
                counts[identity] += value.split(":", 1)[0] != "./."

    target = output / "basecalls.vcf"
    with gzip.open(calls, "rt") as source, target.open("w") as destination:
        for line in source:
            if line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            current = int(fields[1])
            if (len(fields) != 9 + len(identities) or fields[8] != "GT:GQ:DP:PL:AD"
                    or fields[0] != str(summary["chromosome"]) or current < 1
                    or (position is not None and current < position)):
                raise ValueError("Invalid raw sequencing locus or evidence axes")
            if pending and current != position:
                emit(destination)
                pending = []
            position = current
            for value in fields[9:]:
                gt, gq, dp, pl, ad = value.split(":")
                if (gt, gq, dp) != expected_call(pl, ad, protocol):
                    raise ValueError("Independent PL/AD genotype rederivation differs")
                cells += 1
            pending.append(line)
            original_sites += 1
        if pending:
            emit(destination)
    if original_sites != summary["sites"] or original_sites != retained + removed:
        raise ValueError("Raw sequencing record accounting differs")
    fractions = {identity: counts[identity] / retained for identity in identities}
    if min(fractions.values()) < protocol["minimum_called_fraction"]:
        raise ValueError("Retained evidence fails unchanged callability policy")
    subprocess.run(["bgzip", str(target)], check=True)
    target = output / "basecalls.vcf.gz"
    subprocess.run(["bcftools", "index", str(target)], check=True)
    write_json(output / "duplicate_evidence.json", evidence)
    write_json(output / "verified.json", dict(
        status="every_observed_GT_GQ_DP_independently_rederived_from_PL_AD",
        chromosome=summary["chromosome"], original_sites=original_sites,
        retained_sites=retained, removed_strict_redundant_records=removed,
        genotype_cells_independently_rederived=cells, sample_ids=identities,
        called_fraction=fractions, input_sha256=digest(calls), output_sha256=digest(target),
        release_imputed_GT_DS_GP_used=False))


def verify_resources(project, source, target, chain, output):
    protocol = json.loads((project / "config/dense_public_support.json").read_text())
    paths = dict(source_fasta=source, source_index=Path(str(source) + ".fai"),
                 target_fasta=target, target_index=Path(str(target) + ".fai"), chain=chain)
    require_hashes({str(path): protocol["normalization_resource_sha256"][name]
                    for name, path in paths.items()})
    write_json(output, dict(status="original_normalization_resources_reverified",
                            resource_sha256=protocol["normalization_resource_sha256"],
                            ancestry_models_authorized=False))
    print("Original source/target FASTAs, indexes and chain identities verified.", flush=True)


def density(project, selection_path, normalized, verification_path, chromosome, control_path, output):
    output.mkdir()
    protocol = json.loads((project / "config/dense_public_support.json").read_text())
    selected = json.loads(selection_path.read_text())
    previous = json.loads((project / protocol["reusable_normalization"] /
                           "summary/normalization_complete.json").read_text())
    verification = json.loads(verification_path.read_text())
    control = json.loads(control_path.read_text())
    calling_policy = json.loads((project / "config/sgdp_public_import.json").read_text())
    if (verification["status"] != "coordinate_allele_and_genotype_preservation_verified"
            or verification["normalized_sha256"] != digest(normalized)
            or verification["source_vcf_sha256"] != control["basecalls_sha256"]
            or control["concordance"] < calling_policy["raw_control_min_concordance"]
            or control["n_compared"] < calling_policy["raw_control_min_shared_calls"]):
        raise ValueError("Independent normalization verification required before density audit")
    reused_path = project / protocol["reusable_normalization"] / f"normalized/chr{chromosome}/normalized.vcf.gz"
    require_hashes({str(path): previous["published_sha256"][str(path)]
                    for path in (reused_path, Path(str(reused_path) + ".csi"))})
    masks = {}
    for name in ("lonni", "frederik"):
        path = project / f"results_v3_phasefree_array_validation_01/{name}.observed_keys.tsv"
        masks[name] = {
            (fields[0], int(fields[1]), fields[2], fields[3])
            for line in path.read_text().splitlines()
            if (fields := line.split()) and fields[0] == str(chromosome)}
        if not masks[name]:
            raise ValueError(f"Empty real mask on chromosome {chromosome}: {name}")
    region = output / "mask_positions.tsv"
    region.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({
        key[1] for keys in masks.values() for key in keys})))
    counts = {person["sample_id"]: Counter() for person in selected["people"]}
    for path, identities in ((normalized, selected["new_ids"]),
                             (reused_path, selected["reused_ids"])):
        seen = set()
        for key, genotypes in genotype_rows(path, region, identities):
            if key in seen:
                raise ValueError("Ambiguous exact allele in mask-support audit")
            seen.add(key)
            for name, keys in masks.items():
                if key in keys:
                    for identity, genotype in zip(identities, genotypes):
                        counts[identity][name] += bool(re.fullmatch(r"[01]/[01]", genotype))
    people = [dict(
        person=person,
        mask_support={name: dict(called=counts[person["sample_id"]][name], total=len(keys),
                                fraction=counts[person["sample_id"]][name] / len(keys))
                      for name, keys in masks.items()})
        for person in selected["people"]]
    passing = [row for row in people if all(
        record["fraction"] >= protocol["minimum_real_mask_coverage"]
        for record in row["mask_support"].values())]
    source_counts = dict(Counter(row["person"]["source"] for row in passing))
    report = dict(
        status="observed_mask_support_measured_not_genetically_frozen",
        chromosome=chromosome, people=people, eligible_by_coverage=source_counts,
        minimum_real_mask_coverage=protocol["minimum_real_mask_coverage"],
        normalized_sha256=digest(normalized), reused_normalized_sha256=digest(reused_path),
        original_call_control_concordance=control["concordance"],
        independent_normalization_verified=True, ancestry_models_authorized=False)
    write_json(output / "density.json", report)
    print(f"chr{chromosome} >=90% real input support, both masks: {source_counts}", flush=True)
    if set(source_counts) != set(SOURCES) or min(source_counts.values()) < 2:
        raise ValueError(f"Dense support still lacks two eligible people per source: {source_counts}")


def finalize(run, post_cleanup):
    from verify_sgdp_full_manifest import source_hashes
    source_bindings = source_hashes(run.parent / "source_sha256.txt")
    require_hashes(source_bindings)
    manifest = run / "published_verification.json"
    require_hashes(json.loads((run / "metadata/selected.json").read_text())["input_sha256"])
    if post_cleanup:
        record = json.loads(manifest.read_text())
        require_hashes(record["published_sha256"])
        write_json(run / "post_cleanup_verification.json", dict(
            status="dense_observed_support_hashes_reverified_after_exact_session_cleanup",
            manifest_sha256=digest(manifest), ancestry_models_authorized=False))
        print("Dense observed evidence and source bindings verified after cache cleanup.", flush=True)
        return
    reports = sorted((run / "density").glob("chr*/density/density.json"))
    expected = [22] if run.name == "smoke" else list(range(1, 22))
    data = [json.loads(path.read_text()) for path in reports]
    if sorted(item["chromosome"] for item in data) != expected:
        raise ValueError("Dense source scope incomplete")
    for item in data:
        chromosome = item["chromosome"]
        normalized = run / f"normalized/chr{chromosome}/normalized.vcf.gz"
        observed = run / f"raw_observed/chr{chromosome}/observed/basecalls.vcf.gz"
        raw_report = json.loads((observed.parent / "verified.json").read_text())
        verified = json.loads((run / f"verification/chr{chromosome}.verified.json").read_text())
        if (digest(normalized) != item["normalized_sha256"]
                or verified["normalized_sha256"] != item["normalized_sha256"]
                or digest(str(normalized) + ".csi") != verified["normalized_index_sha256"]
                or digest(normalized.parent / "normalization.json") != verified["normalization_summary_sha256"]
                or digest(observed) != raw_report["output_sha256"]
                or verified["source_vcf_sha256"] != raw_report["output_sha256"]
                or raw_report["status"] != "every_observed_GT_GQ_DP_independently_rederived_from_PL_AD"):
            raise ValueError("Published dense observed/normalized copies differ from independent checks")
    if run.name != "smoke":
        smoke_manifest = json.loads((run.parent / "smoke/published_verification.json").read_text())
        require_hashes(smoke_manifest["published_sha256"])
        if digest(run / "metadata/selected.json") != digest(run.parent / "smoke/metadata/selected.json"):
            raise ValueError("Dense candidate selection changed between smoke and all-autosome recovery")
        data += [json.loads((run.parent / "smoke/density/chr22/density/density.json").read_text())]
    eligible = None
    for item in data:
        current = {row["person"]["sample_id"] for row in item["people"]
                   if all(record["fraction"] >= item["minimum_real_mask_coverage"]
                          for record in row["mask_support"].values())}
        eligible = current if eligible is None else eligible & current
    labels = {row["person"]["sample_id"]: row["person"]["source"] for row in data[0]["people"]}
    coverage_counts = dict(Counter(labels[identity] for identity in eligible))
    if set(coverage_counts) != set(SOURCES) or min(coverage_counts.values()) < 2:
        write_json(run / "coverage_blocked.json", dict(eligible_source_counts=coverage_counts))
        raise ValueError("No two-per-source cohort meets real input coverage across this entire scope")
    published = {} if run.name == "smoke" else {
        **smoke_manifest["published_sha256"],
        str(run.parent / "smoke/published_verification.json"):
            digest(run.parent / "smoke/published_verification.json")}
    for folder in ("resources", "metadata", "raw_observed", "normalized", "verification", "raw_control", "density"):
        for path in (run / folder).rglob("*"):
            if path.is_file():
                published[str(path.resolve())] = digest(path)
    write_json(manifest, dict(
        status="dense_observed_support_verified_not_genetically_frozen",
        chromosomes=[22] if run.name == "smoke" else list(range(1, 23)),
        new_chromosomes=expected, published_sha256=published,
        eligible_by_real_mask_coverage=sorted(eligible), eligible_source_counts=coverage_counts,
        source_sha256=source_bindings, ancestry_models_authorized=False))
    require_hashes(published)
    print(f"Verified durable observed support for {expected}; no ancestry models authorized.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("metadata", "resources", "verify-raw", "density", "finalize"))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--target", type=Path)
    parser.add_argument("--chain", type=Path)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--calls", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--normalized", type=Path)
    parser.add_argument("--verification", type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--raw-control", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--post-cleanup", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.stage == "resources":
        verify_resources(args.project, args.source, args.target, args.chain, args.output)
    elif args.stage == "metadata":
        metadata(args.project, args.protocol)
    elif args.stage == "verify-raw":
        verify_raw(args.calls, args.summary, args.selection, args.protocol, args.output)
    elif args.stage == "density":
        density(args.project, args.selection, args.normalized, args.verification,
                args.chromosome, args.raw_control, args.output)
    else:
        finalize(args.run, args.post_cleanup)
