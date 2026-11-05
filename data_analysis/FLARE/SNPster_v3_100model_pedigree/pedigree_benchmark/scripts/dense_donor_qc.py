#!/usr/bin/env python3
"""Screen dense observed donors without using ancestry predictions or assigning roles."""
import argparse
from collections import Counter
import copy
import itertools
import json
from pathlib import Path
import subprocess

from audit_production_ancestry import SOURCES, digest, write_json
from sgdp_genetic_qc import audit_pairs, genotype_rows, query_keys, sample_ids
from verify_sgdp_full_manifest import require_hashes, source_hashes


def prepare(project, dense_run):
    dense_manifest = dense_run / "full/published_verification.json"
    dense = json.loads(dense_manifest.read_text())
    post = json.loads((dense_run / "full/post_cleanup_verification.json").read_text())
    if (dense["chromosomes"] != list(range(1, 23))
            or dense["ancestry_models_authorized"] is not False
            or post["manifest_sha256"] != digest(dense_manifest)):
        raise ValueError("Completed independently verified all-autosome dense inputs required")
    require_hashes(dense["published_sha256"])
    require_hashes(dense["source_sha256"])
    selected_path = dense_run / "full/metadata/selected.json"
    selected = json.loads(selected_path.read_text())
    require_hashes(selected["input_sha256"])
    protocol = json.loads((project / "config/dense_public_support.json").read_text())
    previous_root = project / protocol["reusable_normalization"]
    previous_path = previous_root / "summary/normalization_complete.json"
    previous = json.loads(previous_path.read_text())
    original_selection = project / "results_v3_public_sgdp_smoke_20261021_203540/metadata/selected.json"
    original = json.loads(original_selection.read_text())
    historical = previous["sample_ids"]
    new_ids = selected["new_ids"]
    public_ids = new_ids + historical
    candidates = [person["sample_id"] for person in selected["people"]]
    if (len(public_ids) != len(set(public_ids))
            or not set(candidates) <= set(public_ids)
            or set(historical) != set(previous["genetic_identity_screen"]["candidate_ids"])):
        raise ValueError("Incomplete or duplicated current/historical public comparison universe")
    preflight_path = project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
    preflight = json.loads(preflight_path.read_text())
    fixed_ids = preflight["references"][0]["sample_ids"]
    union = set(fixed_ids)
    groups = []
    inputs = project / "results_v3_model_family_screen_05/inputs"
    for role, filename in (("ancestry_reference", "reference.retained.vcf.gz"),
                           ("pure_training", "pure_train.canonical.vcf.gz"),
                           ("pure_development", "pure_pure.canonical.vcf.gz")):
        path = inputs / filename
        identities = sample_ids(path)
        additional = [identity for identity in identities if identity not in union]
        union.update(identities)
        groups.append(dict(role=role, path=str(path), ids=additional, all_ids=identities))
    if union & set(public_ids):
        raise ValueError("Literal public/reference identity overlap must be resolved before panel construction")
    jobs = {}
    for chromosome in range(1, 23):
        scope = dense_run / ("smoke" if chromosome == 22 else "full")
        new = scope / f"normalized/chr{chromosome}/normalized.vcf.gz"
        old = previous_root / f"normalized/chr{chromosome}/normalized.vcf.gz"
        verification = scope / f"verification/chr{chromosome}.verified.json"
        new_files = [new, Path(str(new) + ".csi"), verification]
        old_files = [old, Path(str(old) + ".csi")]
        bindings = {str(path): dense["published_sha256"][str(path)] for path in new_files}
        bindings.update({str(path): previous["published_sha256"][str(path)] for path in old_files})
        jobs[str(chromosome)] = dict(new=str(new), old=str(old), verification=str(verification),
                                     input_sha256=bindings)
    labels = {person["sample_id"]: person["source"] for person in original["people"] + selected["people"]}
    historical_roles = {}
    for group in groups:
        for identity in group["all_ids"]:
            historical_roles.setdefault(identity, []).append(group["role"])
    bindings = {str(path): digest(path) for path in
                (dense_manifest, selected_path, previous_path, original_selection, preflight_path,
                 project / "config/sgdp_normalization.json", project / "config/dense_public_support.json")}
    bindings.update({group["path"]: digest(group["path"]) for group in groups})
    policy = copy.deepcopy(json.loads((project / "config/sgdp_normalization.json").read_text()))
    policy["relatedness"]["comparison"] = (
        f"All{len(candidates)} current public candidates plus all{len(public_ids) - len(candidates)} "
        "historical excluded public controls versus the "
        "fixed imputation reference, all historical fitting/development identities and each other.")
    Path("policy/config").mkdir(parents=True)
    write_json("policy/config/sgdp_normalization.json", policy)
    write_json("plan.json", dict(
        candidate_ids=candidates, public_ids=public_ids, new_ids=new_ids, historical_ids=historical,
        historical_control_ids=sorted(set(historical) - set(candidates)), sources=labels,
        recalled_control_id="S_Chane-1.DG", minimum_control_shared_calls=1000,
        minimum_control_concordance=0.995, coverage_eligible_ids=dense["eligible_by_real_mask_coverage"],
        fixed_ids=fixed_ids, groups=groups, historical_roles=historical_roles,
        all_ids=fixed_ids + [identity for group in groups for identity in group["ids"]] + public_ids,
        jobs=jobs, marker_vcf=groups[0]["path"], input_sha256=bindings,
        production_preprocessing_changed=False, ancestry_models_authorized=False))
    print(f"Genetic screen: {len(candidates)} candidates + {len(public_ids) - len(candidates)} "
          f"historical controls; {len(union) + len(public_ids)} comparison people.", flush=True)


def checked_sources(plan, chromosome, new, old, verification):
    job = plan["jobs"][str(chromosome)]
    require_hashes({str(new): job["input_sha256"][job["new"]],
                    str(old): job["input_sha256"][job["old"]],
                    str(new) + ".csi": job["input_sha256"][job["new"] + ".csi"],
                    str(old) + ".csi": job["input_sha256"][job["old"] + ".csi"],
                    str(verification): job["input_sha256"][job["verification"]]})
    record = json.loads(verification.read_text())
    if (record["status"] != "coordinate_allele_and_genotype_preservation_verified"
            or record["normalized_sha256"] != digest(new)
            or record["chromosome"] != chromosome
            or sample_ids(new) != record["sample_ids"]
            or sample_ids(old) != plan["historical_ids"]
            or set(record["sample_ids"]) != set(plan["new_ids"]) | {plan["recalled_control_id"]}):
        raise ValueError("Independent normalization/sample axes differ from frozen input plan")


def shared_keys(plan, chromosome, new, old, positions):
    markers = query_keys(Path(plan["marker_vcf"]), chromosome)
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({k[1] for k in markers})))
    keys = set(markers) & set(query_keys(new, chromosome, positions)) & set(query_keys(old, chromosome, positions))
    if not keys or len(keys) != len({key[1] for key in keys}):
        raise ValueError("Empty or allele-ambiguous common identity markers")
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({k[1] for k in keys})))
    return keys


def block(path, positions, identities, keys):
    result = {}
    for key, values in genotype_rows(path, positions, identities):
        if key in keys:
            if key in result:
                raise ValueError("Duplicated exact genotype key")
            values = [value.replace("|", "/") for value in values]
            if any(value not in ("0/0", "0/1", "1/0", "1/1", "./.") for value in values):
                raise ValueError("Non-biallelic diploid observed genotype")
            result[key] = values
    if set(result) != keys:
        raise ValueError("Incomplete exact-key genotype block")
    return result


def merge_candidates(plan, chromosome, new, old, verification, output):
    output.mkdir()
    checked_sources(plan, chromosome, new, old, verification)
    positions = output / "positions.tsv"
    keys = shared_keys(plan, chromosome, new, old, positions)
    new_values = block(new, positions, plan["new_ids"], keys)
    old_values = block(old, positions, plan["historical_ids"], keys)
    path = output / "candidates.vcf"
    with path.open("w") as handle:
        handle.write(f'##fileformat=VCFv4.2\n##contig=<ID={chromosome}>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Observed unphased GT">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                     "\t".join(plan["public_ids"]) + "\n")
        for chrom, pos, ref, alt in sorted(keys):
            key = (chrom, pos, ref, alt)
            handle.write(f"{chrom}\t{pos}\t{chrom}:{pos}:{ref}:{alt}\t{ref}\t{alt}\t.\t.\t.\tGT\t" +
                         "\t".join(new_values[key] + old_values[key]) + "\n")
    subprocess.run(["bgzip", str(path)], check=True)
    subprocess.run(["bcftools", "index", str(path) + ".gz"], check=True)
    write_json(output / "merge.json", dict(chromosome=chromosome, sites=len(keys),
        sample_ids=plan["public_ids"], purpose="identity_screen_only_not_ancestry_training_truth",
        ancestry_models_authorized=False))


def verify_merge(plan, chromosome, new, old, verification, folder, output):
    checked_sources(plan, chromosome, new, old, verification)
    positions = output.with_suffix(".positions.tsv")
    keys = shared_keys(plan, chromosome, new, old, positions)
    calls = folder / "candidates.vcf.gz"
    if sample_ids(calls) != plan["public_ids"]:
        raise ValueError("Merged candidate sample axis changed")
    if (set(query_keys(calls, chromosome)) != keys
            or int(subprocess.check_output(["bcftools", "index", "-n", str(calls)], text=True)) != len(keys)):
        raise ValueError("Merged candidate loci differ from the complete exact-key intersection")
    actual = block(calls, positions, plan["public_ids"], keys)
    start = 0
    for source, identities in ((new, plan["new_ids"]), (old, plan["historical_ids"])):
        for key, values in block(source, positions, identities, keys).items():
            if actual[key][start:start + len(identities)] != values:
                raise ValueError("Merged observed genotypes differ from source")
        start += len(identities)
    control = plan["recalled_control_id"]
    first = block(new, positions, [control], keys)
    second = block(old, positions, [control], keys)
    pairs = [(first[key][0], second[key][0]) for key in keys
             if first[key][0] != "./." and second[key][0] != "./."]
    matches = sum(sorted(a.split("/")) == sorted(b.split("/")) for a, b in pairs)
    if (len(pairs) < plan["minimum_control_shared_calls"]
            or matches / len(pairs) < plan["minimum_control_concordance"]):
        raise ValueError("Recalled/reused Chane identity-control concordance failed")
    write_json(output, dict(
        status="coordinate_allele_and_genotype_preservation_verified",
        verification_kind="exact_identity_subset_merge_of_independently_normalized_observed_sources",
        chromosome=chromosome, sites=len(keys), sample_ids=plan["public_ids"],
        genotype_cells_verified=len(keys) * len(plan["public_ids"]),
        normalized_sha256=digest(calls), normalized_index_sha256=digest(str(calls) + ".csi"),
        control_shared_calls=len(pairs), control_concordance=matches / len(pairs),
        scientific_ancestry_launch_authorized=False))
    print(f"chr{chromosome}: {len(keys):,} observed identity markers; Chane control "
          f"{matches}/{len(pairs)}; no missing calls filled.", flush=True)


def verify_panel(plan, chromosome, folder, merged, production_root, preflight, output):
    report = json.loads((folder / "panel.json").read_text())
    calls = folder / "panel.vcf.gz"
    if (report["chromosome"] != chromosome or report["sample_ids"] != plan["all_ids"]
            or report["candidate_ids"] != plan["public_ids"] or sample_ids(calls) != plan["all_ids"]
            or digest(calls) != report["panel_sha256"]):
        raise ValueError("Identity comparison universe/published panel changed")
    reference = next(row for row in preflight["references"] if row["chromosome"] == chromosome)
    fixed = production_root / "beagle_references" / Path(reference["reference"]).name
    require_hashes({path: item["sha256"] for path, item in reference["files"].items()})
    positions = output.with_suffix(".positions.tsv")
    keys = set(query_keys(Path(plan["marker_vcf"]), chromosome))
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({k[1] for k in keys})))
    sources = [(fixed, plan["fixed_ids"])] + [
        (Path(group["path"]), group["ids"]) for group in plan["groups"] if group["ids"]]
    sources += [(merged / "candidates.vcf.gz", plan["public_ids"])]
    for path, _ in sources:
        keys &= set(query_keys(path, chromosome, positions))
    if (len(keys) != report["sites"] or set(query_keys(calls, chromosome)) != keys
            or int(subprocess.check_output(["bcftools", "index", "-n", str(calls)], text=True)) != len(keys)):
        raise ValueError("Identity panel dropped or invented exact comparison markers")
    positions.write_text("".join(f"{chromosome}\t{pos}\t{pos}\n" for pos in sorted({k[1] for k in keys})))
    called = Counter()
    checked = 0
    for source, identities in sources:
        expected = ((key, values) for key, values in genotype_rows(source, positions, identities) if key in keys)
        actual = genotype_rows(calls, positions, identities)
        for a, b in itertools.zip_longest(actual, expected):
            if a is None or b is None or a[0] != b[0]:
                raise ValueError("Independent identity-panel row/sample coverage changed")
            values = [value.replace("|", "/") for value in b[1]]
            if a[1] != values:
                raise ValueError("Identity-panel genotype changed from original comparison source")
            for identity, value in zip(identities, values):
                called[identity] += value != "./."
            checked += len(identities)
    if dict(called) != report["called"] or checked != len(keys) * len(plan["all_ids"]):
        raise ValueError("Independent panel completeness/callability accounting differs")
    write_json(output, dict(
        status="every_identity_panel_genotype_compared_to_original_source",
        chromosome=chromosome, sites=len(keys), sample_ids=plan["all_ids"],
        public_ids=plan["public_ids"], genotype_cells_verified=checked,
        panel_sha256=digest(calls), panel_index_sha256=digest(str(calls) + ".csi"),
        panel_report_sha256=digest(folder / "panel.json"), ancestry_models_authorized=False))
    print(f"chr{chromosome}: independently checked all {checked:,} comparison genotype cells.", flush=True)


def classify(plan, related, policy):
    candidates, fixed = set(plan["candidate_ids"]), set(plan["fixed_ids"])
    direct_fixed, conflicts_excluded, bridges, conflicts = set(), set(), [], []
    parents = {}

    def find(identity):
        parents.setdefault(identity, identity)
        if parents[identity] != identity:
            parents[identity] = find(parents[identity])
        return parents[identity]

    for identity in candidates:
        find(identity)
    for pair in related:
        first, second = pair["first"], pair["second"]
        parents[find(first)] = find(second)
        if {first, second} & fixed:
            direct_fixed.update({first, second} & candidates)
        if pair["genotype_concordance"] >= policy["duplicate_concordance_threshold"]:
            for candidate, partner in ((first, second), (second, first)):
                if candidate in candidates and partner not in candidates and partner not in fixed:
                    bridges.append(dict(candidate=candidate, historical_identity=partner,
                                        historical_roles=plan["historical_roles"].get(partner, []),
                                        evidence=pair, must_be_absent_from_new_fitting_roles=True))
            if first in candidates and second in candidates and plan["sources"][first] != plan["sources"][second]:
                conflicts.append(pair)
                conflicts_excluded.update((first, second))
    fixed_roots = {find(identity) for identity in sorted(parents) if identity in fixed}
    fixed_excluded = {identity for identity in candidates if find(identity) in fixed_roots}
    eligible = (candidates & set(plan["coverage_eligible_ids"])) - fixed_excluded - conflicts_excluded
    components = {}
    for identity in sorted(parents):
        components.setdefault(find(identity), []).append(identity)
    groups = sorted(components.values())
    capacity = {source: len({find(identity) for identity in eligible if plan["sources"][identity] == source})
                for source in SOURCES}
    return dict(
        fixed_reference_direct_related_exclusions=sorted(direct_fixed),
        fixed_reference_connected_exclusions=sorted(fixed_excluded),
        identity_label_conflict_exclusions=sorted(conflicts_excluded),
        eligible_candidate_ids=sorted(eligible), eligible_source_counts=dict(Counter(plan["sources"][i] for i in eligible)),
        independent_component_capacity_by_source=capacity, family_identity_components=groups,
        historical_identity_bridges=bridges, conflicting_candidate_identity_labels=conflicts,
        at_least_two_components_per_source=all(number >= 2 for number in capacity.values()),
        historical_fitting_roles_not_authorized_for_reuse=True, donor_roles_frozen=False,
        ancestry_models_authorized=False)


def collect(project, plan_path, kinship, panels, proofs, output):
    plan = json.loads(plan_path.read_text())
    require_hashes(plan["input_sha256"])
    folders = {json.loads((folder / "panel.json").read_text())["chromosome"]: folder for folder in panels}
    records = [json.loads(path.read_text()) for path in proofs]
    if (len(folders) != 22 or len(records) != 22 or sorted(folders) != list(range(1, 23))
            or sorted(row["chromosome"] for row in records) != list(range(1, 23))):
        raise ValueError("Exactly one verified comparison panel/report per autosome required")
    for chromosome in range(1, 23):
        folder = folders[chromosome]
        record = next(row for row in records if row["chromosome"] == chromosome)
        require_hashes({str(folder / "panel.vcf.gz"): record["panel_sha256"],
                        str(folder / "panel.vcf.gz.csi"): record["panel_index_sha256"],
                        str(folder / "panel.json"): record["panel_report_sha256"]})
        if (record["chromosome"] != chromosome or record["sample_ids"] != plan["all_ids"]
                or record["public_ids"] != plan["public_ids"]):
            raise ValueError("Incomplete all-autosome panel verification")
    policy = json.loads((project / "config/sgdp_normalization.json").read_text())["relatedness"]
    with (kinship / "relatedness.kin0").open() as handle:
        names = handle.readline().lstrip("#").split()
        pairs = [dict(zip(names, line.split())) for line in handle if line.strip()]
    independent = audit_pairs(pairs, plan["all_ids"], plan["public_ids"], policy)
    original = json.loads((kinship / "genetic_eligibility.json").read_text())
    if any(original[key] != value for key, value in independent.items()):
        raise ValueError("Independent KING pair-coverage audit differs")
    result = classify(plan, independent["candidate_involving_related_pairs"], policy)
    result.update(status="all_autosome_genetic_screen_complete_new_cohort_freeze_required",
                  chromosomes=list(range(1, 23)), comparison_people=len(plan["all_ids"]),
                  current_candidates=len(plan["candidate_ids"]), historical_controls=plan["historical_control_ids"],
                  pairs_checked=independent["pairs_checked"], pruned_markers=original["pruned_markers"],
                  minimum_pair_called_markers=independent["minimum_pair_called_markers"],
                  genotype_cells_verified=sum(row["genotype_cells_verified"] for row in records),
                  all_related_pairs=independent["candidate_involving_related_pairs"],
                  input_sha256=plan["input_sha256"], source_purity_certified=False)
    write_json(output, result)
    print("All22 genetic screen complete. Independent component capacity: "
          f"{result['independent_component_capacity_by_source']}; no ancestry models authorized.", flush=True)


def finalize(run, post_cleanup):
    require_hashes(source_hashes(run.parent / "source_sha256.txt"))
    manifest = run / "published_verification.json"
    if post_cleanup:
        record = json.loads(manifest.read_text())
        require_hashes(record["published_sha256"])
        write_json(run / "post_cleanup_verification.json", dict(
            status="genetic_qc_hashes_reverified_after_exact_session_cleanup",
            manifest_sha256=digest(manifest), ancestry_models_authorized=False))
        return
    plan = json.loads((run / "preparation/plan.json").read_text())
    require_hashes(plan["input_sha256"])
    chromosomes = [22] if run.name == "smoke" else list(range(1, 22))
    hashes = {}
    for chromosome in chromosomes:
        folder = run / f"identity_panels/chr{chromosome}"
        record = json.loads((run / f"identity_panels/chr{chromosome}.verified.json").read_text())
        if record["chromosome"] != chromosome or record["sample_ids"] != plan["all_ids"]:
            raise ValueError("Unexpected verified panel scope")
        require_hashes({str(folder / "panel.vcf.gz"): record["panel_sha256"],
                        str(folder / "panel.vcf.gz.csi"): record["panel_index_sha256"],
                        str(folder / "panel.json"): record["panel_report_sha256"]})
    if run.name == "full":
        smoke = run.parent / "smoke"
        if digest(run / "preparation/plan.json") != digest(smoke / "preparation/plan.json"):
            raise ValueError("Genetic screen plan changed after real task smoke")
        record = json.loads((smoke / "published_verification.json").read_text())
        require_hashes(record["published_sha256"])
        hashes.update(record["published_sha256"])
        hashes[str(smoke / "published_verification.json")] = digest(smoke / "published_verification.json")
        summary = json.loads((run / "summary/genetic_screen.json").read_text())
        if summary["chromosomes"] != list(range(1, 23)) or summary["ancestry_models_authorized"] is not False:
            raise ValueError("All-autosome genetic screen summary required")
    for name in ("preparation", "candidate_merges", "identity_panels", "kinship", "summary"):
        for path in (run / name).rglob("*"):
            if path.is_file():
                hashes[str(path.resolve())] = digest(path)
    write_json(manifest, dict(status="genetic_screen_evidence_verified_not_donor_roles_frozen",
                              published_sha256=hashes, ancestry_models_authorized=False))
    require_hashes(hashes)
    print(f"Durable genetic evidence verified for {chromosomes}.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "merge", "verify-merge", "verify-panel", "collect", "finalize"))
    for name in ("project", "dense-run", "plan", "new", "old", "verification", "folder",
                 "merged", "production-root", "preflight", "run", "kinship", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--panels", type=Path, nargs="+")
    parser.add_argument("--reports", type=Path, nargs="+")
    parser.add_argument("--post-cleanup", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.project, args.dense_run)
    elif args.stage == "finalize":
        finalize(args.run, args.post_cleanup)
    else:
        plan = json.loads(args.plan.read_text())
        if args.stage == "merge":
            merge_candidates(plan, args.chromosome, args.new, args.old, args.verification, args.output)
        elif args.stage == "verify-merge":
            verify_merge(plan, args.chromosome, args.new, args.old, args.verification, args.folder, args.output)
        elif args.stage == "verify-panel":
            verify_panel(plan, args.chromosome, args.folder, args.merged, args.production_root,
                         json.loads(args.preflight.read_text()), args.output)
        else:
            collect(args.project, args.plan, args.kinship, args.panels, args.reports, args.output)
