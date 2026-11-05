#!/usr/bin/env python3
"""Freeze the requested dense development split without consulting predictions."""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import subprocess

from audit_production_ancestry import SOURCES, digest, rows, write_json
from freeze_production_cohort import psam_ids
from pilot_matched_imputation import ancestor_donors
from sgdp_genetic_qc import audit_pairs, run, sample_ids
from verify_sgdp_full_manifest import require_hashes


def related_to(identities, pairs):
    identities = set(identities)
    result = set(identities)
    for pair in pairs:
        first, second = pair["first"], pair["second"]
        if first in identities:
            result.add(second)
        if second in identities:
            result.add(first)
    return result


def complete_families(families, retained):
    retained = set(retained)
    return {family for family, donors in families.items() if set(donors) <= retained}


def genotype_digest(path, identities, log):
    value = hashlib.sha256()
    with log.open("w") as error:
        process = subprocess.Popen(
            ["bcftools", "query", "-S", str(identities), "-f",
             "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", str(path)],
            stdout=subprocess.PIPE, stderr=error)
        try:
            with process.stdout:
                for block in iter(lambda: process.stdout.read(1024 * 1024), b""):
                    value.update(block)
            if process.wait():
                raise ValueError(f"Reference genotype verification failed: {log}")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    return value.hexdigest()


def purge_reference(cohort_path, chromosome, output):
    cohort = json.loads(cohort_path.read_text())
    source = Path("/srv/dependencies/imputation_runner/imputer/beagle_references") / (
        f"hgdp1kgp_chr{chromosome}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz")
    if digest(source) != cohort["template_sha256"][str(source)]:
        raise ValueError("Frozen synthetic genome template changed")
    if set(sample_ids(source)) != set(cohort["original_imputation_reference_ids"]):
        raise ValueError("Dense reference sample axis changes across chromosomes")
    output.mkdir(exist_ok=True)
    selected = output / f"chr{chromosome}.keep.ids"
    selected.write_text("".join(identity + "\n" for identity in cohort["imputation_reference_ids"]))
    target = output / source.name
    run(["bcftools", "view", "-S", selected, "-Oz", "-o", target, source], output / "purge")
    subprocess.run(["bcftools", "index", str(target)], check=True)
    if sample_ids(target) != cohort["imputation_reference_ids"]:
        raise ValueError("Benchmark imputation reference has an incorrect sample axis")
    before = genotype_digest(source, selected, output / "original_query.log")
    after = genotype_digest(target, selected, output / "copy_query.log")
    if before != after:
        raise ValueError("Reference membership change altered allele keys or retained phased genotypes")
    write_json(output / f"chr{chromosome}.reference.json", dict(
        status="all_retained_reference_alleles_and_phased_genotypes_preserved",
        chromosome=chromosome, source_sha256=digest(source), copy_sha256=digest(target),
        index_sha256=digest(Path(str(target) + ".csi")), genotype_stream_sha256=after,
        retained_ids=cohort["imputation_reference_ids"], excluded_ids=cohort["imputation_excluded"],
        cohort_sha256=digest(cohort_path), production_resources_modified=False))
    print(f"chr{chromosome}: donor-purged benchmark reference and every retained genotype verified.", flush=True)


def local_truth(plan, bp, sample, sources, span):
    mapping = {}
    for founder, row in plan["founders"].items():
        donor = row["vcf_sample"]
        for hap_index, hap in enumerate(("h0", "h1")):
            identities = {part[2] for part in bp[founder][hap]}
            if len(identities) != 1 or next(iter(identities)) in mapping:
                raise ValueError("Ambiguous or duplicated founder homolog identity")
            mapping[next(iter(identities))] = (donor, sources.get(donor), hap_index)
    segments = []
    amounts = {source: 0 for source in SOURCES}
    for hap_index, hap in enumerate(("h0", "h1")):
        previous = span[0]
        for start, inclusive_end, identity in bp[sample][hap]:
            first, end = max(start, span[0]), min(inclusive_end + 1, span[1])
            if end <= first:
                continue
            if first != previous or identity not in mapping:
                raise ValueError("Pedigree ancestry has a gap, overlap or unknown homolog")
            donor, source, donor_haplotype = mapping[identity]
            if source not in SOURCES:
                raise ValueError(f"Inherited donor lacks an authoritative native source: {donor}")
            segments.append(dict(haplotype=hap_index, start=first, end=end,
                                 donor=donor, source=source, donor_haplotype=donor_haplotype))
            amounts[source] += end - first
            previous = end
        if previous != span[1]:
            raise ValueError("Pedigree ancestry does not cover the complete physical span")
    denominator = 2 * (span[1] - span[0])
    return segments, [amounts[source] / denominator for source in SOURCES]


def freeze(project, output):
    from realized_ancestry import load_bp

    output.mkdir()
    policy_path = project / "config/dense_split_protocol.json"
    policy = json.loads(policy_path.read_text())
    base_path = project / "config" / policy["base_protocol"]
    protocol = copy.deepcopy(json.loads(base_path.read_text()))
    selected_path = project / "results_v3_observed_dense_01/selection.json"
    selected = json.loads(selected_path.read_text())
    legacy_path = project / "results_v3_phasefree_array_validation_01/experiment.json"
    legacy = json.loads(legacy_path.read_text())
    original_path = project / "results_v3_model_family_screen_05/cohort.json"
    original = json.loads(original_path.read_text())
    authority_path = Path("/srv/dependencies/ancestry/flare2_test/source_model/ancestry_role_audit.tsv")
    authority = {}
    for row in rows(authority_path):
        if row["source_group"]:
            label = (row["source_group"], row["fine"])
            if row["sample"] in authority and authority[row["sample"]] != label:
                raise ValueError("Conflicting authoritative native source labels")
            authority[row["sample"]] = label
    roles = {}
    evidence = [policy_path, base_path, selected_path, legacy_path, original_path, authority_path]
    for role, filename, key in (("reference", "reference_labels.tsv", "reference_ids"),
                                ("training", "pure_train_labels.tsv", "pure_ids")):
        path = project / "results_v3_observed_dense_01/inputs" / filename
        evidence.append(path)
        records = rows(path)
        identities = [row["sample_id"] for row in records]
        if len(identities) != len(set(identities)) or set(identities) != set(selected[key]):
            raise ValueError("Stored requested split and source-label identities disagree")
        roles[role] = []
        for row in records:
            if authority.get(row["sample_id"]) != (row["original_source_group"], row["fine"]):
                raise ValueError(f"Authoritative source mismatch: {row['sample_id']}")
            roles[role].append(dict(sample_id=row["sample_id"], source=row["original_source_group"],
                                    fine=row["fine"], public=False))
    requested = (len(roles["reference"]), len(roles["training"]), selected["n_safe_training_families"],
                 selected["n_training_descendants"], len(legacy["samples"]))
    if requested != (3033, 240, 25, 76, 14):
        raise ValueError(f"Requested legacy split changed: {requested}")
    evaluation_ids = set(selected["evaluation_donors"])
    if (set(selected["reference_ids"]) & set(selected["pure_ids"])
            or evaluation_ids & (set(selected["reference_ids"]) | set(selected["pure_ids"]))):
        raise ValueError("Requested reference/training/evaluation roles overlap")
    all_founders = sorted(evaluation_ids | set(selected["pure_ids"]))
    genetic_root = project / "results_v3_dense_donor_qc_20261022_192734/full"
    genetic = genetic_root / "kinship"
    genetic_manifest = genetic_root / "published_verification.json"
    published = json.loads(genetic_manifest.read_text())["published_sha256"]
    genetic_files = [genetic / name for name in ("all.psam", "all.pgen", "all.pvar", "prune.prune.in")]
    require_hashes({str(path): published[str(path)] for path in genetic_files})
    evidence.append(genetic_manifest)
    axis = psam_ids(genetic / "all.psam")
    if not set(all_founders) <= set(axis):
        raise ValueError("Requested founders are absent from the verified all-autosome genetic panel")
    (output / "founders.id").write_text("".join(f"0\t{identity}\n" for identity in all_founders))
    run(["plink2", "--pfile", genetic / "all", "--extract", genetic / "prune.prune.in",
         "--make-king-table", "cols=id,nsnp,hethet,ibs0,ibs1,kinship",
         "--king-table-require", output / "founders.id", "--threads", "1", "--memory", "4096",
         "--out", output / "founders"], output / "genetic_screen")
    with (output / "founders.kin0").open() as handle:
        fields = handle.readline().lstrip("#").split()
        pairs = [dict(zip(fields, line.split())) for line in handle if line.strip()]
    relatedness_path = project / "config/sgdp_normalization.json"
    audit = audit_pairs(pairs, axis, all_founders, json.loads(relatedness_path.read_text())["relatedness"])
    forbidden = related_to(evaluation_ids, audit["candidate_involving_related_pairs"])
    for role in roles:
        roles[role] = [row for row in roles[role] if row["sample_id"] not in forbidden]
        if {row["source"] for row in roles[role]} != set(SOURCES):
            write_json(output / "blocked.json", dict(reason="native_source_lost_after_family_purge",
                       role=role, source_counts=dict(Counter(row["source"] for row in roles[role]))))
            raise ValueError(f"Relatedness purge removes native-source coverage from {role}")
    families_path = project / "results_v3_phasefree_rich_01/donors.json"
    families = json.loads(families_path.read_text())
    safe = complete_families(families, [row["sample_id"] for row in roles["training"]])
    train_root = project / "results_v3_pedigree_training_screen_01/cohort/pedigree"
    eval_root = project / "results_v3_pedigree_pilot_10/pedigree"
    people, plans = [], {}
    for row in roles["training"]:
        people.append(dict(sample_id=row["sample_id"], kind="pure", role="training",
                           source=row["source"], donors=[row["sample_id"]],
                           family_id=f"training:{row['sample_id']}"))
    for batch in selected["pedigree_train"]:
        name = batch["batch"]
        samples_path = project / "results_v3_observed_dense_01/inputs" / Path(batch["samples_tsv"]).name
        plan_path = train_root / f"plans/plans/{name}/pedigrees.json"
        evidence += [samples_path, plan_path]
        plan = json.loads(plan_path.read_text())
        for row in rows(samples_path):
            if row["family_id"] in safe:
                people.append(dict(sample_id=row["sample_id"], kind="pedigree", role="training",
                                   donors=sorted(ancestor_donors(plan, row["sample_id"])),
                                   family_id=row["family_id"], plan_key=f"training:{name}"))
        plans[f"training:{name}"] = (plan_path, train_root / f"chromosomes/pdb_{name}_chr{{chromosome}}.bp")
    for row in legacy["samples"]:
        if row["track"] == "pure":
            people.append(dict(sample_id=row["sample_id"], kind="pure", role="evaluation",
                               source=authority[row["sample_id"]][0], donors=[row["sample_id"]],
                               family_id=row["family_id"]))
        else:
            name = row["batch"]
            plan_path = eval_root / f"plans/plans/{name}/pedigrees.json"
            plans[f"evaluation:{name}"] = (plan_path, eval_root / f"chromosomes/pdb_{name}_chr{{chromosome}}.bp")
            plan = json.loads(plan_path.read_text())
            people.append(dict(sample_id=row["sample_id"], kind="pedigree", role="evaluation",
                               donors=sorted(ancestor_donors(plan, row["sample_id"])),
                               family_id=row["family_id"], plan_key=f"evaluation:{name}"))
    sea = policy["additional_pure_evaluation_donor"]
    if sea not in evaluation_ids or authority[sea][0] != "SEA":
        raise ValueError("Additional SEA control must already be a held-out SEA founder")
    people.append(dict(sample_id=sea, kind="pure", role="evaluation", source="SEA", donors=[sea],
                       family_id=f"pure:{sea}"))
    if len({row["sample_id"] for row in people}) != len(people):
        raise ValueError("Duplicate simulated training/evaluation person identity")
    for mask_index, mask in enumerate(protocol["masks"]):
        for index, row in enumerate(people):
            row.setdefault("production_ids", {})[mask] = str(100000 + mask_index * 10000 + index)
    source_by_id = {identity: authority[identity][0] for identity in all_founders}
    local, truth = {}, {}
    for chromosome in range(1, 23):
        span = original["chromosome_spans"][str(chromosome)]
        decoded = {}
        for name, (plan_path, bp_template) in plans.items():
            bp_path = Path(str(bp_template).format(chromosome=chromosome))
            evidence += [plan_path, bp_path]
            decoded[name] = (json.loads(plan_path.read_text()), load_bp(bp_path, chromosome))
        local[str(chromosome)], truth[str(chromosome)] = {}, {}
        for row in people:
            if row["kind"] == "pure":
                segments = [dict(haplotype=hap, start=span[0], end=span[1],
                                 donor=row["donors"][0], source=row["source"],
                                 donor_haplotype=hap) for hap in (0, 1)]
                value = [float(source == row["source"]) for source in SOURCES]
            else:
                plan, bp = decoded[row["plan_key"]]
                segments, value = local_truth(plan, bp, row["sample_id"], source_by_id, span)
            local[str(chromosome)][row["sample_id"]] = segments
            truth[str(chromosome)][row["sample_id"]] = value
    retained_founders = evaluation_ids | {row["sample_id"] for row in roles["training"]}
    imputation_excluded = related_to(retained_founders, audit["candidate_involving_related_pairs"])
    dense_root = Path("/srv/dependencies/imputation_runner/imputer")
    fixed_ids = sample_ids(dense_root / "beagle_references" /
                          "hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz")
    if not retained_founders <= set(fixed_ids):
        raise ValueError("A requested synthetic founder is absent from the dense template")
    external_path = genetic_root / "summary/genetic_screen.json"
    external = json.loads(external_path.read_text())["eligible_candidate_ids"]
    if len(external) != 67 or set(external) & retained_founders:
        raise ValueError("The67-person external observed-data holdout changed or overlaps this split")
    protocol["evaluation"] = dict(interpretation=policy["truth_interpretation"], people=15,
                                  independent_final_test=False, selection=policy["evaluation_policy"])
    protocol["training"] = dict(reference=policy["split_policy"], pure_people=len(roles["training"]),
                                pedigree_people=sum(row["role"] == "training" and row["kind"] == "pedigree"
                                                    for row in people), no_calibration=True)
    protocol["production"]["founder_phase"] = "Frozen synthetic template homologs define the genome; do not rephase synthetic truth. Raw microarray uploads use unchanged production phasing and imputation."
    evidence += [families_path, relatedness_path, external_path, genetic / "all.psam",
                 genetic / "all.pgen", genetic / "all.pvar", genetic / "prune.prune.in"]
    masks = {mask: str(project / f"results_v3_phasefree_array_validation_01/{mask}.observed_keys.tsv")
             for mask in protocol["masks"]}
    evidence += [Path(path) for path in masks.values()]
    record = dict(status="cohort_and_candidate_plan_frozen_before_predictions", protocol=protocol,
                  truth_design=policy["truth_design"], split_policy=policy,
                  reference=roles["reference"], training_founders=roles["training"],
                  evaluation_founders=[dict(sample_id=identity, source=authority[identity][0],
                                           fine=authority[identity][1]) for identity in sorted(evaluation_ids)],
                  people=people, masks=masks, chromosome_spans=original["chromosome_spans"],
                  local_ancestry=local, truth_by_chromosome=truth,
                  template_sha256={str(dense_root / "beagle_references" / Path(path).name): value
                                   for path, value in legacy["reference_sha256"].items()},
                  imputation_reference_ids=sorted(set(fixed_ids) - imputation_excluded),
                  original_imputation_reference_ids=fixed_ids,
                  imputation_excluded=sorted(set(fixed_ids) & imputation_excluded),
                  fitting_related_exclusions=sorted(forbidden - evaluation_ids),
                  requested_counts=dict(reference=3033, pure_training=240, training_families=25,
                                        training_descendants=76, evaluation=15),
                  retained_counts=dict(reference=len(roles["reference"]), pure_training=len(roles["training"]),
                                       training_families=len(safe),
                                       training_descendants=protocol["training"]["pedigree_people"], evaluation=15),
                  reserved_external_ids=external, genetic_audit=audit,
                  independent_final_test=False, production_preprocessing_changed=False,
                  input_sha256={str(path): digest(path) for path in set(evidence)})
    write_json(output / "cohort.json", record)
    print(f"Frozen dense synthetic development split: {record['retained_counts']}; "
          f"benchmark imputation reference {len(record['imputation_reference_ids'])}/{len(fixed_ids)}; "
          f"{audit['pairs_checked']:,} genetic pairs verified;67 external people reserved.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--reference-chromosome", type=int)
    args = parser.parse_args()
    if args.reference_chromosome is None:
        freeze(args.project, args.output)
    else:
        purge_reference(args.cohort, args.reference_chromosome, args.output)
