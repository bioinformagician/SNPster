#!/usr/bin/env python3
"""Freeze prediction-independent donors and audit legacy as well as public identities."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import subprocess

from audit_production_ancestry import SOURCES, digest, rows, write_json
from sgdp_genetic_qc import audit_pairs, run
from verify_sgdp_full_manifest import require_hashes


def psam_ids(path):
    with path.open() as handle:
        fields = handle.readline().lstrip("#").split()
        if "IID" not in fields:
            raise ValueError("PLINK sample metadata has no IID column")
        index = fields.index("IID")
        identities = [line.split()[index] for line in handle if line.strip()]
    if not identities or len(identities) != len(set(identities)):
        raise ValueError("Empty or duplicated PLINK sample identities")
    return identities


def freeze(project, output):
    output.mkdir()
    protocol_path = project / "config/production_matched_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if tuple(protocol["sources"]) != SOURCES or len(protocol["candidates"]) != 10:
        raise ValueError("Require the declared eight-source, ten-configuration plan")
    old_path = project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
    normalized_path = project / "results_v3_sgdp_normalization_full_v2_20261022_095828/summary/normalization_complete.json"
    old, normalized = [json.loads(path.read_text()) for path in (old_path, normalized_path)]
    require_hashes(normalized["published_sha256"])
    for relative, expected in json.loads((project / "config/ancestry_benchmark_contract_v2.json").read_text())[
            "production_preprocessing"]["source_sha256"].items():
        if digest(project.parents[2] / "data pipeline" / relative) != expected:
            raise ValueError(f"Production source changed: {relative}")
    authority_path = Path("/srv/dependencies/ancestry/flare2_test/source_model/ancestry_role_audit.tsv")
    if digest(authority_path) != old["input_sha256"][str(authority_path)]:
        raise ValueError("Authoritative ancestry labels changed")
    authority = {}
    for record in rows(authority_path):
        if record["source_group"]:
            label = (record["source_group"], record["fine"])
            if record["sample"] in authority and authority[record["sample"]] != label:
                raise ValueError("Conflicting authoritative source labels")
            authority[record["sample"]] = label
    selection_path = project / "results_v3_public_sgdp_smoke_20261021_203540/metadata/selected.json"
    public = json.loads(selection_path.read_text())["people"]
    public_ids = set(normalized["genetic_identity_screen"]["eligible_candidate_ids"])
    candidates = [dict(row, public=False) for row in old["cohorts"]["pure_development"]["direct_id_disjoint_people"]]
    candidates += [dict(sample_id=row["sample_id"], source=row["source"], fine="SGDP_" + row["population"],
                        public=True) for row in public if row["sample_id"] in public_ids]
    candidate_ids = [row["sample_id"] for row in candidates]
    base = normalized_path.parents[1] / "kinship"
    sample_axis = psam_ids(base / "all.psam")
    if len(candidate_ids) != len(set(candidate_ids)) or not set(candidate_ids) <= set(sample_axis):
        raise ValueError("Unknown or duplicate evaluation candidate identity")
    candidate_file = output / "candidates.id"
    candidate_file.write_text("".join(f"0\t{identity}\n" for identity in candidate_ids))
    run(["plink2", "--pfile", base / "all", "--extract", base / "prune.prune.in",
         "--make-king-table", "cols=id,nsnp,hethet,ibs0,ibs1,kinship",
         "--king-table-require", candidate_file, "--threads", "1", "--memory", "4096",
         "--out", output / "evaluation"], output / "genetic_screen")
    with (output / "evaluation.kin0").open() as handle:
        fields = handle.readline().lstrip("#").split()
        pairs = [dict(zip(fields, line.split())) for line in handle if line.strip()]
    policy = json.loads((project / "config/sgdp_normalization.json").read_text())["relatedness"]
    audit = audit_pairs(pairs, sample_axis, candidate_ids, policy)
    fixed = set(old["references"][0]["sample_ids"])
    fixed_related, fitting_excluded = set(), set()
    for pair in audit["candidate_involving_related_pairs"]:
        first, second = pair["first"], pair["second"]
        for identity, partner in ((first, second), (second, first)):
            if identity in candidate_ids:
                if partner in fixed:
                    fixed_related.add(identity)
                elif partner not in candidate_ids:
                    fitting_excluded.add(partner)
    usable = [row for row in candidates if row["sample_id"] not in fixed_related]
    by_source = defaultdict(list)
    for row in usable:
        if row["source"] not in SOURCES:
            raise ValueError("Unknown operational source label")
        if not row["public"] and authority.get(row["sample_id"]) != (row["source"], row["fine"]):
            raise ValueError("Legacy evaluation source/fine label differs from authority")
        by_source[row["source"]].append(row)
    if set(by_source) != set(SOURCES):
        write_json(output / "blocked.json", dict(
            status="blocked_missing_genetically_independent_source", source_counts=dict(Counter(
                row["source"] for row in usable)), audit=audit, fixed_reference_related=sorted(fixed_related)))
        raise ValueError("All eight sources require genetically screened reference-disjoint founders")
    founders = [min(by_source[source], key=lambda row: (not row["public"], row["sample_id"]))
                for source in SOURCES]
    inputs = project / "results_v3_model_family_screen_05/inputs"
    roles = {}
    for name, filename in (("reference", "reference_labels.tsv"), ("training", "pure_train_labels.tsv")):
        records = rows(inputs / filename)
        selected = []
        for row in records:
            identity = row["sample_id"]
            label = (row["original_source_group"], row["fine"])
            if authority.get(identity) != label:
                raise ValueError(f"Authoritative fitting label mismatch: {identity}")
            if identity in fixed and identity not in fitting_excluded:
                selected.append(dict(sample_id=identity, source=label[0], fine=label[1], public=False))
        if name == "training":
            selected = [row for source in SOURCES for row in
                        sorted((r for r in selected if r["source"] == source), key=lambda r: r["sample_id"])[:2]]
        if set(row["source"] for row in selected) != set(SOURCES):
            raise ValueError(f"Fitting role lacks eight-source coverage: {name}")
        if name == "training" and Counter(row["source"] for row in selected) != Counter({source: 2 for source in SOURCES}):
            raise ValueError("Require exactly two training founders per source")
        roles[name] = selected
    groups = [set(row["sample_id"] for row in records) for records in (founders, *roles.values())]
    if any(first & second for i, first in enumerate(groups) for second in groups[i + 1:]):
        raise ValueError("Evaluation/reference/training donor overlap")
    evaluation = [dict(sample_id=f"pure_{source}", kind="pure", donors=[row["sample_id"]],
                       family_id="evaluation_founders", source=source)
                  for source, row in zip(SOURCES, founders)]
    for index, source in enumerate(SOURCES):
        other = (index + 1) % len(SOURCES)
        evaluation.append(dict(sample_id=f"f1_{source}_{SOURCES[other]}", kind="F1",
                               donors=[founders[index]["sample_id"], founders[other]["sample_id"]],
                               family_id="evaluation_founders"))
    training = [dict(sample_id=f"train_{index}", kind="pure", donors=[row["sample_id"]],
                     source=row["source"], family_id=f"training:{row['sample_id']}")
                for index, row in enumerate(roles["training"])]
    people = [dict(row, role=role) for role, records in (("evaluation", evaluation), ("training", training))
              for row in records]
    for mask_index, mask in enumerate(protocol["masks"]):
        for index, row in enumerate(people):
            row.setdefault("production_ids", {})[mask] = str(100000 + mask_index * 1000 + index)
    original_cohort = json.loads((project / "results_v3_model_family_screen_05/cohort.json").read_text())
    evidence = [protocol_path, old_path, normalized_path, authority_path, selection_path,
                inputs / "reference_labels.tsv", inputs / "pure_train_labels.tsv",
                inputs / "pure_pure.canonical.vcf.gz", inputs / "pure_train.canonical.vcf.gz",
                project / "results_v3_model_family_screen_05/cohort.json"]
    masks = {mask: str(project / f"results_v3_phasefree_array_validation_01/{mask}.observed_keys.tsv")
             for mask in protocol["masks"]}
    evidence += [Path(path) for path in masks.values()]
    record = dict(
        status="cohort_and_candidate_plan_frozen_before_predictions",
        protocol=protocol, evaluation_founders=founders, reference=roles["reference"],
        training_founders=roles["training"], people=people, masks=masks,
        chromosome_spans=original_cohort["chromosome_spans"],
        normalized_vcfs=normalized["normalized_vcfs"],
        fixed_reference_related_candidates=sorted(fixed_related),
        fitting_related_exclusions=sorted(fitting_excluded),
        genetic_audit=audit, independent_final_test=False,
        production_preprocessing_changed=False,
        input_sha256={str(path): digest(path) for path in evidence},
    )
    write_json(output / "cohort.json", record)
    print(f"Frozen eight-source cohort: {len(evaluation)} evaluation people, {len(training)} training people, "
          f"{len(roles['reference'])} ancestry-reference people; {audit['pairs_checked']} genetic pairs.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    freeze(args.project, args.output)
