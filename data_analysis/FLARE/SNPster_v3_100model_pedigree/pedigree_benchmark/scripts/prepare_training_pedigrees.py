#!/usr/bin/env python3
"""Build explicitly training-only pedigrees without weakening held-out safeguards."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from make_plans import generate
from prepare_founders import (
    SOURCES, load_label_metadata, parse_tsv, read_reference_panel, sample_ids,
    validate_labels,
)
from validate_frozen_candidates import sha256_file


def validate_pool(training, validation, labels, vcf_ids, reference_ids, evaluation_donors):
    validate_labels(training, labels, "Training founders")
    validate_labels(validation, labels, "Pure validation")
    if not training or set(training) != set(vcf_ids):
        raise ValueError("Training VCF and metadata IDs must match exactly")
    forbidden = set(validation) | set(reference_ids) | set(evaluation_donors)
    overlaps = set(training) & forbidden
    if overlaps:
        raise ValueError(f"Training founders overlap reference/development donors: {sorted(overlaps)}")
    if {label[0] for label in training.values()} != set(SOURCES):
        raise ValueError("Training founders must cover all eight sources")


def run(training_meta, validation_meta, label_meta, founder_vcf, reference_vcf,
        reference_panel, evaluation_plan, outdir, families=108, batch_size=10, seed=20261014):
    training, validation = parse_tsv(training_meta), parse_tsv(validation_meta)
    labels = load_label_metadata(label_meta)
    reference_labels = read_reference_panel(reference_panel)
    reference_ids = sample_ids(reference_vcf)
    if set(reference_labels) != reference_ids:
        raise ValueError("Actual reference VCF and panel IDs differ")
    for donor, fine in reference_labels.items():
        if donor not in labels or labels[donor][1] != fine:
            raise ValueError(f"Reference-panel label mismatch: {donor}")
    with Path(evaluation_plan).open() as stream:
        evaluation = list(csv.DictReader(stream, delimiter="\t"))
    evaluation_donors = {r["founder_donor"] for r in evaluation
                         if r["is_founder"].lower() in ("true", "1")}
    if not evaluation_donors:
        raise ValueError("Missing development founder assignments")
    validate_pool(training, validation, labels, sample_ids(founder_vcf),
                  reference_ids, evaluation_donors)
    outdir = Path(outdir)
    if outdir.exists():
        raise ValueError("Training plans already exist; do not overwrite a cohort")
    generate(training_meta, [], outdir, families, batch_size, seed, True,
             design="full", family_prefix="trainfam")
    with (outdir / "expected_ancestry.tsv").open() as stream:
        nodes = list(csv.DictReader(stream, delimiter="\t"))
    used = {r["founder_donor"] for r in nodes if r["is_founder"].lower() in ("true", "1")}
    if set(r["sample_id"] for r in nodes) & set(r["sample_id"] for r in evaluation):
        raise ValueError("Training and development node IDs collide")
    audit = {
        "cohort_role": "decoder_training", "pedigree_design": "full",
        "reference_disjoint": True, "labels_validated": True,
        "development_donor_overlap": 0, "pure_validation_donor_overlap": 0,
        "founders_are_existing_pure_training_donors": True,
        "uses_old_mosaic_training_examples": False,
        "assay": "phased_panel_ready_founder_VCF",
        "n_pool_donors": len(training), "n_used_donors": len(used),
        "pool_per_source": dict(Counter(group for group, fine in training.values())),
        "used_per_source": dict(Counter(training[d][0] for d in used)),
        "used_donor_ids": sorted(used), "development_donor_ids": sorted(evaluation_donors),
        "n_families": families, "n_nodes": len(nodes),
        "n_descendants": sum(r["is_founder"].lower() in ("false", "0") for r in nodes),
        "chromosomes": list(range(1, 23)),
        "input_sha256": {str(Path(p).resolve()): sha256_file(p) for p in (
            training_meta, validation_meta, label_meta, founder_vcf, reference_vcf,
            reference_panel, evaluation_plan, Path(__file__), Path(__file__).with_name("make_plans.py"),
        )},
    }
    (outdir / "training_cohort_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, indent=2))
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("training-meta", "validation-meta", "label-meta", "founder-vcf",
                 "reference-vcf", "reference-panel", "evaluation-plan", "outdir"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--families", type=int, default=108)
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261014)
    args = parser.parse_args()
    run(**vars(args))


if __name__ == "__main__":
    main()
