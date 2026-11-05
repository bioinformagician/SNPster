#!/usr/bin/env python3
"""Audit unused donors and production-QC coverage without changing a frozen model."""

import argparse
import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path
import pickle
import runpy
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd

import model_family_gnomix as adapter
from screen_dosage_models import genotype_dosage


PROJECT = Path(__file__).resolve().parents[2]


def unused_donors(pools, roles, dense):
    used, group_counts, fine_sources = set(), {}, {}
    for label, records, key in pools:
        identities = {row[key] for row in records}
        if len(identities) != len(records):
            raise ValueError(f"Duplicate identities in {label}")
        if key == "sample_id" and identities & used:
            raise ValueError(f"Fitted/development identity overlap in {label}")
        group_counts[label] = len(identities)
        used.update(identities)
        if key != "sample_id":
            continue
        for row in records:
            source = row["source_group"]
            source = "EAS_SEA" if source in ("EAS", "SEA") else source
            if source not in adapter.SOURCES:
                raise ValueError(f"Unknown current source: {source}")
            if fine_sources.setdefault(row["fine"], source) != source:
                raise ValueError(f"Inconsistent current population labels: {row['fine']}")
    unused = {}
    for row in roles:
        if row["sample"] in used:
            continue
        value = dict(
            sample_id=row["sample"], fine=row["fine"], historical_role=row["role"],
            current_source=fine_sources.get(row["fine"], ""),
            dense_available=row["sample"] in dense,
            eligible_for_known_source_truth=(
                row["role"] == "source" and row["fine"] in fine_sources),
        )
        if row["sample"] in unused and unused[row["sample"]] != value:
            raise ValueError(f"Inconsistent unused-donor metadata: {row['sample']}")
        unused[row["sample"]] = value
    records = sorted(unused.values(), key=lambda row: row["sample_id"])
    eligible = [row for row in records if row["eligible_for_known_source_truth"]]
    report = dict(
        fitted_and_inspected_identity_union=len(used), group_counts=group_counts,
        dense_reference_n=len(dense), unused_annotated_identities=len(records),
        unused_annotated_dense_identities=sum(row["dense_available"] for row in records),
        unused_current_known_source_identities=len(eligible),
        unused_current_known_source_dense_identities=sum(row["dense_available"] for row in eligible),
        unused_role_counts={role: sum(row["historical_role"] == role for row in records)
                            for role in sorted({row["historical_role"] for row in records})},
        unannotated_dense_identities=len(dense - {row["sample"] for row in roles} - used),
        limitation="Recorded sample identities only; not kinship or biological-alias certification.",
    )
    return report, records


def load_production_qc(directory):
    config = runpy.run_path(str(directory / "config.py"))
    thresholds = SimpleNamespace(
        gp_min=config["GP_MIN"], ds_tol=config["DS_TOL"],
        snps_only=config["SNPS_ONLY"], biallelic_only=config["BIALLELIC_ONLY"])
    tree = ast.parse((directory / "qc_classes.py").read_text())
    container = next(node for node in tree.body
                     if isinstance(node, ast.ClassDef) and node.name == "ImputedDataContainer")
    method = next(node for node in container.body
                  if isinstance(node, ast.FunctionDef) and node.name == "qc_imputed_data")
    # Execute the unchanged QC method without importing the unrelated Polars VCF loader.
    namespace = {"np": np}
    exec(compile(ast.Module(body=[method], type_ignores=[]),
                 str(directory / "qc_classes.py"), "exec"), namespace)
    return namespace["qc_imputed_data"], thresholds


def qc_kept_indices(frame, method, thresholds):
    case = SimpleNamespace(imputed_data=frame, qc_thresholds=thresholds, qc_status=False)
    with contextlib.redirect_stdout(io.StringIO()):
        method(case)
    return set(case.imputed_data.index) if case.qc_status else set()


def qc_frame(path, sample, keys):
    requested = set(keys)
    if len(requested) != len(keys) or not keys:
        raise ValueError("Empty/duplicate frozen marker axis")
    command = ["bcftools", "query", "-f",
               "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT\t%DS\t%GP]\n", str(path)]
    text = subprocess.check_output(command, text=True)
    values = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 8:
            raise ValueError(f"Expected one sample and GT/DS/GP fields: {path}")
        key = (fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3])
        if key not in requested:
            continue
        if key in values or fields[4] != sample:
            raise ValueError(f"Duplicate marker or incorrect sample in {path}")
        genotype_dosage(fields[5])
        ds, gp = float(fields[6]), np.array(fields[7].split(","), dtype=float)
        if (not np.isfinite(ds) or not 0 <= ds <= 2 or gp.shape != (3,)
                or not np.all(np.isfinite(gp)) or np.any(gp < 0) or np.any(gp > 1)
                or abs(gp.sum() - 1) > .020001):
            raise ValueError(f"Invalid biallelic DS/GP at {key}")
        values[key] = [key[2], key[3], fields[5], ds, *gp]
    if len(values) != len(keys):
        raise ValueError(f"Missing {len(keys) - len(values)} frozen markers in {path}")
    return pd.DataFrame([values[key] for key in keys],
                        columns=["REF", "ALT", "GT", "DS", "GP_00", "GP_01", "GP_11"])


def audit(args):
    if (args.output / "donor_audit.json").exists() or (args.output / "readiness.json").exists():
        raise ValueError("Readiness audit requires a fresh output")
    args.output.mkdir(parents=True, exist_ok=True)
    stack = PROJECT / "results_v3_phasefree_stack_01"
    manifest_path = stack / "selected_candidate.json"
    manifest = json.loads(manifest_path.read_text())
    model_path = stack / manifest["model_path"]
    model_hash = adapter.digest(model_path)
    if model_hash != manifest["model_sha256"]:
        raise ValueError("Frozen candidate digest changed")
    inputs = {}

    def record(path):
        value = adapter.digest(path)
        if str(path) in inputs and inputs[str(path)] != value:
            raise ValueError(f"Input changed during audit: {path}")
        inputs[str(path)] = value
        return path

    pools = []
    for label, path, key in (
        ("reference_fit", args.screen / "inputs/reference_labels.tsv", "sample_id"),
        ("pure_training_and_pedigree_founder_pool",
         args.screen / "inputs/pure_train_labels.tsv", "sample_id"),
        ("pure_development", args.screen / "inputs/pure.samples.tsv", "sample_id"),
        ("eligible_development_founder_pool",
         PROJECT / "results_v3_pedigree_pilot_10/pedigree/inputs/eligible_heldout_founders.tsv",
         "sample"),
        ("historically_inspected_additional_founders",
         PROJECT / "config/pedigree_additional_founders.tsv", "sample"),
    ):
        pools.append((label, adapter.records(record(path)), key))
    roles = adapter.records(record(args.source_metadata / "ancestry_role_audit.tsv"))
    dense_text = subprocess.check_output(["bcftools", "query", "-l", str(args.dense_vcf)], text=True)
    dense = dense_text.splitlines()
    if len(set(dense)) != len(dense):
        raise ValueError("Duplicate dense-reference identities")
    donor_report, unused = unused_donors(pools, roles, set(dense))
    (args.output / "dense_sample_ids.txt").write_text(dense_text)
    pd.DataFrame(unused, columns=[
        "sample_id", "fine", "historical_role", "current_source",
        "dense_available", "eligible_for_known_source_truth",
    ]).to_csv(args.output / "unused_annotated_identities.tsv", sep="\t", index=False)
    adapter.write_json(args.output / "donor_audit.json", donor_report)

    record(manifest_path)
    record(Path(__file__))
    record(Path(adapter.__file__))
    record(Path(__file__).with_name("screen_dosage_models.py"))
    record(args.qc_module / "config.py")
    record(args.qc_module / "qc_classes.py")
    method, thresholds = load_production_qc(args.qc_module)
    with model_path.open("rb") as stream:
        bundle = pickle.load(stream)
    keys = bundle["marker_keys"]
    if len(keys) != manifest["n_markers"] or len(set(keys)) != len(keys):
        raise ValueError("Frozen marker count mismatch")
    pilot = json.loads(record(args.matched / "pilot.json").read_text())
    samples = [row["sample_id"] for row in pilot["samples"] if row["role"] == "evaluation"]
    if len(samples) != 5 or len(set(samples)) != 5:
        raise ValueError("Expected the five previously evaluated array controls")
    sample_rows, chromosome_rows = [], []
    for sample in samples:
        total_kept = 0
        for chrom in range(1, 23):
            relative = Path(f"chr{chrom}/{sample}/truth_panel.vcf.gz")
            path = args.matched / relative
            previous = json.loads(record(args.matched / f"chr{chrom}/audit.json").read_text())
            if previous["status"] != "complete":
                raise ValueError(f"Incomplete preserved chromosome: {chrom}")
            record(path)
            if inputs[str(path)] != previous["output_sha256"][str(Path("/pilot") / relative)]:
                raise ValueError(f"Preserved imputation output changed: {path}")
            chrom_keys = [key for key in keys if key[0] == str(chrom)]
            frame = qc_frame(path, sample, chrom_keys)
            kept = len(qc_kept_indices(frame, method, thresholds))
            total_kept += kept
            chromosome_rows.append(dict(sample_id=sample, chrom=chrom,
                                        required=len(frame), kept=kept, removed=len(frame) - kept))
        sample_rows.append(dict(sample_id=sample, required=len(keys), kept=total_kept,
                                removed=len(keys) - total_kept,
                                retained_fraction=total_kept / len(keys),
                                compatible_with_frozen_complete_gt_input=total_kept == len(keys)))
        print(json.dumps(sample_rows[-1]), flush=True)
    pd.DataFrame(chromosome_rows).to_csv(
        args.output / "qc_coverage_by_chromosome.tsv", sep="\t", index=False)
    pd.DataFrame(sample_rows).to_csv(args.output / "qc_coverage.tsv", sep="\t", index=False)
    if adapter.digest(model_path) != model_hash:
        raise ValueError("Model changed during the audit")
    for path, expected in inputs.items():
        if adapter.digest(path) != expected:
            raise ValueError(f"Input changed during audit: {path}")
    report = dict(
        status="complete_readiness_audit_not_independent_validation",
        candidate=manifest["candidate"], model_sha256=model_hash,
        donor_audit=donor_report, production_qc_thresholds=vars(thresholds),
        production_qc_coverage=sample_rows,
        complete_gt_input_compatible_after_qc=all(
            row["compatible_with_frozen_complete_gt_input"] for row in sample_rows),
        production_qc_scope="Current module defaults in an isolated runtime; no deployed environment audited.",
        new_validation_cohort_started=False, model_fitting=False,
        new_accuracy_results=False, production_promotion=False, pgs_calibration_changed=False,
        dense_vcf=str(args.dense_vcf),
        dense_sample_axis_sha256=hashlib.sha256(dense_text.encode()).hexdigest(),
        input_sha256={**inputs, str(model_path): model_hash},
        output_sha256={str(args.output / name): adapter.digest(args.output / name) for name in (
            "dense_sample_ids.txt", "unused_annotated_identities.tsv", "donor_audit.json",
            "qc_coverage.tsv", "qc_coverage_by_chromosome.tsv")},
    )
    adapter.write_json(args.output / "readiness.json", report)
    print("Readiness audit complete; no models or acceptance gates changed.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--matched", type=Path, required=True)
    parser.add_argument("--source-metadata", type=Path, required=True)
    parser.add_argument("--dense-vcf", type=Path, required=True)
    parser.add_argument("--qc-module", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    audit(parser.parse_args())
