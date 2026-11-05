#!/usr/bin/env python3
"""Frozen training-only FLARE copying/mismatch diagnostic on matched array GTs."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import zipfile

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_dense_observed as dense
import benchmark_observed_ancestry as common
import model_family_gnomix as adapter
import pilot_matched_imputation as shared
from screen_model_families import physical_predictions, read_model, write_model
from verify_dense_observed import check_case, local_physical


PROJECT = Path(__file__).resolve().parents[2]
PROTOCOL = PROJECT / "config/flare_native_diagnostic.json"
BASE = Path("/benchmark")
STOCK = Path("/opt/flare/flare.jar")


def select_controls(rows, forbidden, seed):
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise ValueError("Duplicate training identity")
    selected = []
    for source in adapter.SOURCES:
        choices = [row for row in rows if row["source_group"] == source]
        if not choices:
            raise ValueError(f"Missing training source: {source}")
        row = min(choices, key=lambda r: hashlib.sha256(
            f"{seed}:flare-native:{r['sample_id']}".encode()).digest())
        if row["sample_id"] in forbidden:
            raise ValueError("Training-only diagnostic overlaps reference/evaluation identities")
        selected.append(row)
    return selected


def anchored_model(template, fraction, theta=None):
    names, panels = template["ancestries"], template["panels"]
    if (not 0 <= fraction <= 1 or len(set(names)) != len(names) or len(set(panels)) != len(panels)
            or set(names) != set(adapter.SOURCES) or set(panels) != set(names)):
        raise ValueError("Anchoring requires identical named seven-source ancestry/panel axes")
    identity = np.asarray([[float(name == panel) for panel in panels] for name in names])
    result = dict(template, p=(1 - fraction) * template["p"] + fraction * identity)
    if theta is not None:
        if not 0 < theta < .5:
            raise ValueError("Invalid controlled mismatch probability")
        result["theta"] = np.full_like(template["theta"], theta)
    return result


def source_audit(reference_path, labels, training, controls):
    ids, _ = adapter.vcf_header(reference_path)
    if set(ids) != {row["sample_id"] for row in labels} or len(ids) != len(labels):
        raise ValueError("Reference VCF/metadata identity coverage differs")
    authority = pd.read_csv("/data/flare2_test/source_model/ancestry_role_audit.tsv", sep="\t")
    authority = authority[authority.dataset.eq("reference")].set_index("sample")
    panel = dict(line.split() for line in (BASE / "inputs/reference.panel").read_text().splitlines())
    if set(panel) != set(ids):
        raise ValueError("Reference panel identity mapping changed")
    for row in labels:
        record = authority.loc[row["sample_id"]]
        group = "EAS_SEA" if record.source_group in ("EAS", "SEA") else record.source_group
        if group != row["source_group"] or panel[row["sample_id"]] != group or record.fine != row["fine"]:
            raise ValueError("Reference authority/source/fine label differs")
    forbidden = set(ids) | set(common.read("/legacy/experiment.json")["donor_ids"])
    if forbidden & {row["sample_id"] for row in training}:
        raise ValueError("Reference/evaluation donors overlap diagnostic training pool")
    if not {r["sample_id"] for r in controls} <= {r["sample_id"] for r in training}:
        raise ValueError("Selected positive control is not training-only")
    return dict(n_reference=len(ids), authoritative_labels_verified=True, donor_overlap=0)


def jar_check(precision):
    with zipfile.ZipFile(STOCK) as a, zipfile.ZipFile(precision) as b:
        left = {name: hashlib.sha256(a.read(name)).hexdigest() for name in a.namelist() if name.endswith(".class")}
        right = {name: hashlib.sha256(b.read(name)).hexdigest() for name in b.namelist() if name.endswith(".class")}
    changed = sorted(name for name in set(left) | set(right) if left.get(name) != right.get(name))
    if changed != ["admix/EstimatedGlobalAncProportions.class"]:
        raise ValueError(f"Unexpected scientific jar class differences: {changed}")
    return dict(stock_sha256=adapter.digest(STOCK), precision_sha256=adapter.digest(precision),
                changed_classes=changed)


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Fresh diagnostic output required")
    parent = common.read(BASE / "experiment.json")
    shared.verify_hashes(parent["implementation_sha256"])
    training = adapter.records(BASE / "inputs/pure_train_labels.tsv")
    reference = adapter.records(BASE / "inputs/reference_labels.tsv")
    forbidden = {r["sample_id"] for r in reference} | set(common.read("/legacy/experiment.json")["donor_ids"])
    protocol = common.read(PROTOCOL)
    controls = select_controls(training, forbidden, protocol["selection_seed"])
    inputs = [PROTOCOL, BASE / "experiment.json", BASE / "selection.json", BASE / "inputs/pure_train_labels.tsv",
              BASE / "inputs/reference_labels.tsv", BASE / "inputs/reference.panel",
              Path("/data/flare2_test/source_model/ancestry_role_audit.tsv"), STOCK]
    audits = []
    for mask in protocol["masks"]:
        task = BASE / f"tasks/{mask}_chr22"
        checkpoint = common.read(task / "inputs_complete.json")
        shared.verify_hashes(checkpoint["output_sha256"])
        cohort = common.read(task / "cohort.json")
        audits.append(source_audit(cohort["reference"]["vcf"], reference, training, controls))
        inputs += [task / "inputs_complete.json", task / "cohort.json", task / "markers.tsv",
                   Path(cohort["reference"]["vcf"]), Path(cohort["pure_train"]["vcf"]),
                   Path(cohort["flare_jar"]), Path(cohort["genetic_map"])]
    cohort = common.read(BASE / "tasks/lonni_chr22/cohort.json")
    jar = jar_check(cohort["flare_jar"])
    sources = [Path(__file__), Path(dense.__file__), Path(common.__file__), Path(adapter.__file__),
               Path(shared.__file__), Path(__file__).with_name("screen_model_families.py"),
               Path(__file__).with_name("verify_dense_observed.py"),
               PROJECT / "run_flare_native_diagnostic.sh"]
    adapter.write_json(out / "experiment.json", dict(
        protocol=protocol, controls=controls, source_audits=audits, jar_audit=jar,
        input_sha256=shared.fingerprint_files(inputs), implementation_sha256=shared.fingerprint_files(sources),
        evaluation_people_included=False, independent_final_test=False, production_promotion=False))
    pd.DataFrame(controls).to_csv(out / "controls.tsv", sep="\t", index=False)
    adapter.write_json(out / "state.json", dict(status="prepared_training_only"))


def command(cohort, target, prefix, em=False, model=None, stock=False):
    result = common.flare_command(cohort, target, prefix, dict(array=True, em=em))
    if model is not None:
        result.append(f"model={model}")
    if stock:
        result[result.index("-jar") + 1] = str(STOCK)
    return result


def run_mask(out, mask):
    config = common.read(out / "experiment.json")
    controls = config["controls"]
    task = BASE / f"tasks/{mask}_chr22"
    cohort = common.read(task / "cohort.json")
    markers = pd.read_csv(task / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(r.chrom, int(r.pos), r.ref, r.alt) for r in markers.itertuples(index=False)]
    span = cohort["chromosome_spans"]["22"]
    directory = out / mask
    directory.mkdir()
    cases = []
    for row in controls:
        name = row["sample_id"]
        work = directory / name
        work.mkdir()
        inherited = adapter.load_haplotypes(cohort["pure_train"]["vcf"], [name], 22, keys)
        case = dense.phase_case(work, name, 22, keys, inherited, cohort)
        n, switches, pairs = check_case((case, keys))
        if n != len(keys):
            raise ValueError("Independent phase-only GT verification failed")
        case.update(source_group=row["source_group"], phase_switches=switches, phase_pairs=pairs)
        cases.append(case)
    template_prefix = directory / "template"
    shared.checked(command(cohort, cases[0]["inherited_phase"], template_prefix), directory, "template")
    template = read_model(str(template_prefix) + ".model")
    template_file = directory / "native_fixed.model"
    write_model(template_file, template)
    models = {"native_fixed": str(template_file), "native_em": str(template_file)}
    for name, policy in config["protocol"]["variants"].items():
        if policy["model"] == "anchored":
            path = directory / f"{name}.model"
            write_model(path, anchored_model(template, policy["anchor_fraction"], policy.get("theta")))
            models[name] = str(path)
    rows, stock_audit = [], []
    for variant, policy in config["protocol"]["variants"].items():
        work = directory / variant
        work.mkdir()
        for case in cases:
            for stage in config["protocol"]["stages"]:
                name = case["sample_id"]
                prefix = work / f"{name}_{stage}"
                seconds = shared.checked(command(
                    cohort, case[stage], prefix, policy["em"], models[variant]), work, prefix.name)
                model = read_model(str(prefix) + ".model")
                ids, p = physical_predictions(str(prefix) + ".anc.vcf.gz", model["ancestries"], {"22": span}, markers)
                if ids != [name]:
                    raise ValueError("FLARE diagnostic sample-axis mismatch")
                independent = local_physical(str(prefix) + ".anc.vcf.gz", model["ancestries"], span, keys)
                np.testing.assert_allclose(p[0], independent, rtol=0, atol=1e-10)
                target = np.array([float(s == case["source_group"]) for s in adapter.SOURCES])
                rows.append(dict(mask=mask, variant=variant, stage=stage, sample_id=name,
                                 source_group=case["source_group"], chromosome=22,
                                 tv=float(np.abs(p[0] - target).sum() / 2), seconds=seconds,
                                 **{f"pred_{s}": float(p[0, i]) for i, s in enumerate(adapter.SOURCES)}))
        print(f"{mask}: {variant} complete on seven training controls / both phases", flush=True)
    sas = next(case for case in cases if case["source_group"] == "SAS")
    for stage in config["protocol"]["stages"]:
        prefix = directory / f"stock_SAS_{stage}"
        shared.checked(command(cohort, sas[stage], prefix, model=template_file, stock=True), directory, prefix.name)
        model = read_model(str(prefix) + ".model")
        _, p = physical_predictions(str(prefix) + ".anc.vcf.gz", model["ancestries"], {"22": span}, markers)
        precision = next(row for row in rows if row["variant"] == "native_fixed"
                         and row["sample_id"] == sas["sample_id"] and row["stage"] == stage)
        error = float(np.max(np.abs(p[0] - [precision[f"pred_{s}"] for s in adapter.SOURCES])))
        if error > 1e-12:
            raise ValueError("Stock and precision-export build differ in local ancestry")
        stock_audit.append(dict(stage=stage, max_fraction_error=error))
    frame = pd.DataFrame(rows)
    frame.to_csv(directory / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(cases).to_json(directory / "cases.json", orient="records", indent=2)
    adapter.write_json(directory / "verification.json", dict(
        status="complete_training_only_chr22_diagnostic", stock_precision_agree=True,
        stock_checks=stock_audit, independent_local_aggregations=len(rows),
        no_changed_gt=True, evaluation_people_included=False,
        output_sha256=shared.fingerprint_files([directory / "predictions.tsv", directory / "cases.json"])))
    return mask


def run(out):
    config = common.read(out / "experiment.json")
    shared.verify_hashes(config["input_sha256"])
    shared.verify_hashes(config["implementation_sha256"])
    adapter.write_json(out / "state.json", dict(status="running_training_only"))
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(run_mask, out, mask) for mask in config["protocol"]["masks"]]
        for future in as_completed(futures):
            print(f"Verified training-only {future.result()} diagnostics", flush=True)
    frame = pd.concat([pd.read_csv(out / mask / "predictions.tsv", sep="\t", float_precision="round_trip")
                       for mask in config["protocol"]["masks"]], ignore_index=True)
    frame.to_csv(out / "predictions.tsv", sep="\t", index=False)
    frame.groupby(["mask", "variant", "stage"]).agg(
        mean_tv=("tv", "mean"), max_tv=("tv", "max"), n_people=("sample_id", "size")
    ).reset_index().to_csv(out / "comparison.tsv", sep="\t", index=False)
    adapter.write_json(out / "completion.json", dict(
        status="complete_training_only_chr22_diagnostic", independent_final_test=False,
        production_promotion=False, n_unique_people=7,
        output_sha256=shared.fingerprint_files([out / "predictions.tsv", out / "comparison.tsv"])))
    adapter.write_json(out / "state.json", dict(status="complete_training_only_chr22_diagnostic"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["prepare", "run"], required=True)
    parser.add_argument("--outdir", type=Path, default=Path("/diagnostic"))
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.stage == "prepare":
            prepare(args.outdir)
        else:
            try:
                run(args.outdir)
            except Exception as error:
                adapter.write_json(args.outdir / "state.json", dict(status="failed", error=str(error)))
                raise
