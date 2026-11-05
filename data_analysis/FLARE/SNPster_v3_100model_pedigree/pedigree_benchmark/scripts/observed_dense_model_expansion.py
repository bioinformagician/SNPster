#!/usr/bin/env python3
"""Expand the verified observed-only benchmark without reintroducing imputation."""
import argparse
from collections import Counter, defaultdict
import itertools
from pathlib import Path
import pickle
import shutil
import subprocess
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_dense_observed as dense
import benchmark_observed_ancestry as common
import model_family_gnomix as native
import pilot_matched_imputation as shared
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import physical_predictions, read_model


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(native.SOURCES)


def read(path):
    return common.read(Path(path))


def rebase_retained(value, source):
    if isinstance(value, dict):
        return {rebase_retained(key, source): rebase_retained(item, source)
                for key, item in value.items()}
    if isinstance(value, list):
        return [rebase_retained(item, source) for item in value]
    if isinstance(value, str) and Path(value).is_relative_to("/benchmark"):
        return str(source / Path(value).relative_to("/benchmark"))
    if isinstance(value, str) and Path(value).is_relative_to("/workspace"):
        return str(PROJECT / Path(value).relative_to("/workspace"))
    if isinstance(value, str) and Path(value).is_relative_to("/results"):
        return str(PROJECT / "results_v3_model_family_screen_05" /
                   Path(value).relative_to("/results"))
    return value


def candidate_grid(plan):
    candidates = {}
    flare = plan["flare"]
    for reference, em, generations in itertools.product(
            flare["reference_panels"], flare["em"], flare["generations"]):
        name = f"flare_{reference}_array_{'em' if em else 'fixed'}_g{generations:03d}"
        candidates[name] = dict(
            family="flare", reference=reference, array=True, em=em, gen=generations)
    gnomix = plan["gnomix"]
    for window, c2 in itertools.product(gnomix["window_sizes_cM"], gnomix["crf_l2"]):
        name = f"gnomix_fast_{round(window * 100):03d}_l2_{c2:g}".replace(".", "p")
        candidates[name] = dict(
            family="gnomix", architecture="fast", reference="full",
            window_size_cM=window, smooth_size=1, context_ratio=0.5, crf_c2=c2)
    counts = Counter(policy["family"] for policy in candidates.values())
    if len(candidates) != 96 or counts != Counter(flare=64, gnomix=32):
        raise ValueError("Observed-only expansion must contain 64 FLARE and 32 Gnomix settings")
    return candidates


def freeze(source, plan_path, output):
    source_protocol = read(dense.PROTOCOL)
    completion = read(source / "completion.json")
    recovery = read(source / "recovery_state.json")
    if (completion["status"] != "complete_native_observed_only_pilot"
            or completion["imputation_performed"]
            or recovery["status"] != "complete_native_observed_only_pilot_after_explicit_recovery"
            or recovery["imputation_performed"]):
        raise ValueError("Expansion requires the complete recovered observed-only benchmark")
    retained_hashes = {
        str(source / Path(path).relative_to("/benchmark")) if Path(path).is_relative_to("/benchmark") else path: expected
        for path, expected in completion["output_sha256"].items()
    }
    shared.verify_hashes(retained_hashes)
    plan = read(plan_path)
    protocol = dict(source_protocol)
    protocol["candidates"] = candidate_grid(plan)
    protocol["model_expansion"] = dict(
        source=str(source), source_protocol_sha256=native.digest(dense.PROTOCOL),
        plan_sha256=native.digest(plan_path), candidate_selection_uses_performance=False,
        imputation_performed=False,
        omitted_controls=plan["phasefree_controls"],
        omitted_controls_reason="The approved observed-only protocol forbids imputed inputs and calibration controls")
    native.write_json(output, protocol)
    print("Frozen observed-only expansion: 64 FLARE and 32 Gnomix settings; no imputation controls.", flush=True)


def reference_records(cohort):
    rows = native.records(cohort["reference"]["labels_tsv"])
    if len(rows) != 3033 or {row["source_group"] for row in rows} != set(SOURCES):
        raise ValueError("Observed-only reference identities or source coverage changed")
    return rows


def prepare_references(source_task, cohort, output):
    rows = reference_records(cohort)
    result = {"full": (Path(cohort["reference"]["vcf"]), Path(cohort["reference"]["broad_panel"]))}
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["fine"]].append(row)
    for name, cap in (("finecap8", 8), ("finecap", 16), ("finecap32", 32)):
        selected = [row for fine in sorted(grouped)
                    for row in sorted(grouped[fine], key=lambda item: item["sample_id"])[:cap]]
        if {row["source_group"] for row in selected} != set(SOURCES):
            raise ValueError(f"{name} loses observed-only source coverage")
        ids = output / f"{name}.ids"
        panel = output / f"{name}.panel"
        vcf = output / f"{name}.vcf.gz"
        ids.write_text("".join(f"{row['sample_id']}\n" for row in selected))
        panel.write_text("".join(f"{row['sample_id']}\t{row['source_group']}\n" for row in selected))
        subprocess.run(["bcftools", "view", "-S", str(ids), "-Oz", "-o", str(vcf),
                        str(cohort["reference"]["vcf"])], check=True)
        subprocess.run(["bcftools", "index", str(vcf)], check=True)
        result[name] = (vcf, panel)
    return result


def flare_command(cohort, target, reference, panel, prefix, policy, nthreads=1):
    return [
        "java", "-Xmx6g", "-jar", cohort["flare_jar"], f"gt={target}", f"ref={reference}",
        f"ref-panel={panel}", f"map={cohort['genetic_map']}", f"out={prefix}",
        f"array={str(policy['array']).lower()}", f"em={str(policy['em']).lower()}",
        f"gen={policy['gen']}", "update-p=false",
        "probs=true", "min-maf=0", "min-mac=0", f"nthreads={nthreads}", f"seed={native.SEED}"
    ]


def prepare(source_task, protocol_path, output):
    output.mkdir()
    protocol = read(protocol_path)
    inputs = read(source_task / "inputs_complete.json")
    source = source_task.parents[1]
    shared.verify_hashes(rebase_retained(inputs["output_sha256"], source))
    cohort = rebase_retained(read(source_task / "cohort.json"), source)
    cases = rebase_retained(read(source_task / "cases.json"), source)
    frame = pd.read_csv(source_task / "markers.tsv", sep="\t", dtype={"chrom": str})
    chrom = int(frame.chrom.iloc[0])
    mask = read(source_task / "coverage.json")["mask"]
    if not cases or set(protocol["candidates"]) != set(candidate_grid(read(PROJECT / "config/dense_model_expansion.json"))):
        raise ValueError("Prepared axis does not match the frozen observed-only expansion")
    prepare_references(source_task, cohort, output)
    native.write_json(output / "cohort.json", cohort)
    native.write_json(output / "cases.json", cases)
    shutil.copy2(source_task / "markers.tsv", output / "markers.tsv")
    shutil.copy2(protocol_path, output / "protocol.json")
    native.write_json(output / "prepared.json", dict(
        status="complete_observed_only_model_inputs", chromosome=chrom, mask=mask,
        source_inputs_sha256=native.digest(source_task / "inputs_complete.json"),
        protocol_sha256=native.digest(protocol_path), candidates=len(protocol["candidates"]),
        imputation_performed=False))


def adopt_candidate(recovery, candidate, cases, protocol_path, output):
    if recovery is None:
        return None
    certificate = read(recovery / "complete.json")
    if (certificate["status"] != "complete_observed_only_model_expansion_task"
            or certificate["protocol_sha256"] != native.digest(protocol_path)
            or certificate["predictions_sha256"] != native.digest(recovery / "predictions.tsv")
            or certificate["imputation_performed"]):
        raise ValueError("Candidate recovery lacks a matching complete task certificate")
    source = recovery / candidate
    if not source.is_dir():
        raise ValueError(f"Candidate recovery is missing {candidate}")
    predictions = pd.read_csv(recovery / "predictions.tsv", sep="\t", float_precision="round_trip")
    predictions = predictions[predictions.candidate == candidate].copy()
    expected = {(case["sample_id"], stage) for case in cases
                for stage in ("inherited_phase", "array_phase_only")}
    actual = set(predictions[["sample_id", "stage"]].itertuples(index=False, name=None))
    if actual != expected or len(predictions) != len(expected):
        raise ValueError(f"Recovered candidate predictions are incomplete: {candidate}")
    shutil.copytree(source, output / "model")
    return predictions


def evaluate_candidate(prepared, protocol_path, candidate, output, recovery=None):
    output.mkdir()
    if recovery is not None and not recovery.is_dir():
        raise ValueError(f"Recovery path is not a directory: {recovery}")
    prepared_record = read(prepared / "prepared.json")
    protocol = read(protocol_path)
    if (prepared_record["status"] != "complete_observed_only_model_inputs"
            or prepared_record["protocol_sha256"] != native.digest(protocol_path)
            or prepared_record["imputation_performed"]):
        raise ValueError("Prepared model inputs do not match the frozen protocol")
    if candidate not in protocol["candidates"]:
        raise ValueError(f"Unknown frozen candidate: {candidate}")
    cohort, cases = read(prepared / "cohort.json"), read(prepared / "cases.json")
    frame = pd.read_csv(prepared / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(row.chrom, int(row.pos), row.ref, row.alt) for row in frame.itertuples(index=False)]
    chrom, mask = prepared_record["chromosome"], prepared_record["mask"]
    policy = protocol["candidates"][candidate]
    recovered = adopt_candidate(recovery, candidate, cases, protocol_path, output)
    rows = []
    if recovered is None:
        work = output / "model"
        work.mkdir()
        span = cohort["chromosome_spans"][str(chrom)]
        length = span[1] - span[0]
        started = time.monotonic()
        bundle = restored = None
        if policy["family"] == "gnomix":
            config = dict(policy, candidate=candidate, calibrate=False, retrain_base=False,
                          seed=native.SEED, sources=SOURCES)
            bundle = native.fit_chromosome(cohort, config, chrom, work, 1, 64, native.load_upstream())
            if bundle["keys"] != keys:
                raise ValueError("Gnomix changed the observed exact marker panel")
            with (work / "models" / f"chr{chrom}.pkl").open("rb") as stream:
                restored = pickle.load(stream)
        fit_seconds = time.monotonic() - started
        for case in cases:
            for stage in ("inherited_phase", "array_phase_only"):
                begin = time.monotonic()
                target = case[stage]
                prefix = work / f"{case['sample_id']}_{stage}"
                if policy["family"] == "flare":
                    name = policy["reference"]
                    reference = (Path(cohort["reference"]["vcf"]) if name == "full"
                                 else prepared / f"{name}.vcf.gz")
                    panel = (Path(cohort["reference"]["broad_panel"]) if name == "full"
                             else prepared / f"{name}.panel")
                    shared.checked(
                        flare_command(cohort, target, reference, panel, prefix, policy),
                        work, prefix.name)
                    model = read_model(str(prefix) + ".model")
                    ids, values = physical_predictions(
                        str(prefix) + ".anc.vcf.gz", model["ancestries"], {str(chrom): span}, frame)
                    if ids != [case["sample_id"]]:
                        raise ValueError("FLARE changed the evaluation sample axis")
                    prediction = values[0]
                else:
                    haplotypes = native.load_haplotypes(target, [case["sample_id"]], chrom, keys)
                    before = native.ordered_probabilities(
                        bundle["model"], bundle["model"].base.predict_proba(haplotypes))
                    after = native.ordered_probabilities(
                        restored["model"], restored["model"].base.predict_proba(haplotypes))
                    np.testing.assert_allclose(before, after, rtol=0, atol=1e-10)
                    _, physical = native.aggregate_probabilities(after, bundle["layout"])
                    prediction = physical[0] / length
                if (not np.isfinite(prediction).all() or np.any(prediction < 0)
                        or not np.isclose(prediction.sum(), 1)):
                    raise ValueError("Invalid observed-only ancestry fractions")
                rows.append(dict(
                    candidate=candidate, mask=mask, chromosome=chrom, span_bp=length,
                    sample_id=case["sample_id"], stage=stage, track=case["track"],
                    family_id=case["family_id"], fit_seconds=fit_seconds,
                    inference_seconds=time.monotonic() - begin,
                    **{f"pred_{source}": float(prediction[index])
                       for index, source in enumerate(SOURCES)}))
        predictions = pd.DataFrame(rows)
    else:
        predictions = recovered
    predictions_path = output / "predictions.tsv"
    predictions.to_csv(predictions_path, sep="\t", index=False)
    native.write_json(output / "complete.json", dict(
        status="complete_observed_only_model_candidate", chromosome=chrom, mask=mask,
        candidate=candidate, prepared_sha256=native.digest(prepared / "prepared.json"),
        protocol_sha256=native.digest(protocol_path),
        predictions_sha256=native.digest(predictions_path), recovered=recovered is not None,
        imputation_performed=False))


def evaluate(source_task, protocol_path, output, recovery=None):
    output.mkdir()
    if recovery is not None and not recovery.is_dir():
        raise ValueError(f"Recovery path is not a directory: {recovery}")
    protocol = read(protocol_path)
    inputs = read(source_task / "inputs_complete.json")
    source = source_task.parents[1]
    shared.verify_hashes(rebase_retained(inputs["output_sha256"], source))
    cohort = rebase_retained(read(source_task / "cohort.json"), source)
    cases = rebase_retained(read(source_task / "cases.json"), source)
    frame = pd.read_csv(source_task / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(row.chrom, int(row.pos), row.ref, row.alt) for row in frame.itertuples(index=False)]
    chrom = int(frame.chrom.iloc[0])
    mask = read(source_task / "coverage.json")["mask"]
    span = cohort["chromosome_spans"][str(chrom)]
    length = span[1] - span[0]
    references = prepare_references(source_task, cohort, output)
    rows = []
    for candidate, policy in protocol["candidates"].items():
        work = output / candidate
        recovered = recovery / candidate if recovery is not None and policy["family"] == "flare" else None
        if recovered and recovered.is_dir():
            expected_files = {
                f"{case['sample_id']}_{stage}{suffix}"
                for case in cases for stage in ("inherited_phase", "array_phase_only")
                for suffix in (".anc.vcf.gz", ".model", ".command.json")
            }
            actual_files = {path.name for path in recovered.iterdir() if path.is_file()}
            if not expected_files <= actual_files:
                raise ValueError(f"Incomplete FLARE recovery candidate: {candidate}")
            shutil.copytree(recovered, work)
        else:
            work.mkdir()
        started = time.monotonic()
        bundle = restored = None
        if policy["family"] == "gnomix":
            config = dict(policy, candidate=candidate, calibrate=False, retrain_base=False,
                          seed=native.SEED, sources=SOURCES)
            bundle = native.fit_chromosome(cohort, config, chrom, work, 6, 64, native.load_upstream())
            if bundle["keys"] != keys:
                raise ValueError("Gnomix changed the observed exact marker panel")
            with (work / "models" / f"chr{chrom}.pkl").open("rb") as stream:
                restored = pickle.load(stream)
        fit_seconds = time.monotonic() - started
        for case in cases:
            for stage in ("inherited_phase", "array_phase_only"):
                begin = time.monotonic()
                target = case[stage]
                prefix = work / f"{case['sample_id']}_{stage}"
                if policy["family"] == "flare":
                    reference, panel = references[policy["reference"]]
                    command = flare_command(cohort, target, reference, panel, prefix, policy)
                    if recovered:
                        saved = read(work / f"{prefix.name}.command.json")
                        if saved != command:
                            raise ValueError(f"Recovered FLARE command changed: {candidate}/{prefix.name}")
                    else:
                        shared.checked(command, work, prefix.name)
                    model = read_model(str(prefix) + ".model")
                    ids, values = physical_predictions(
                        str(prefix) + ".anc.vcf.gz", model["ancestries"], {str(chrom): span}, frame)
                    if ids != [case["sample_id"]]:
                        raise ValueError("FLARE changed the evaluation sample axis")
                    prediction = values[0]
                else:
                    haplotypes = native.load_haplotypes(target, [case["sample_id"]], chrom, keys)
                    before = native.ordered_probabilities(
                        bundle["model"], bundle["model"].base.predict_proba(haplotypes))
                    after = native.ordered_probabilities(
                        restored["model"], restored["model"].base.predict_proba(haplotypes))
                    np.testing.assert_allclose(before, after, rtol=0, atol=1e-10)
                    _, physical = native.aggregate_probabilities(after, bundle["layout"])
                    prediction = physical[0] / length
                if (not np.isfinite(prediction).all() or np.any(prediction < 0)
                        or not np.isclose(prediction.sum(), 1)):
                    raise ValueError("Invalid observed-only ancestry fractions")
                rows.append(dict(
                    candidate=candidate, mask=mask, chromosome=chrom, span_bp=length,
                    sample_id=case["sample_id"], stage=stage, track=case["track"],
                    family_id=case["family_id"], fit_seconds=fit_seconds,
                    inference_seconds=time.monotonic() - begin,
                    **{f"pred_{source}": float(prediction[index])
                       for index, source in enumerate(SOURCES)}))
        print(f"{mask} chr{chrom}: {candidate} verified", flush=True)
    predictions = output / "predictions.tsv"
    pd.DataFrame(rows).to_csv(predictions, sep="\t", index=False)
    native.write_json(output / "complete.json", dict(
        status="complete_observed_only_model_expansion_task", chromosome=chrom, mask=mask,
        source_inputs_sha256=native.digest(source_task / "inputs_complete.json"),
        protocol_sha256=native.digest(protocol_path), predictions_sha256=native.digest(predictions),
        imputation_performed=False))


def summarize(source, protocol_path, results, scope_chromosomes, masks, output):
    output.mkdir()
    protocol = read(protocol_path)
    frames = []
    expected = {(mask, chrom, candidate) for mask in masks for chrom in scope_chromosomes
                for candidate in protocol["candidates"]}
    actual = set()
    for result in results:
        record = read(result / "complete.json")
        if (record["status"] != "complete_observed_only_model_candidate"
                or record["protocol_sha256"] != native.digest(protocol_path)
                or record["predictions_sha256"] != native.digest(result / "predictions.tsv")
                or record["imputation_performed"]):
            raise ValueError("Invalid observed-only expansion task certificate")
        actual.add((record["mask"], record["chromosome"], record["candidate"]))
        frames.append(pd.read_csv(result / "predictions.tsv", sep="\t", float_precision="round_trip"))
    if actual != expected:
        raise ValueError(f"Incomplete observed-only expansion axes: {expected - actual}")
    chromosomes = pd.concat(frames, ignore_index=True)
    legacy = read(PROJECT / "results_v3_phasefree_array_validation_01/experiment.json")
    aggregate = common.aggregate(
        chromosomes, scope_chromosomes, protocol["candidates"],
        [row["sample_id"] for row in legacy["samples"]], masks)
    truth = {row["sample_id"]: row["truth"] for row in legacy["samples"]}
    targets = read(PROJECT / "config/phasefree_accuracy_targets.json")
    summaries = []
    for axes, group in aggregate.groupby(["candidate", "mask", "stage", "track"]):
        candidate, mask, stage, track = axes
        expected_values = np.asarray([truth[name] for name in group.sample_id])
        predictions = group[[f"pred_{source}" for source in SOURCES]].to_numpy(float)
        metrics, errors = accuracy_summary(expected_values, predictions, targets)
        summaries.append(dict(candidate=candidate, mask=mask, stage=stage, track=track, **metrics))
        aggregate.loc[group.index, "tv"] = errors["tv"]
    aggregate.to_csv(output / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(output / "comparison.tsv", sep="\t", index=False)
    chromosomes.to_csv(output / "chromosome_predictions.tsv", sep="\t", index=False)
    native.write_json(output / "completion.json", dict(
        status="complete_observed_only_96_model_expansion",
        configurations=len(protocol["candidates"]), imputation_performed=False,
        independent_final_test=False, production_promotion=False,
        output_sha256=shared.fingerprint_files([
            output / "predictions.tsv", output / "comparison.tsv",
            output / "chromosome_predictions.tsv"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "prepare", "evaluate", "evaluate-candidate", "summarize"))
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--source-task", type=Path)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--candidate")
    parser.add_argument("--results", type=Path, nargs="*")
    parser.add_argument("--chromosomes", type=int, nargs="+")
    parser.add_argument("--masks", nargs="+")
    parser.add_argument("--recovery", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.stage == "freeze":
            freeze(args.source, args.plan, args.output)
        elif args.stage == "prepare":
            prepare(args.source_task, args.protocol, args.output)
        elif args.stage == "evaluate":
            evaluate(args.source_task, args.protocol, args.output, args.recovery)
        elif args.stage == "evaluate-candidate":
            evaluate_candidate(
                args.prepared, args.protocol, args.candidate, args.output, args.recovery)
        else:
            summarize(
                args.source, args.protocol, args.results, args.chromosomes, args.masks, args.output)


if __name__ == "__main__":
    main()
