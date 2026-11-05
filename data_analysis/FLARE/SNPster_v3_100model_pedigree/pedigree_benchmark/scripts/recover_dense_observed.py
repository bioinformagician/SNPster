#!/usr/bin/env python3
"""Explicit recovery of the single native-JVM failure; frozen producers unchanged."""

import argparse
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_dense_observed as dense
import benchmark_observed_ancestry as common
import model_family_gnomix as native
import pilot_matched_imputation as shared
from screen_model_families import physical_predictions, read_model


FAILED = ("flare_nonarray_fixed", "fam0000031_g2-b1-i1", "array_phase_only")


def output_prefix(work, candidate, sample, stage):
    prefix = work / f"{sample}_{stage}"
    if (candidate, sample, stage) == FAILED:
        return prefix.with_name(prefix.name + "_native_retry")
    return prefix


def recover(out):
    experiment = common.read(out / "experiment.json")
    shared.verify_hashes(experiment["implementation_sha256"])
    shared.verify_hashes(experiment["input_sha256"])
    directory = out / "tasks/lonni_chr4"
    inputs = common.read(directory / "inputs_complete.json")
    shared.verify_hashes(inputs["output_sha256"])
    expected = dict(protocol_sha256=native.digest(dense.PROTOCOL),
                    selection_sha256=native.digest(out / "selection.json"),
                    chromosome=4, mask="lonni", smoke=False)
    if inputs["identity"] != expected or (directory / "complete.json").exists():
        raise ValueError("Recovery requires exactly the inspected incomplete Lonni chr4 task")
    missing = []
    for chrom in experiment["protocol"]["chromosomes"]:
        for mask in experiment["protocol"]["masks"]:
            task = out / "tasks" / f"{mask}_chr{chrom}"
            marker = task / "complete.json"
            if not marker.exists():
                missing.append(task.name)
                continue
            identity = dict(inputs_sha256=native.digest(task / "inputs_complete.json"),
                            protocol_sha256=native.digest(dense.PROTOCOL))
            common.complete(marker, identity)
    if missing != ["lonni_chr4"]:
        raise ValueError(f"Recovery scope changed: {missing}")
    crash = out / "hs_err_pid122181.log"
    if "StubRoutines::jbyte_disjoint_arraycopy" not in crash.read_text():
        raise ValueError("Expected retained native JVM crash evidence")
    annex = out / "recovery.json"
    if annex.exists():
        raise ValueError("Inspect previous recovery before retrying")
    native.write_json(annex, dict(
        reason="Native Java SIGSEGV in arraycopy during BGZIP input, not an ancestry assertion",
        failed_task="lonni_chr4", original_state=common.read(out / "state.json"),
        original_experiment_sha256=native.digest(out / "experiment.json"),
        implementation_sha256=shared.fingerprint_files([Path(__file__)]),
        crash_sha256=native.digest(crash), reuse_complete_tasks=43,
        scientific_parameters_changed=False, imputation_performed=False,
        failed_command_policy="Exact same command and seed; distinct output prefix preserves failed files"))
    cohort, cases = common.read(directory / "cohort.json"), common.read(directory / "cases.json")
    frame = pd.read_csv(directory / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(r.chrom, int(r.pos), r.ref, r.alt) for r in frame.itertuples(index=False)]
    span = cohort["chromosome_spans"]["4"]
    length, rows, reused, executed = span[1] - span[0], [], 0, 0
    for candidate, policy in experiment["protocol"]["candidates"].items():
        work = directory / candidate
        bundle = restored = None
        started = time.monotonic()
        if policy["family"] == "gnomix":
            work.mkdir()
            config = dict(policy, candidate=candidate, calibrate=False, retrain_base=False,
                          seed=native.SEED, sources=dense.SOURCES)
            bundle = native.fit_chromosome(cohort, config, 4, work, 6, 64, native.load_upstream())
            if bundle["keys"] != keys:
                raise ValueError("Recovered Gnomix marker axis differs")
            with (work / "models/chr4.pkl").open("rb") as stream:
                restored = pickle.load(stream)
        fit_seconds = time.monotonic() - started
        for case in cases:
            for stage in ("inherited_phase", "array_phase_only"):
                begin = time.monotonic()
                target = case[stage]
                if policy["family"] == "flare":
                    prefix = output_prefix(work, candidate, case["sample_id"], stage)
                    anc, model_file = Path(str(prefix) + ".anc.vcf.gz"), Path(str(prefix) + ".model")
                    if anc.exists() and model_file.exists():
                        old = common.read(work / f"{prefix.name}.command.json")
                        if old != common.flare_command(cohort, target, prefix, policy):
                            raise ValueError("Saved FLARE command differs from fixed protocol")
                        reused += 1
                    elif anc.exists() or model_file.exists():
                        raise ValueError(f"Unexpected partial result: {prefix}")
                    else:
                        shared.checked(common.flare_command(cohort, target, prefix, policy), work, prefix.name)
                        executed += 1
                    model = read_model(model_file)
                    ids, p = physical_predictions(anc, model["ancestries"], {"4": span}, frame)
                    if ids != [case["sample_id"]]:
                        raise ValueError("Recovered FLARE sample axis differs")
                    prediction = p[0]
                else:
                    x = native.load_haplotypes(target, [case["sample_id"]], 4, keys)
                    a = native.ordered_probabilities(bundle["model"], bundle["model"].base.predict_proba(x))
                    b = native.ordered_probabilities(restored["model"], restored["model"].base.predict_proba(x))
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)
                    _, physical = native.aggregate_probabilities(b, bundle["layout"])
                    prediction = physical[0] / length
                if (not np.isfinite(prediction).all() or np.any(prediction < 0)
                        or not np.isclose(prediction.sum(), 1)):
                    raise ValueError("Invalid recovered ancestry fractions")
                rows.append(dict(candidate=candidate, mask="lonni", chromosome=4, span_bp=length,
                                 sample_id=case["sample_id"], stage=stage, track=case["track"],
                                 family_id=case["family_id"], fit_seconds=fit_seconds,
                                 inference_seconds=time.monotonic() - begin,
                                 **{f"pred_{s}": float(prediction[i]) for i, s in enumerate(dense.SOURCES)}))
        print(f"Recovered lonni chr4 {candidate}; preserved outputs={reused}, new FLARE calls={executed}", flush=True)
    pd.DataFrame(rows).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    outputs = [directory / "predictions.tsv"] + list(directory.glob("*/models/*")) + list(directory.glob("*/*.model"))
    native.write_json(directory / "complete.json", dict(
        status="complete", identity=dict(inputs_sha256=native.digest(directory / "inputs_complete.json"),
                                        protocol_sha256=native.digest(dense.PROTOCOL)),
        recovery_sha256=native.digest(annex), output_sha256=shared.fingerprint_files(outputs)))
    dense.summarize(out)
    native.write_json(out / "recovery_state.json", dict(
        status="complete_native_observed_only_pilot_after_explicit_recovery",
        complete_tasks=44, reused_tasks=43, reused_flare_calls=reused, executed_flare_calls=executed,
        failed_native_call_retried_once=True, original_failed_state_preserved=True,
        imputation_performed=False, production_promotion=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("/benchmark"))
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        recover(args.outdir)
