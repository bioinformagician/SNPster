#!/usr/bin/env python3
"""Single-task implementations; scheduling and concurrency belong to Nextflow."""

import argparse
from collections import defaultdict
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_observed_ancestry as common
import model_family_gnomix as native
import pilot_matched_imputation as shared
import screen_flare_reference_density as parent
import verify_flare_reference_density as independent
from screen_model_families import physical_predictions, read_model, write_model
from verify_dense_observed import check_case, local_physical


PROJECT = Path(__file__).resolve().parents[2]
OLD = Path("/screen")
PROTOCOL = PROJECT / "config/flare_nextflow_followup.json"
SOURCES = tuple(native.SOURCES)
AXES = ["panel", "chromosome", "candidate", "stage", "sample_id"]


def crossed_theta(template, donor):
    if (set(template["ancestries"]) != set(donor["ancestries"])
            or set(template["panels"]) != set(donor["panels"])
            or set(template["ancestries"]) != set(SOURCES)):
        raise ValueError("Theta transfer requires matching named seven-source axes")
    result = dict(template)
    i = [donor["ancestries"].index(name) for name in template["ancestries"]]
    j = [donor["panels"].index(name) for name in template["panels"]]
    result["theta"] = donor["theta"][np.ix_(i, j)].copy()
    return result


def current_cohort(panel, chrom, reference):
    cohort = common.read(Path("/benchmark/tasks") / f"{panel}_chr{chrom}/cohort.json")
    return {**cohort, "reference": dict(
        vcf=(cohort["reference"]["vcf"] if reference == "full" else
             str(OLD / f"tasks/{panel}_chr{chrom}/{reference}/reference.vcf.gz")),
        broad_panel=str(OLD / "inputs" / f"{reference}.panel"))}


def row_key(row):
    return tuple(row[key] for key in AXES)


def build_jobs(selection, old_rows, protocol, smoke=False):
    people = (parent.smoke_people(selection, 20261018) if smoke else selection["people"])
    existing = {row_key(row) for row in old_rows}
    jobs = []
    for person, candidate, stage in itertools.product(
            people, protocol["recovery"]["candidates"], protocol["stages"]):
        key = ("frederik", 1, candidate, stage, person["sample_id"])
        if key in existing:
            raise ValueError("Recovery would rerun an already verified prediction")
        jobs.append(dict(kind="recovery", panel="frederik", chromosome=1, candidate=candidate,
                         reference=candidate.rsplit("_", 1)[0], em=candidate.endswith("_em"),
                         stage=stage, sample_id=person["sample_id"]))
    for panel, person, reference, stage in itertools.product(
            protocol["theta_controls"]["panels"], people,
            protocol["theta_controls"]["references"], protocol["stages"]):
        jobs.append(dict(kind="theta", panel=panel, chromosome=12,
                         candidate=f"{reference}_crossed_theta_em", reference=reference,
                         em=True, stage=stage, sample_id=person["sample_id"]))
    for job in jobs:
        job["id"] = "_".join(str(job[k]) for k in ("kind", "panel", "chromosome", "candidate", "stage", "sample_id"))
    if len({job["id"] for job in jobs}) != len(jobs):
        raise ValueError("Duplicate task identities")
    return jobs


def prepare(output, mode, source_manifest):
    output.mkdir()
    (output / "models").mkdir()
    source_hashes = common.read(source_manifest)
    shared.verify_hashes(source_hashes)
    original = common.read(OLD / "experiment.json")
    shared.verify_hashes(original["input_sha256"])
    shared.verify_hashes(original["implementation_sha256"])
    audit = common.read(OLD / "report_20261021/verification.json")
    if audit["n_audited_predictions"] != 31600 or audit["complete_sweep"]:
        raise ValueError("Expected explicit independently audited partial failure")
    shared.verify_hashes(audit["output_sha256"])
    selection, protocol = common.read(OLD / "selection.json"), common.read(PROTOCOL)
    pure = [row for row in selection["people"] if row["track"] == "pure"]
    reference_labels = native.records("/benchmark/inputs/reference_labels.tsv")
    forbidden = set(common.read("/benchmark/selection.json")["evaluation_donors"])
    forbidden.update(row["sample_id"] for row in reference_labels)
    if (len(selection["people"]) != 316 or len(pure) != 240 or
            forbidden & {row["sample_id"] for row in selection["people"]}):
        raise ValueError("Training-only cohort contract violated")
    original_rows = pd.read_csv(OLD / "report_20261021/audited_predictions.tsv",
                               sep="\t", float_precision="round_trip")
    if len(original_rows) != 31600 or original_rows.duplicated(AXES).any():
        raise ValueError("Prior audit row coverage changed")
    truth_count, truth_hashes = independent.verify_truth(selection, original["protocol"])
    jobs = build_jobs(selection, original_rows.to_dict("records"), protocol, mode == "smoke")
    metadata = {row["sample_id"]: row for row in selection["people"]}
    inputs = [OLD / "experiment.json", OLD / "selection.json",
              OLD / "report_20261021/verification.json", Path(source_manifest), PROTOCOL]
    inputs += [Path(path) for path in audit["output_sha256"]]
    blocks = {(job["panel"], job["chromosome"]) for job in jobs}
    model_deltas, genotype_count = {}, 0
    for panel, chrom in sorted(blocks):
        task = OLD / f"tasks/{panel}_chr{chrom}"
        cohort = common.read(Path("/benchmark/tasks") / f"{panel}_chr{chrom}/cohort.json")
        keys = parent.markers(task)
        cases = common.read(task / "cases.json")
        desired = {job["sample_id"] for job in jobs if (job["panel"], job["chromosome"]) == (panel, chrom)}
        paths = defaultdict(list)
        batches = {b["batch"]: b["vcf"] for b in cohort["pedigree_train"]["batches"]}
        for identity in desired:
            person = metadata[identity]
            path = cohort["pure_train"]["vcf"] if person["track"] == "pure" else batches[person["batch"]]
            paths[path].append(identity)
        original_calls = {}
        for path, identities in paths.items():
            original_calls.update(independent.parent_calls(path, sorted(identities), keys))
        for case in cases:
            if case["sample_id"] not in desired:
                continue
            original_check = task / "cases" / case["sample_id"] / "verified.json"
            shared.verify_hashes(common.read(original_check)["output_sha256"])
            n, _, _ = check_case((case, keys))
            genotype_count += n
            inherited = list(independent.calls(case["inherited_phase"]))[1:]
            if any(original_calls[case["sample_id"]][i] != 2 * int(gt[1][1][0]) + int(gt[1][1][2])
                   for i, gt in enumerate(inherited)):
                raise ValueError("Reused inherited phase differs from original person")
            inputs += [Path(case["inherited_phase"]), Path(case["array_phase_only"]), original_check]
        for reference in ("full", "finecap16", "size_control"):
            template = task / reference / "initialization.fixed.model"
            donor_reference = "finecap16" if reference == "full" else "full"
            donor = task / donor_reference / "initialization.fixed.model"
            model = read_model(template)
            if chrom == 12:
                updated = crossed_theta(model, read_model(donor))
                name = f"{panel}_{chrom}_{reference}_crossed.model"
                write_model(output / "models" / name, updated)
                restored = read_model(output / "models" / name)
                for field in ("p", "rho", "mu", "generation"):
                    np.testing.assert_array_equal(restored[field], model[field])
                np.testing.assert_array_equal(restored["theta"], updated["theta"])
                model_deltas[f"{panel}_{reference}"] = dict(
                    native_theta=float(model["theta"][0, 0]), crossed_theta=float(updated["theta"][0, 0]))
            current = current_cohort(panel, chrom, reference)
            inputs += [template, donor, Path(current["reference"]["vcf"]),
                       Path(current["reference"]["broad_panel"]), Path(current["flare_jar"]),
                       Path(current["genetic_map"]), task / "markers.tsv", task / "cases.json"]
        for job in jobs:
            if (job["panel"], job["chromosome"]) != (panel, chrom):
                continue
            job["case"] = next(case for case in cases if case["sample_id"] == job["sample_id"])
            job["model"] = (str(task / job["reference"] / "initialization.fixed.model")
                            if job["kind"] == "recovery" else
                            f"models/{panel}_{chrom}_{job['reference']}_crossed.model")
    characterization_task = OLD / "tasks/frederik_chr1"
    inputs += list((characterization_task / "finecap16").glob("initialization*.model"))
    inputs += [characterization_task / "finecap16/initialization.anc.vcf.gz"]
    original_rows.to_csv(output / "reused_predictions.tsv", sep="\t", index=False)
    native.write_json(output / "jobs.json", jobs)
    native.write_json(output / "selection.json", selection)
    native.write_json(output / "experiment.json", dict(
        mode=mode, protocol=protocol, n_jobs=len(jobs), expected_new_jobs=6320 if mode == "full" else 40,
        genotype_calls_reverified=genotype_count, raw_bp_truth_vectors=truth_count,
        raw_bp_plan_sha256=truth_hashes, theta_contrasts=model_deltas,
        source_sha256=source_hashes, input_sha256=shared.fingerprint_files(inputs)))
    if len(jobs) != (6320 if mode == "full" else 40):
        raise ValueError("Job coverage differs from locked design")


def characterize(prepared, output):
    output.mkdir()
    task = OLD / "tasks/frederik_chr1"
    keys = parent.markers(task)
    cohort = current_cohort("frederik", 1, "finecap16")
    span = cohort["chromosome_spans"]["1"]
    case = common.read(task / "cases.json")[0]
    implicit = read_model(task / "finecap16/initialization.model")
    explicit = read_model(task / "finecap16/initialization.fixed.model")
    for field in ("p", "rho", "theta", "generation"):
        np.testing.assert_array_equal(implicit[field], explicit[field])
    values = {}
    policies = dict(implicit_repeat=None,
                    explicit_repeat1=task / "finecap16/initialization.fixed.model",
                    explicit_repeat2=task / "finecap16/initialization.fixed.model",
                    rounded_prior=task / "finecap16/initialization.model")
    for name, model_path in policies.items():
        prefix = output / name
        shared.checked(parent.command(cohort, case["inherited_phase"], prefix, False, model_path), output, name)
        model = read_model(str(prefix) + ".model")
        values[name] = local_physical(str(prefix) + ".anc.vcf.gz", model["ancestries"], span, keys)
    values["original_implicit"] = local_physical(
        task / "finecap16/initialization.anc.vcf.gz", implicit["ancestries"], span, keys)
    deltas = {f"{a}_vs_{b}": float(np.max(np.abs(values[a] - values[b])))
              for a, b in itertools.combinations(values, 2)}
    repeat_error = float(np.max(np.abs(values["explicit_repeat1"] - values["explicit_repeat2"])))
    if repeat_error > 1e-10 or max(deltas.values()) > 1e-6:
        raise ValueError(f"Initialization not reproducible within the new declared numeric contract: {deltas}")
    native.write_json(output / "verification.json", dict(
        status="characterized_with_explicit_numeric_contract", original_assertion_unchanged=True,
        original_tolerance=1e-10, new_max_initialization_fraction_delta=1e-6,
        independent_posterior_aggregation_tolerance=1e-10, repeat_explicit_error=repeat_error,
        native_mu=implicit["mu"].tolist(), explicit_mu=explicit["mu"].tolist(), measured_deltas=deltas,
        exact_numerical_cause_proven=False,
        prediction_vectors={key: value.tolist() for key, value in values.items()}))


def infer(prepared, characterization, identity, output):
    check = common.read(characterization / "verification.json")
    if check["status"] != "characterized_with_explicit_numeric_contract":
        raise ValueError("Characterization gate did not pass")
    job = next(row for row in common.read(prepared / "jobs.json") if row["id"] == identity)
    selection = common.read(prepared / "selection.json")
    person = next(row for row in selection["people"] if row["sample_id"] == job["sample_id"])
    output.mkdir()
    task = OLD / f"tasks/{job['panel']}_chr{job['chromosome']}"
    keys = parent.markers(task)
    frame = pd.DataFrame(keys, columns=["chrom", "pos", "ref", "alt"])
    cohort = current_cohort(job["panel"], job["chromosome"], job["reference"])
    span = cohort["chromosome_spans"][str(job["chromosome"])]
    model_path = Path(job["model"]) if job["kind"] == "recovery" else prepared / job["model"]
    prefix = output / "native"
    target = Path(job["case"][job["stage"]])
    input_hashes = shared.fingerprint_files([target, model_path, Path(cohort["reference"]["vcf"]),
                                            Path(cohort["reference"]["broad_panel"])])
    command = parent.command(cohort, target, prefix, job["em"], model_path)
    native.write_json(output / "command.json", command)
    seconds = shared.checked(command, output, "inference")
    model = read_model(str(prefix) + ".model")
    ids, p = physical_predictions(str(prefix) + ".anc.vcf.gz", model["ancestries"],
                                   {str(job["chromosome"]): span}, frame)
    if ids != [job["sample_id"]]:
        raise ValueError("Native inference sample axis changed")
    independently = local_physical(str(prefix) + ".anc.vcf.gz", model["ancestries"], span, keys)
    np.testing.assert_allclose(p[0], independently, rtol=0, atol=1e-10)
    if not np.isfinite(p).all() or (p < 0).any() or abs(p.sum() - 1) > 1e-10:
        raise ValueError("Invalid native fractions")
    truth = np.array(selection["truth"][str(job["chromosome"])][job["sample_id"]])
    error = np.abs(p[0] - truth)
    row = dict(panel=job["panel"], chromosome=job["chromosome"], candidate=job["candidate"],
               stage=job["stage"], sample_id=job["sample_id"], kind=job["kind"],
               track=person["track"], family_id=person["family_id"], fine=person["fine"],
               source_group=person["source_group"], scenario=person["scenario"],
               n_markers=len(keys), span_bp=span[1] - span[0], seconds=seconds,
               tv=float(error.sum() / 2), max_component_error=float(error.max()),
               off_target_mass=float(p[0][truth == 0].sum()),
               **{f"pred_{s}": float(p[0, i]) for i, s in enumerate(SOURCES)},
               **{f"truth_{s}": float(truth[i]) for i, s in enumerate(SOURCES)})
    native.write_json(output / "prediction.json", dict(
        row=row, input_sha256=input_hashes, independent_local_posterior_aggregation=True,
        output_sha256=shared.fingerprint_files([
            Path(str(prefix) + ".anc.vcf.gz"), Path(str(prefix) + ".model")])))


def metric_tables(frame, output, name):
    summaries = []
    axes = ["panel", "chromosome", "candidate", "stage", "track"]
    for key, group in frame.groupby(axes):
        summaries.append(dict(zip(axes, key), n_people=len(group), mean_tv=group.tv.mean(),
                              p95_tv=group.tv.quantile(.95), max_tv=group.tv.max(),
                              coverage10=(group.max_component_error <= .1 + 1e-12).mean(),
                              coverage5=(group.max_component_error <= .05 + 1e-12).mean(),
                              family_balanced_tv=group.groupby("family_id").tv.mean().mean(),
                              off_target_mass=group.off_target_mass.mean()))
    frame.to_csv(output / f"{name}_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(output / f"{name}_comparison.tsv", sep="\t", index=False)
    frame[frame.track.eq("pure")].groupby(axes[:-1] + ["source_group"]).agg(
        n_people=("sample_id", "size"), mean_tv=("tv", "mean"),
        p95_tv=("tv", lambda values: values.quantile(.95))).reset_index().to_csv(
            output / f"{name}_source_errors.tsv", sep="\t", index=False)


def aggregate(prepared, characterization, output):
    jobs = common.read(prepared / "jobs.json")
    records, spot_checks = [], 0
    selection = common.read(prepared / "selection.json")
    mode = common.read(prepared / "experiment.json")["mode"]
    for job in jobs:
        path = Path(job["id"])
        record = common.read(path / "prediction.json")
        shared.verify_hashes(record["output_sha256"])
        row = record["row"]
        if row_key(row) != row_key(job):
            raise ValueError("Per-task output identity differs from frozen job")
        p = np.asarray([row[f"pred_{s}"] for s in SOURCES])
        truth = np.array(selection["truth"][str(row["chromosome"])][row["sample_id"]])
        np.testing.assert_allclose(truth, [row[f"truth_{s}"] for s in SOURCES], rtol=0, atol=1e-12)
        if abs(np.abs(p - truth).sum() / 2 - row["tv"]) > 1e-12:
            raise ValueError("Incorrect persisted per-task TV")
        if not np.isfinite(p).all() or (p < 0).any() or abs(p.sum() - 1) > 1e-10:
            raise ValueError("Invalid persisted per-task fractions")
        model = read_model(path / "native.model")
        keys = parent.markers(OLD / f"tasks/{row['panel']}_chr{row['chromosome']}")
        cohort = current_cohort(row["panel"], row["chromosome"], job["reference"])
        value = local_physical(path / "native.anc.vcf.gz", model["ancestries"],
                               cohort["chromosome_spans"][str(row["chromosome"])], keys)
        np.testing.assert_allclose(p, value, rtol=0, atol=1e-10)
        records.append(row)
        spot_checks += 1
    output.mkdir()
    new = pd.DataFrame(records)
    if len(new) != len(jobs) or new.duplicated(AXES).any():
        raise ValueError("Incomplete/duplicate Nextflow task coverage")
    old = pd.read_csv(prepared / "reused_predictions.tsv", sep="\t", float_precision="round_trip")
    recovered = pd.concat([old, new[new.kind.eq("recovery")]], ignore_index=True)
    if recovered.duplicated(AXES).any() or (mode == "full" and len(recovered) != 34128):
        raise ValueError("Recovered sweep coverage differs from the original Cartesian design")
    if mode == "full":
        expected = set(itertools.product(
            ("lonni", "frederik", "frederik_thin"), (22, 12, 1),
            [name for name, _, _ in parent.candidates(common.read(OLD / "experiment.json")["protocol"])],
            ("inherited_phase", "array_phase_only"), [r["sample_id"] for r in selection["people"]]))
        if {row_key(row) for row in recovered.to_dict("records")} != expected:
            raise ValueError("Recovered global prediction axes incomplete")
    native_baselines = old[old.chromosome.eq(12) & old.panel.isin(["lonni", "frederik"])
                           & old.candidate.isin(["full_em", "finecap16_em", "size_control_em"])]
    if mode == "smoke":
        native_baselines = native_baselines[native_baselines.sample_id.isin(new.sample_id)]
    crossed = pd.concat([native_baselines, new[new.kind.eq("theta")]], ignore_index=True)
    metric_tables(recovered, output, "recovered_screen")
    metric_tables(crossed, output, "crossed_theta")
    native.write_json(output / "verification.json", dict(
        status="independently_reaggregated_smoke" if mode == "smoke" else "independently_reaggregated_complete_followup",
        mode=mode, new_tasks=len(new), independent_local_reaggregations=spot_checks,
        n_recovered_screen_predictions=len(recovered), expected_full_screen_predictions=34128,
        n_crossed_theta_predictions=int(new.kind.eq("theta").sum()),
        original_failure_preserved=True, original_assertion_unchanged=True,
        whole_genome_score=False, independent_final_test=False,
        characterization=common.read(characterization / "verification.json"),
        input_sha256=shared.fingerprint_files([
            prepared / "experiment.json", prepared / "jobs.json", prepared / "reused_predictions.tsv"]),
        output_sha256=shared.fingerprint_files(list(output.glob("*.tsv")))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "characterize", "infer", "aggregate"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--characterization", type=Path)
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--mode", choices=("smoke", "full"))
    parser.add_argument("--job")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.stage == "prepare":
            prepare(args.output, args.mode, args.source_manifest)
        elif args.stage == "characterize":
            characterize(args.prepared, args.output)
        elif args.stage == "infer":
            infer(args.prepared, args.characterization, args.job, args.output)
        else:
            aggregate(args.prepared, args.characterization, args.output)
