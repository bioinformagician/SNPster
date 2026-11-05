#!/usr/bin/env python3
"""Shared observed-only, phase-only native FLARE/Gnomix development benchmark."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import pickle
import subprocess
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as native
import pilot_matched_imputation as shared
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import physical_predictions, read_model
from validate_phasefree_arrays import write_truth_vcf


PROJECT = Path(__file__).resolve().parents[2]
SCREEN = Path("/results")
LEGACY = Path("/legacy")
PROTOCOL = PROJECT / "config/observed_ancestry_protocol.json"
BEAGLE = PROJECT / "results_v3_matched_global_01/beagle.jar"
SOURCES = list(native.SOURCES)


def read(path):
    return json.loads(Path(path).read_text())


def typed_keys(path):
    keys = []
    for line in Path(path).read_text().splitlines():
        chrom, pos, ref, alt = line.split("\t")
        key = (chrom.removeprefix("chr"), int(pos), ref, alt)
        if key[0] not in {str(c) for c in range(1, 23)} or len(ref) != 1 or len(alt) != 1:
            raise ValueError(f"Expected autosomal biallelic SNP key: {key}")
        keys.append(key)
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate exact observed marker")
    return set(keys)


def validate_phase(original, inferred, expected_ids, actual_ids, expected_keys, actual_keys):
    if (expected_ids != actual_ids or expected_keys != actual_keys
            or original.shape != inferred.shape or original.ndim != 2
            or original.shape[0] != 2 * len(expected_ids)
            or original.shape[1] != len(expected_keys)
            or not np.isin(original, [0, 1]).all()
            or not np.isin(inferred, [0, 1]).all()):
        raise ValueError("Phase-only processing changed axes or introduced missing/invalid alleles")
    a = original.reshape(-1, 2, original.shape[1]).sum(axis=1)
    b = inferred.reshape(-1, 2, inferred.shape[1]).sum(axis=1)
    if not np.array_equal(a, b):
        raise ValueError("Phase-only processing changed supplied diploid GTs")
    return int(a.size)


def phasing_command(target, reference, genetic_map, prefix):
    return ["java", "-Xmx4g", "-jar", str(BEAGLE), f"gt={target}", f"ref={reference}",
            f"map={genetic_map}", f"out={prefix}", "impute=false", "nthreads=6",
            f"seed={native.SEED}"]


def flare_command(cohort, target, prefix, policy):
    return ["java", "-Xmx8g", "-jar", cohort["flare_jar"], f"gt={target}",
            f"ref={cohort['reference']['vcf']}", f"ref-panel={cohort['reference']['broad_panel']}",
            f"map={cohort['genetic_map']}", f"out={prefix}",
            f"array={str(policy['array']).lower()}", f"em={str(policy['em']).lower()}",
            "update-p=false", "probs=true", "min-maf=0", "min-mac=0",
            "nthreads=6", f"seed={native.SEED}"]


def implementation():
    paths = [Path(__file__), PROTOCOL, BEAGLE, Path(native.__file__), Path(shared.__file__),
             Path(__file__).with_name("screen_model_families.py"),
             Path(__file__).with_name("validate_phasefree_arrays.py"),
             Path(__file__).with_name("report_phasefree_accuracy.py"),
             PROJECT / "config/phasefree_accuracy_targets.json"]
    return shared.fingerprint_files(paths)


def context():
    protocol, cohort, legacy = read(PROTOCOL), read(SCREEN / "cohort.json"), read(LEGACY / "experiment.json")
    if protocol["sources"] != SOURCES or cohort["sources"] != SOURCES:
        raise ValueError("Incompatible source axis")
    native.validate_splits(cohort)
    reference = native.sample_records(cohort["reference"]["labels_tsv"], True)
    train = native.sample_records(cohort["pure_train"]["labels_tsv"], True)
    ref_ids, train_ids = {r["sample_id"] for r in reference}, {r["sample_id"] for r in train}
    donors = set(legacy["donor_ids"])
    if ref_ids & train_ids or donors & (ref_ids | train_ids):
        raise ValueError("Reference/training/evaluation founder overlap")
    if len(legacy["samples"]) != 14 or len(ref_ids) != 3447 or len(train_ids) != 307:
        raise ValueError("Changed pilot/reference/training identities")
    return protocol, cohort, legacy


def panel_keys(cohort, observed, chrom):
    frame = pd.read_csv(cohort["markers_tsv"], sep="\t", header=None,
                        names=["chrom", "pos", "ref", "alt"], dtype={"chrom": str})
    keys = [(r.chrom, int(r.pos), r.ref, r.alt) for r in frame.itertuples(index=False)
            if r.chrom == str(chrom) and (r.chrom, int(r.pos), r.ref, r.alt) in observed]
    if len(keys) < 2 or len({key[1] for key in keys}) != len(keys):
        raise ValueError(f"Empty/single-marker or ambiguous panel on chr{chrom}")
    return keys


def complete(path, identity):
    if not path.exists():
        return False
    value = read(path)
    if value["identity"] != identity or value["status"] != "complete":
        raise ValueError(f"Incompatible checkpoint: {path}")
    shared.verify_hashes(value["output_sha256"])
    return True


def subset(source, ids, keys, directory, name):
    output = directory / f"{name}.vcf.gz"
    regions = directory / f"{name}.positions.tsv"
    regions.write_text("".join(f"{key[0]}\t{key[1]}\n" for key in keys))
    command = ["bcftools", "view", "-R", regions, "-Oz", "-o", output]
    if ids is not None:
        selected = directory / f"{name}.samples.txt"
        selected.write_text("".join(identity + "\n" for identity in ids))
        command += ["-S", selected]
    shared.checked(command + [source], directory, f"{name}_extract")
    shared.checked(["bcftools", "index", "-t", output], directory, f"{name}_index")
    actual = [(key, gt) for key, gt in native.variants(output, keys[0][0], ids)]
    if [key for key, _ in actual] != keys:
        raise ValueError(f"Observed-only extraction changed exact marker axis: {name}")
    if any(not np.isin(native.phased_gt(gt), [0, 1]).all()
           for _, values in actual for gt in values):
        raise ValueError(f"Missing/incomplete GT in {name}; no filling permitted")
    return str(output)


def compress(path, directory, name):
    output = Path(str(path) + ".gz")
    shared.checked(["bcftools", "view", "-Oz", "-o", output, path], directory, name)
    shared.checked(["bcftools", "index", "-t", output], directory, name + "_index")
    path.unlink()
    return str(output)


def build_inputs(out, mask, chrom, smoke=False):
    protocol, original, legacy = context()
    directory = out / f"{mask}_chr{chrom}"
    identity = dict(protocol_sha256=native.digest(PROTOCOL), mask=mask, chrom=chrom, smoke=smoke)
    marker = directory / "inputs_complete.json"
    if complete(marker, identity):
        return read(directory / "cohort.json"), read(directory / "cases.json")
    if directory.exists():
        raise ValueError(f"Inspect incomplete observed-only inputs rather than overwriting: {directory}")
    directory.mkdir(parents=True)
    observed = typed_keys(legacy["masks"][mask])
    keys = panel_keys(original, observed, chrom)
    cohort = json.loads(json.dumps(original))
    ref_ids = [r["sample_id"] for r in native.sample_records(original["reference"]["labels_tsv"], True)]
    cohort["reference"]["vcf"] = subset(original["reference"]["vcf"], ref_ids, keys, directory, "reference")
    train_ids = [r["sample_id"] for r in native.sample_records(original["pure_train"]["labels_tsv"], True)]
    cohort["pure_train"]["vcf"] = subset(original["pure_train"]["vcf"], train_ids, keys, directory, "pure_train")
    cohort.pop("markers_tsv")
    cohort["reference_filters"] = dict(min_maf=0, min_mac=0)
    for i, batch in enumerate(cohort["pedigree_train"]["batches"]):
        ids = [r["sample_id"] for r in native.sample_records(batch["samples_tsv"])]
        batch["vcf"] = subset(batch["vcf"], ids, keys, directory, f"pedigree_train_{i}")
    samples = legacy["samples"]
    if smoke:
        samples = [next(r for r in samples if r["track"] == track) for track in ("pure", "balanced")]
    cases = []
    for sample in samples:
        name = sample["sample_id"]
        work = directory / "cases" / name
        work.mkdir(parents=True)
        inherited = native.load_haplotypes(legacy["batches"][sample["batch"]]["vcf"], [name], chrom, keys)
        if not np.isin(inherited, [0, 1]).all():
            raise ValueError("Missing observed GT; imputation is forbidden")
        known = work / "inherited.vcf"
        write_truth_vcf(known, chrom, keys, [name], inherited.reshape(1, 2, -1))
        inherited_path = compress(known, work, "inherited")
        target = work / "typed.unphased.vcf"
        write_truth_vcf(target, chrom, keys, [name], inherited.reshape(1, 2, -1), unphased=True)
        unphased_path = compress(target, work, "typed")
        prefix = work / "phase_only"
        seconds = shared.checked(phasing_command(
            unphased_path, cohort["reference"]["vcf"], cohort["genetic_map"], prefix), work, "phase_only")
        raw_phase = str(prefix) + ".vcf.gz"
        contig_header = work / "phase.contig.hdr"
        contig_header.write_text(f"##contig=<ID={chrom}>\n")
        phased = str(prefix) + ".normalized.vcf.gz"
        shared.checked(["bcftools", "annotate", "-h", contig_header, "-Oz", "-o", phased, raw_phase],
                       work, "phase_only_contig")
        shared.checked(["bcftools", "index", "-t", phased], work, "phase_only_index")
        actual = list(native.variants(phased, chrom, [name]))
        inferred = np.asarray([native.phased_gt(values[0]) for _, values in actual], np.int8).T
        n = validate_phase(inherited, inferred, [name], list(native.vcf_header(phased)[0]),
                           keys, [key for key, _ in actual])
        phase_info = shared.phase_metrics(inherited, inferred)
        row = dict(sample_id=name, mask=mask, track="pure" if sample["track"] == "pure" else "pedigree",
                   family_id=sample["family_id"],
                   inherited_phase=inherited_path, array_phase_only=phased,
                   n_typed_gt_preserved=n, phase_seconds=seconds, phase_diagnostic=phase_info)
        native.write_json(work / "phase_verification.json", row)
        cases.append(row)
    native.write_json(directory / "cohort.json", cohort)
    native.write_json(directory / "cases.json", cases)
    pd.DataFrame(keys, columns=["chrom", "pos", "ref", "alt"]).to_csv(
        directory / "markers.tsv", sep="\t", index=False)
    native.write_json(directory / "coverage.json", dict(
        mask=mask, chromosome=chrom, observed_sites=sum(k[0] == str(chrom) for k in observed),
        retained_sites=len(keys), common_panel_limited=True, full_array_benchmark=False))
    outputs = [directory / name for name in ("cohort.json", "cases.json", "markers.tsv", "coverage.json")]
    outputs += list(directory.rglob("*.vcf.gz")) + list(directory.rglob("*verification.json"))
    native.write_json(marker, dict(status="complete", identity=identity,
                                  output_sha256=shared.fingerprint_files(outputs)))
    print(f"{mask} chr{chrom}: {len(keys)} typed sites; {len(cases)} singleton phases verified", flush=True)
    return cohort, cases


def evaluate_task(out, mask, chrom, smoke=False):
    protocol = read(PROTOCOL)
    cohort, cases = build_inputs(out, mask, chrom, smoke)
    directory = out / f"{mask}_chr{chrom}"
    keys_frame = pd.read_csv(directory / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(r.chrom, int(r.pos), r.ref, r.alt) for r in keys_frame.itertuples(index=False)]
    identity = dict(inputs_sha256=native.digest(directory / "inputs_complete.json"),
                    protocol_sha256=native.digest(PROTOCOL))
    finished = directory / "complete.json"
    if complete(finished, identity):
        return str(directory)
    if (directory / "predictions.tsv").exists():
        raise ValueError("Incomplete predictions require inspection")
    rows = []
    span = cohort["chromosome_spans"][str(chrom)]
    length = span[1] - span[0]
    for candidate, policy in protocol["candidates"].items():
        work = directory / candidate
        if work.exists():
            raise ValueError(f"Incomplete candidate requires inspection: {work}")
        work.mkdir()
        bundle = reloaded = None
        started = time.monotonic()
        if policy["family"] == "gnomix":
            config = dict(policy, candidate=candidate, calibrate=False, retrain_base=False,
                          seed=native.SEED, sources=SOURCES)
            bundle = native.fit_chromosome(cohort, config, chrom, work, 6, 64, native.load_upstream())
            if bundle["keys"] != keys:
                raise ValueError("Gnomix training changed the shared exact observed marker panel")
            with (work / "models" / f"chr{chrom}.pkl").open("rb") as stream:
                reloaded = pickle.load(stream)
        fit_seconds = time.monotonic() - started
        for case in cases:
            for stage in ("inherited_phase", "array_phase_only"):
                begin = time.monotonic()
                target = case[stage]
                prefix = work / f"{case['sample_id']}_{stage}"
                if policy["family"] == "flare":
                    shared.checked(flare_command(cohort, target, prefix, policy), work, prefix.name)
                    model = read_model(str(prefix) + ".model")
                    ids, values = physical_predictions(str(prefix) + ".anc.vcf.gz", model["ancestries"],
                                                      {str(chrom): span}, keys_frame)
                    if ids != [case["sample_id"]]:
                        raise ValueError("FLARE sample axis changed")
                    prediction = values[0]
                else:
                    x = native.load_haplotypes(target, [case["sample_id"]], chrom, keys)
                    before = native.ordered_probabilities(bundle["model"], bundle["model"].base.predict_proba(x))
                    after = native.ordered_probabilities(reloaded["model"], reloaded["model"].base.predict_proba(x))
                    np.testing.assert_allclose(before, after, rtol=0, atol=1e-10)
                    _, physical = native.aggregate_probabilities(after, bundle["layout"])
                    prediction = physical[0] / length
                if not np.isfinite(prediction).all() or np.any(prediction < 0) or not np.isclose(prediction.sum(), 1):
                    raise ValueError("Invalid native ancestry fractions")
                rows.append(dict(candidate=candidate, sample_id=case["sample_id"], mask=mask,
                                 stage=stage, track=case["track"], family_id=case["family_id"],
                                 chromosome=chrom, span_bp=length, fit_seconds=fit_seconds,
                                 inference_seconds=time.monotonic() - begin,
                                 **{f"pred_{s}": float(prediction[i]) for i, s in enumerate(SOURCES)}))
        print(f"{mask} chr{chrom}: {candidate} native predictions verified", flush=True)
    path = directory / "predictions.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    outputs = [path] + list(directory.glob("*/models/*")) + list(directory.glob("*/*.model"))
    native.write_json(finished, dict(status="complete", identity=identity,
                                    output_sha256=shared.fingerprint_files(outputs)))
    return str(directory)


def aggregate(frame, chromosomes, candidates, samples, masks):
    columns = [f"pred_{s}" for s in SOURCES]
    required = set(chromosomes)
    expected = {(candidate, mask, stage, sample)
                for candidate in candidates for mask in masks
                for stage in ("inherited_phase", "array_phase_only") for sample in samples}
    groups = frame.groupby(["candidate", "mask", "stage", "sample_id"], sort=True)
    if set(groups.groups) != expected:
        raise ValueError("Incomplete method/mask/phase/person coverage")
    rows = []
    for key, group in groups:
        if (len(group) != len(required) or set(group.chromosome) != required
                or group.chromosome.duplicated().any()):
            raise ValueError("Incomplete/duplicate chromosome coverage; partial scores forbidden")
        values = group[columns].to_numpy(float)
        if (not np.isfinite(values).all() or np.any(values < 0)
                or not np.allclose(values.sum(axis=1), 1, rtol=0, atol=1e-8)
                or np.any(group.span_bp <= 0)):
            raise ValueError("Invalid chromosome predictions or physical spans")
        prediction = np.average(values, axis=0, weights=group.span_bp)
        rows.append(dict(candidate=key[0], mask=key[1], stage=key[2], sample_id=key[3],
                         track=group.track.iloc[0], family_id=group.family_id.iloc[0],
                         **dict(zip(columns, prediction))))
    return pd.DataFrame(rows)


def summarize(out):
    protocol, _, legacy = context()
    paths = [out / f"{mask}_chr{chrom}" / "predictions.tsv"
             for mask in protocol["masks"] for chrom in protocol["chromosomes"]]
    for path in paths:
        value = read(path.parent / "complete.json")
        shared.verify_hashes(value["output_sha256"])
    frame = pd.concat([pd.read_csv(p, sep="\t", float_precision="round_trip") for p in paths])
    predictions = aggregate(frame, protocol["chromosomes"], protocol["candidates"],
                            [r["sample_id"] for r in legacy["samples"]], protocol["masks"])
    truth = {r["sample_id"]: r["truth"] for r in legacy["samples"]}
    targets = read(PROJECT / "config/phasefree_accuracy_targets.json")
    summary = []
    for (candidate, mask, stage, track), group in predictions.groupby(["candidate", "mask", "stage", "track"]):
        a = np.asarray([truth[s] for s in group.sample_id])
        b = group[[f"pred_{s}" for s in SOURCES]].to_numpy(float)
        metrics, errors = accuracy_summary(a, b, targets)
        summary.append(dict(candidate=candidate, mask=mask, stage=stage, track=track, **metrics,
                            off_target_mass=float((b * (a < .01)).sum(axis=1).mean()),
                            **{f"mae_{s}": float(np.abs(a - b)[:, i].mean()) for i, s in enumerate(SOURCES)}))
        predictions.loc[group.index, "tv"] = errors["tv"]
        predictions.loc[group.index, "worst_component_error"] = errors["worst_component_error"]
    predictions.to_csv(out / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out / "comparison.tsv", sep="\t", index=False)
    native.write_json(out / "completion.json", dict(
        status="complete_observed_only_development_pilot", independent_final_test=False,
        full_array_benchmark=False, imputation_performed=False, production_promotion=False,
        n_unique_people=14, n_masks=2, n_autosomes=22,
        output_sha256=shared.fingerprint_files([out / "predictions.tsv", out / "comparison.tsv"])))


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Use a fresh benchmark output")
    protocol, cohort, legacy = context()
    if read(out / "smoke" / "smoke_passed.json")["implementation_sha256"] != implementation():
        raise ValueError("Run the real smoke against current sources before freezing")
    inputs = [SCREEN / "cohort.json", LEGACY / "experiment.json", Path(cohort["reference"]["vcf"]),
              Path(cohort["pure_train"]["vcf"]), Path(cohort["markers_tsv"]),
              Path(cohort["reference"]["labels_tsv"]), Path(cohort["pure_train"]["labels_tsv"]),
              Path(cohort["reference"]["broad_panel"]), Path(cohort["flare_jar"]), Path(cohort["genetic_map"])]
    inputs += [Path(path) for path in legacy["masks"].values()]
    for batch in cohort["pedigree_train"]["batches"]:
        inputs += [Path(batch[key]) for key in ("vcf", "samples_tsv", "segments_tsv")]
    inputs += [Path(value["vcf"]) for value in legacy["batches"].values()]
    coverage = []
    for mask in protocol["masks"]:
        observed = typed_keys(legacy["masks"][mask])
        for chrom in protocol["chromosomes"]:
            keys = panel_keys(cohort, observed, chrom)
            for policy in protocol["candidates"].values():
                if policy["family"] == "gnomix":
                    native.window_layout(keys, native.map_points(cohort["genetic_map"], chrom),
                                         cohort["chromosome_spans"][str(chrom)], policy)
            coverage.append(dict(mask=mask, chromosome=chrom, retained_sites=len(keys),
                                 observed_sites=sum(k[0] == str(chrom) for k in observed)))
    pd.DataFrame(coverage).to_csv(out / "coverage.tsv", sep="\t", index=False)
    native.write_json(out / "experiment.json", dict(
        protocol=protocol, implementation_sha256=implementation(),
        input_sha256=shared.fingerprint_files(inputs), no_imputation=True,
        full_array_benchmark=False, independent_final_test=False, production_promotion=False))
    native.write_json(out / "state.json", dict(status="prepared", imputation_performed=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["smoke", "prepare", "run"], required=True)
    parser.add_argument("--outdir", type=Path, default=Path("/benchmark"))
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        if args.stage == "prepare":
            prepare(out)
        elif args.stage == "smoke":
            smoke = out / "smoke"
            if smoke.exists():
                raise ValueError("Inspect and retain an existing smoke; use a fresh output for a new version")
            for mask in read(PROTOCOL)["masks"]:
                evaluate_task(smoke, mask, 22, smoke=True)
            native.write_json(smoke / "smoke_passed.json", dict(
                status="complete_nonranking_smoke", imputation_performed=False,
                implementation_sha256=implementation()))
        else:
            config = read(out / "experiment.json")
            shared.verify_hashes(config["implementation_sha256"])
            shared.verify_hashes(config["input_sha256"])
            protocol = config["protocol"]
            native.write_json(out / "state.json", dict(status="running_observed_only", imputation_performed=False))
            with ProcessPoolExecutor(max_workers=2) as pool:
                tasks = [pool.submit(evaluate_task, out, mask, chrom)
                         for chrom in protocol["chromosomes"] for mask in protocol["masks"]]
                for index, task in enumerate(as_completed(tasks), 1):
                    print(f"Observed benchmark {index}/{len(tasks)} verified: {task.result()}", flush=True)
            summarize(out)
            native.write_json(out / "state.json", dict(
                status="complete_observed_only_development_pilot", production_promotion=False))


if __name__ == "__main__":
    main()
