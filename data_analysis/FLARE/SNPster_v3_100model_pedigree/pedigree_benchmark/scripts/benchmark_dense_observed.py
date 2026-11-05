#!/usr/bin/env python3
"""Native FLARE/Gnomix benchmark on shared dense-supported observed SNP panels."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import pickle
import re
import subprocess
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_observed_ancestry as common
import model_family_gnomix as native
import pilot_matched_imputation as shared
from realized_ancestry import load_bp
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import physical_predictions, read_model
from validate_phasefree_arrays import write_truth_vcf


PROJECT = Path(__file__).resolve().parents[2]
PROTOCOL = PROJECT / "config/observed_dense_protocol.json"
TRAIN_ROOT = PROJECT / "results_v3_pedigree_training_screen_01/cohort/pedigree"
EVAL_ROOT = PROJECT / "results_v3_pedigree_pilot_10/pedigree"
SOURCES = list(native.SOURCES)


def dense_path(chrom):
    return Path("/dense/beagle_references") / shared.REFERENCE_PATTERN.format(chrom=chrom)


def implementation():
    paths = [Path(__file__), PROTOCOL, Path(common.__file__), Path(native.__file__),
             Path(shared.__file__), common.BEAGLE, Path(__file__).with_name("realized_ancestry.py"),
             Path(__file__).with_name("screen_model_families.py"),
             Path(__file__).with_name("validate_phasefree_arrays.py"),
             Path(__file__).with_name("report_phasefree_accuracy.py"),
             PROJECT / "config/phasefree_accuracy_targets.json",
             PROJECT / "run_dense_observed_benchmark.sh"]
    return shared.fingerprint_files(paths)


def safe_families(donors, available):
    return {family for family, identities in donors.items() if set(identities) <= set(available)}


def decode_calls(text, n_people):
    raw = np.frombuffer((text + "\t").encode("ascii"), dtype=np.uint8)
    if len(raw) != 4 * n_people:
        raise ValueError("Require complete phased biallelic diploid GTs; missing filling is forbidden")
    values = raw.reshape(n_people, 4)
    if (not np.isin(values[:, [0, 2]], [48, 49]).all()
            or np.any(values[:, 1] != ord("|")) or np.any(values[:, 3] != ord("\t"))):
        raise ValueError("Invalid/missing/unphased dense GT; no imputation or coercion")
    return (values[:, [0, 2]] - 48).astype(np.int8)


def select_inputs(out):
    target = out / "selection.json"
    if target.exists():
        return common.read(target)
    protocol, original, legacy = common.read(PROTOCOL), common.read(common.SCREEN / "cohort.json"), common.read(common.LEGACY / "experiment.json")
    native.validate_splits(original)
    available = set(native.vcf_header(str(dense_path(22)))[0])
    reference = [r for r in native.sample_records(original["reference"]["labels_tsv"], True)
                 if r["sample_id"] in available]
    pure = [r for r in native.sample_records(original["pure_train"]["labels_tsv"], True)
            if r["sample_id"] in available]
    ref_ids, train_ids = [r["sample_id"] for r in reference], [r["sample_id"] for r in pure]
    evaluation = set(legacy["donor_ids"])
    if (len(reference) != 3033 or len(pure) != 240 or not evaluation <= available
            or set(ref_ids) & set(train_ids) or evaluation & (set(ref_ids) | set(train_ids))
            or any({r["source_group"] for r in rows} != set(SOURCES) for rows in (reference, pure))):
        raise ValueError("Approved reference/training/source/donor contract changed")
    donor_file = PROJECT / "results_v3_phasefree_rich_01/donors.json"
    donors = common.read(donor_file)
    safe = safe_families(donors, train_ids)
    directory = out / "inputs"
    directory.mkdir()
    reference_labels, pure_labels = directory / "reference_labels.tsv", directory / "pure_train_labels.tsv"
    pd.DataFrame(reference).to_csv(reference_labels, sep="\t", index=False)
    pd.DataFrame(pure).to_csv(pure_labels, sep="\t", index=False)
    panel = directory / "reference.panel"
    pd.DataFrame(reference)[["sample_id", "source_group"]].to_csv(panel, sep="\t", header=False, index=False)
    batches = []
    for batch in original["pedigree_train"]["batches"]:
        rows = [r for r in native.records(batch["samples_tsv"]) if r["family_id"] in safe]
        if not rows:
            continue
        samples = directory / f"{batch['batch']}.samples.tsv"
        pd.DataFrame(rows).to_csv(samples, sep="\t", index=False)
        batches.append(dict(batch=batch["batch"], samples_tsv=str(samples),
                            segments_tsv=batch["segments_tsv"], original_vcf=batch["vcf"]))
    result = dict(reference_ids=ref_ids, pure_ids=train_ids, evaluation_donors=sorted(evaluation),
                  reference_labels=str(reference_labels), pure_labels=str(pure_labels), panel=str(panel),
                  n_safe_training_families=len(safe), n_original_training_families=len(donors),
                  n_training_descendants=sum(len(native.records(b["samples_tsv"])) for b in batches),
                  pedigree_train=batches, selection_uses_errors=False,
                  source_counts={name: pd.DataFrame(rows).source_group.value_counts().to_dict()
                                 for name, rows in (("reference", reference), ("pure_train", pure))})
    native.write_json(target, result)
    return result


def extract_dense(source, ids, wanted, directory):
    regions = directory / "positions.tsv"
    regions.write_text("".join(f"{chrom}\t{pos}\n"
                              for chrom, pos in sorted({(k[0], k[1]) for k in wanted}, key=lambda k: k[1])))
    command = ["bcftools", "query", "-H", "-R", str(regions), "-s", ",".join(ids),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(source)]
    native.write_json(directory / "dense_query.command.json", command)
    keys, columns, seen, reorder = [], [], set(), None
    with (directory / "dense_query.log").open("w") as error:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error, text=True)
        try:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t", 4)
                if line.startswith("#"):
                    names = [re.sub(r"^\[\d+\]", "", field).removesuffix(":GT")
                             for field in fields[4].split("\t")]
                    if len(names) != len(ids) or set(names) != set(ids):
                        raise ValueError("Dense query sample-axis mismatch")
                    reorder = np.asarray([names.index(identity) for identity in ids])
                    continue
                if reorder is None:
                    raise ValueError("Dense query lacks explicit sample header")
                key = (fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3])
                if key not in wanted:
                    continue
                if key in seen:
                    raise ValueError("Duplicate exact dense marker")
                seen.add(key)
                columns.append(decode_calls(fields[4], len(ids))[reorder].reshape(-1))
                keys.append(key)
            if process.wait():
                raise ValueError(f"Dense query failed; inspect {directory / 'dense_query.log'}")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    if not keys or len({k[1] for k in keys}) != len(keys):
        raise ValueError("Empty or ambiguous exact typed marker panel")
    return keys, np.stack(columns, axis=1)


def write_haplotypes(directory, name, chrom, keys, ids, haplotypes, unphased=False):
    path = directory / f"{name}.vcf"
    write_truth_vcf(path, chrom, keys, ids, haplotypes.reshape(len(ids), 2, -1), unphased=unphased)
    return common.compress(path, directory, name)


def original_concordance(path, ids, haplotypes, chrom, keys, canonical):
    indices = [i for i, key in enumerate(keys) if key in canonical]
    if not indices:
        raise ValueError("No legacy common-marker anchors for dense genotype validation")
    expected = native.load_haplotypes(path, ids, chrom, [keys[i] for i in indices])
    if not np.array_equal(expected, haplotypes[:, indices]):
        raise ValueError("Dense haplotypes do not reproduce original canonical genotypes/phase")
    return len(ids) * len(indices)


def phase_case(work, name, chrom, keys, inherited, cohort):
    inherited_path = write_haplotypes(work, "inherited", chrom, keys, [name], inherited)
    typed = write_haplotypes(work, "typed.unphased", chrom, keys, [name], inherited, unphased=True)
    prefix = work / "phase_only"
    seconds = shared.checked(common.phasing_command(
        typed, cohort["reference"]["vcf"], cohort["genetic_map"], prefix), work, "phase_only")
    header = work / "contig.hdr"
    header.write_text(f"##contig=<ID={chrom}>\n")
    phased = str(prefix) + ".normalized.vcf.gz"
    shared.checked(["bcftools", "annotate", "-h", header, "-Oz", "-o", phased, str(prefix) + ".vcf.gz"],
                   work, "phase_header")
    shared.checked(["bcftools", "index", "-t", phased], work, "phase_index")
    rows = list(native.variants(phased, chrom, [name]))
    inferred = np.asarray([native.phased_gt(values[0]) for _, values in rows], np.int8).T
    n = common.validate_phase(inherited, inferred, [name], list(native.vcf_header(phased)[0]),
                              keys, [key for key, _ in rows])
    result = dict(sample_id=name, inherited_phase=inherited_path, array_phase_only=phased,
                  n_typed_gt_preserved=n, phase_seconds=seconds,
                  phase_diagnostic=shared.phase_metrics(inherited, inferred))
    native.write_json(work / "phase_verification.json", result)
    return result


def build_inputs(out, mask, chrom, smoke):
    selection = select_inputs(out)
    original, legacy = common.read(common.SCREEN / "cohort.json"), common.read(common.LEGACY / "experiment.json")
    protocol = common.read(PROTOCOL)
    directory = out / ("smoke" if smoke else "tasks") / f"{mask}_chr{chrom}"
    identity = dict(protocol_sha256=native.digest(PROTOCOL), selection_sha256=native.digest(out / "selection.json"),
                    chromosome=chrom, mask=mask, smoke=smoke)
    if common.complete(directory / "inputs_complete.json", identity):
        return common.read(directory / "cohort.json"), common.read(directory / "cases.json")
    if directory.exists():
        raise ValueError(f"Inspect incomplete dense preparation before retry: {directory}")
    directory.mkdir(parents=True)
    span = original["chromosome_spans"][str(chrom)]
    points = native.map_points(original["genetic_map"], chrom)
    observed = common.typed_keys(legacy["masks"][mask])
    wanted = {k for k in observed if k[0] == str(chrom) and span[0] <= k[1] < span[1]
              and points[0, 0] <= k[1] <= points[-1, 0]}
    source = dense_path(chrom)
    expected = legacy["reference_sha256"][str(Path("/imputer/beagle_references") / source.name)]
    if native.digest(source) != expected:
        raise ValueError("Existing dense resource changed")
    ids = selection["reference_ids"] + selection["pure_ids"] + selection["evaluation_donors"]
    keys, haplotypes = extract_dense(source, ids, wanted, directory)
    nref, ntrain = len(selection["reference_ids"]), len(selection["pure_ids"])
    counts = haplotypes[:2 * nref].sum(axis=0, dtype=np.int64)
    mac = np.minimum(counts, 2 * nref - counts)
    keep = (mac >= protocol["reference_filters"]["min_mac"]) & (
        mac / (2 * nref) >= protocol["reference_filters"]["min_maf"])
    n_available = len(keys)
    keys = [key for key, retain in zip(keys, keep) if retain]
    haplotypes = haplotypes[:, keep]
    if len(keys) < 2:
        raise ValueError("Insufficient shared quality-filtered typed markers")
    canonical_frame = pd.read_csv(original["markers_tsv"], sep="\t", header=None,
                                  names=["chrom", "pos", "ref", "alt"], dtype={"chrom": str})
    canonical = {(r.chrom, int(r.pos), r.ref, r.alt) for r in canonical_frame.itertuples(index=False)}
    cohort = dict(sources=SOURCES, seed=native.SEED, genetic_map=original["genetic_map"],
                  flare_jar=original["flare_jar"], chromosome_spans=original["chromosome_spans"],
                  reference_filters=dict(min_maf=0, min_mac=0),
                  reference=dict(labels_tsv=selection["reference_labels"], broad_panel=selection["panel"]),
                  pure_train=dict(labels_tsv=selection["pure_labels"]), pedigree_train=dict(batches=[]))
    cohort["reference"]["vcf"] = write_haplotypes(directory, "reference", chrom, keys,
                                                selection["reference_ids"], haplotypes[:2 * nref])
    pure_haps = haplotypes[2 * nref:2 * (nref + ntrain)]
    cohort["pure_train"]["vcf"] = write_haplotypes(directory, "pure_train", chrom, keys,
                                                 selection["pure_ids"], pure_haps)
    anchors = original_concordance(original["reference"]["vcf"], selection["reference_ids"],
                                   haplotypes[:2 * nref], chrom, keys, canonical)
    anchors += original_concordance(original["pure_train"]["vcf"], selection["pure_ids"],
                                    pure_haps, chrom, keys, canonical)
    positions = np.asarray([key[1] for key in keys])
    for batch in selection["pedigree_train"]:
        plan = common.read(TRAIN_ROOT / f"plans/plans/{batch['batch']}/pedigrees.json")
        bp = load_bp(TRAIN_ROOT / f"chromosomes/pdb_{batch['batch']}_chr{chrom}.bp", chrom)
        names = [r["sample_id"] for r in native.records(batch["samples_tsv"])]
        values = np.vstack([shared.pedigree_haplotypes(
            plan, bp, name, selection["pure_ids"], pure_haps, positions) for name in names])
        anchors += original_concordance(batch["original_vcf"], names, values, chrom, keys, canonical)
        path = write_haplotypes(directory, batch["batch"], chrom, keys, names, values)
        cohort["pedigree_train"]["batches"].append(dict(
            batch=batch["batch"], vcf=path, samples_tsv=batch["samples_tsv"], segments_tsv=batch["segments_tsv"]))
    evaluation_haps = haplotypes[2 * (nref + ntrain):]
    samples = legacy["samples"]
    if smoke:
        samples = [next(r for r in samples if r["track"] == track) for track in ("pure", "balanced")]
    cases = []
    for sample in samples:
        name = sample["sample_id"]
        if sample["track"] == "pure":
            offset = 2 * selection["evaluation_donors"].index(name)
            values = evaluation_haps[offset:offset + 2]
        else:
            plan = common.read(legacy["batches"][sample["batch"]]["plan"])
            bp = load_bp(EVAL_ROOT / f"chromosomes/pdb_{sample['batch']}_chr{chrom}.bp", chrom)
            values = shared.pedigree_haplotypes(
                plan, bp, name, selection["evaluation_donors"], evaluation_haps, positions)
        anchors += original_concordance(legacy["batches"][sample["batch"]]["vcf"],
                                        [name], values, chrom, keys, canonical)
        work = directory / "cases" / name
        work.mkdir(parents=True)
        case = phase_case(work, name, chrom, keys, values, cohort)
        case.update(mask=mask, family_id=sample["family_id"],
                    track="pure" if sample["track"] == "pure" else "pedigree")
        cases.append(case)
    pd.DataFrame(keys, columns=["chrom", "pos", "ref", "alt"]).to_csv(directory / "markers.tsv", sep="\t", index=False)
    native.write_json(directory / "cohort.json", cohort)
    native.write_json(directory / "cases.json", cases)
    native.write_json(directory / "coverage.json", dict(
        mask=mask, chromosome=chrom, observed_sites=sum(k[0] == str(chrom) for k in observed),
        in_span_and_map=len(wanted), available_exact_dense=n_available, retained_sites=len(keys),
        reference_quality_excluded=n_available - len(keys), canonical_gt_phase_anchors_verified=anchors,
        artificially_capped_to_old_common_panel=False))
    outputs = list(directory.rglob("*.vcf.gz")) + list(directory.rglob("*verification.json"))
    outputs += [directory / name for name in ("markers.tsv", "cohort.json", "cases.json", "coverage.json")]
    native.write_json(directory / "inputs_complete.json", dict(
        status="complete", identity=identity, output_sha256=shared.fingerprint_files(outputs)))
    print(f"{mask} chr{chrom}: {len(keys)}/{len(wanted)} typed sites; {anchors} original anchors; "
          f"{len(cases)} phase-only cases verified", flush=True)
    return cohort, cases


def evaluate(out, mask, chrom, smoke=False):
    cohort, cases = build_inputs(out, mask, chrom, smoke)
    directory = out / ("smoke" if smoke else "tasks") / f"{mask}_chr{chrom}"
    identity = dict(inputs_sha256=native.digest(directory / "inputs_complete.json"),
                    protocol_sha256=native.digest(PROTOCOL))
    if common.complete(directory / "complete.json", identity):
        return str(directory)
    keys_frame = pd.read_csv(directory / "markers.tsv", sep="\t", dtype={"chrom": str})
    keys = [(r.chrom, int(r.pos), r.ref, r.alt) for r in keys_frame.itertuples(index=False)]
    span = cohort["chromosome_spans"][str(chrom)]
    length, rows = span[1] - span[0], []
    for candidate, policy in common.read(PROTOCOL)["candidates"].items():
        work = directory / candidate
        if work.exists():
            raise ValueError(f"Inspect incomplete candidate before retry: {work}")
        work.mkdir()
        begin = time.monotonic()
        bundle = restored = None
        if policy["family"] == "gnomix":
            config = dict(policy, candidate=candidate, calibrate=False, retrain_base=False,
                          seed=native.SEED, sources=SOURCES)
            bundle = native.fit_chromosome(cohort, config, chrom, work, 6, 64, native.load_upstream())
            if bundle["keys"] != keys:
                raise ValueError("Gnomix changed the common exact marker panel")
            with (work / "models" / f"chr{chrom}.pkl").open("rb") as stream:
                restored = pickle.load(stream)
        fit_seconds = time.monotonic() - begin
        for case in cases:
            for stage in ("inherited_phase", "array_phase_only"):
                begin = time.monotonic()
                target = case[stage]
                prefix = work / f"{case['sample_id']}_{stage}"
                if policy["family"] == "flare":
                    shared.checked(common.flare_command(cohort, target, prefix, policy), work, prefix.name)
                    model = read_model(str(prefix) + ".model")
                    actual_ids, p = physical_predictions(
                        str(prefix) + ".anc.vcf.gz", model["ancestries"], {str(chrom): span}, keys_frame)
                    if actual_ids != [case["sample_id"]]:
                        raise ValueError("FLARE changed the sample axis")
                    prediction = p[0]
                else:
                    x = native.load_haplotypes(target, [case["sample_id"]], chrom, keys)
                    a = native.ordered_probabilities(bundle["model"], bundle["model"].base.predict_proba(x))
                    b = native.ordered_probabilities(restored["model"], restored["model"].base.predict_proba(x))
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)
                    _, physical = native.aggregate_probabilities(b, bundle["layout"])
                    prediction = physical[0] / length
                if not np.isfinite(prediction).all() or np.any(prediction < 0) or not np.isclose(prediction.sum(), 1):
                    raise ValueError("Invalid native ancestry fractions")
                rows.append(dict(candidate=candidate, mask=mask, chromosome=chrom, span_bp=length,
                                 sample_id=case["sample_id"], stage=stage, track=case["track"],
                                 family_id=case["family_id"], fit_seconds=fit_seconds,
                                 inference_seconds=time.monotonic() - begin,
                                 **{f"pred_{source}": float(prediction[i]) for i, source in enumerate(SOURCES)}))
        print(f"{mask} chr{chrom}: {candidate} verified", flush=True)
    pd.DataFrame(rows).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    outputs = [directory / "predictions.tsv"] + list(directory.glob("*/models/*")) + list(directory.glob("*/*.model"))
    native.write_json(directory / "complete.json", dict(
        status="complete", identity=identity, output_sha256=shared.fingerprint_files(outputs)))
    return str(directory)


def summarize(out):
    protocol = common.read(PROTOCOL)
    legacy = common.read(common.LEGACY / "experiment.json")
    frames, coverage = [], []
    for chrom in protocol["chromosomes"]:
        for mask in protocol["masks"]:
            directory = out / "tasks" / f"{mask}_chr{chrom}"
            shared.verify_hashes(common.read(directory / "complete.json")["output_sha256"])
            shared.verify_hashes(common.read(directory / "inputs_complete.json")["output_sha256"])
            frames.append(pd.read_csv(directory / "predictions.tsv", sep="\t", float_precision="round_trip"))
            coverage.append(common.read(directory / "coverage.json"))
    chromosomes = pd.concat(frames, ignore_index=True)
    frame = common.aggregate(chromosomes, protocol["chromosomes"], protocol["candidates"],
                             [r["sample_id"] for r in legacy["samples"]], protocol["masks"])
    truth = {row["sample_id"]: row["truth"] for row in legacy["samples"]}
    targets = common.read(PROJECT / "config/phasefree_accuracy_targets.json")
    summaries = []
    for (candidate, mask, stage, track), group in frame.groupby(["candidate", "mask", "stage", "track"]):
        a = np.asarray([truth[name] for name in group.sample_id])
        b = group[[f"pred_{s}" for s in SOURCES]].to_numpy(float)
        metrics, errors = accuracy_summary(a, b, targets)
        summaries.append(dict(candidate=candidate, mask=mask, stage=stage, track=track, **metrics,
                              off_target_mass=float((b * (a < .01)).sum(axis=1).mean()),
                              **{f"mae_{s}": float(np.abs(a - b)[:, i].mean()) for i, s in enumerate(SOURCES)}))
        frame.loc[group.index, "tv"] = errors["tv"]
        frame.loc[group.index, "worst_component_error"] = errors["worst_component_error"]
    frame.to_csv(out / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(out / "comparison.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(out / "coverage.tsv", sep="\t", index=False)
    chromosomes.to_csv(out / "chromosome_predictions.tsv", sep="\t", index=False)
    native.write_json(out / "completion.json", dict(
        status="complete_native_observed_only_pilot", imputation_performed=False,
        independent_final_test=False, production_promotion=False, n_unique_people=14,
        n_masks=2, n_autosomes=22, old_common_panel_cap=False,
        output_sha256=shared.fingerprint_files([out / name for name in (
            "predictions.tsv", "comparison.tsv", "coverage.tsv", "chromosome_predictions.tsv")])))


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Fresh experiment required; frozen producers must not be overwritten")
    if common.read(out / "smoke_passed.json")["implementation_sha256"] != implementation():
        raise ValueError("Real smoke does not match current implementation")
    original, legacy = common.read(common.SCREEN / "cohort.json"), common.read(common.LEGACY / "experiment.json")
    selection = select_inputs(out)
    paths = [common.SCREEN / "cohort.json", common.LEGACY / "experiment.json", out / "selection.json",
             out / "image.txt", PROJECT / "results_v3_phasefree_rich_01/donors.json",
             Path(original["flare_jar"]), Path(original["genetic_map"]), Path(original["markers_tsv"])]
    paths += list((out / "inputs").glob("*"))
    paths += [Path(path) for path in legacy["masks"].values()]
    paths += [Path(original["reference"][key]) for key in ("labels_tsv", "vcf")]
    paths += [Path(original["pure_train"][key]) for key in ("labels_tsv", "vcf")]
    paths += [Path(b["original_vcf"]) for b in selection["pedigree_train"]]
    paths += [Path(b["segments_tsv"]) for b in selection["pedigree_train"]]
    paths += [Path(b["vcf"]) for b in legacy["batches"].values()]
    for root, batches in ((TRAIN_ROOT, [b["batch"] for b in selection["pedigree_train"]]),
                          (EVAL_ROOT, [name for name in legacy["batches"] if name != "pure"])):
        for batch in batches:
            paths.append(root / f"plans/plans/{batch}/pedigrees.json")
            paths += [root / f"chromosomes/pdb_{batch}_chr{chrom}.bp" for chrom in range(1, 23)]
    required = set(selection["reference_ids"] + selection["pure_ids"] + selection["evaluation_donors"])
    references = {}
    for chrom in range(1, 23):
        path = dense_path(chrom)
        if not required <= set(native.vcf_header(str(path))[0]):
            raise ValueError(f"Approved identities unavailable on chromosome {chrom}")
        references[str(path)] = legacy["reference_sha256"][str(Path("/imputer/beagle_references") / path.name)]
    native.write_json(out / "experiment.json", dict(
        protocol=common.read(PROTOCOL), input_sha256=shared.fingerprint_files(paths),
        implementation_sha256=implementation(), reference_sha256=references,
        imputation_performed=False, independent_final_test=False, production_promotion=False))
    native.write_json(out / "state.json", dict(status="prepared_observed_only"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["smoke", "prepare", "run"], required=True)
    parser.add_argument("--outdir", type=Path, default=Path("/benchmark"))
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        if args.stage == "smoke":
            if (out / "smoke").exists():
                raise ValueError("Retain old smoke and use a fresh output for a new version")
            select_inputs(out)
            for mask in common.read(PROTOCOL)["masks"]:
                evaluate(out, mask, 22, smoke=True)
            native.write_json(out / "smoke_passed.json", dict(
                status="complete_nonranking_observed_only_smoke", implementation_sha256=implementation()))
        elif args.stage == "prepare":
            prepare(out)
        else:
            try:
                config = common.read(out / "experiment.json")
                shared.verify_hashes(config["implementation_sha256"])
                shared.verify_hashes(config["input_sha256"])
                native.write_json(out / "state.json", dict(status="running_observed_only", imputation_performed=False))
                with ProcessPoolExecutor(max_workers=2) as pool:
                    jobs = [pool.submit(evaluate, out, mask, chrom)
                            for chrom in config["protocol"]["chromosomes"] for mask in config["protocol"]["masks"]]
                    for index, job in enumerate(as_completed(jobs), 1):
                        print(f"Dense observed benchmark {index}/{len(jobs)}: {job.result()}", flush=True)
                summarize(out)
                native.write_json(out / "state.json", dict(
                    status="complete_native_observed_only_pilot", production_promotion=False))
            except Exception as error:
                native.write_json(out / "state.json", dict(status="failed", error=str(error)))
                raise


if __name__ == "__main__":
    main()
