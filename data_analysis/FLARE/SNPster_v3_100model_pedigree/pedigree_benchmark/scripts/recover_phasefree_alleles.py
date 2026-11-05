#!/usr/bin/env python3
"""Explicit allele-key recovery for the frozen array run; preserve every model/observed variant."""

import argparse
import json
from pathlib import Path
import pickle
import subprocess
import sys

import numpy as np
import pandas as pd

import schedule_phasefree_arrays as scheduler
import validate_phasefree_arrays as base
from verify_phasefree_arrays import one_sample_vcf


def load_exact(path, sample, keys):
    if not keys or len(set(keys)) != len(keys):
        raise ValueError("Exact-allele extraction requires unique nonempty variant keys")
    wanted = set(keys)
    found = {key: value for key, value in one_sample_vcf(path, sample) if key in wanted}
    if set(found) != wanted:
        raise ValueError(f"Missing {len(wanted - set(found))} exact-allele genotypes")
    return np.asarray([found[key] for key in keys], dtype=np.int8)


def annex(out):
    path = out / "allele_recovery.json"
    record = json.loads(path.read_text())
    if record["experiment_sha256"] != base.adapter.digest(out / "experiment.json"):
        raise ValueError("Allele recovery belongs to another frozen experiment")
    base.shared.verify_hashes(record["implementation_sha256"])
    return record


def prepare(out):
    config = base.configuration(out)
    path = out / "allele_recovery.json"
    if path.exists():
        annex(out)
    else:
        archive = scheduler.archive_incomplete_preparation(out, 4)
        base.adapter.write_json(path, dict(
            status="configured_exact_allele_recovery",
            experiment_sha256=base.adapter.digest(out / "experiment.json"),
            reason="chr4:164203107 has required model G/A and observed Frederik G/T; both exist in dense reference",
            policy="Keep both exact alleles; never drop, replace or reinterpret a model or supplied observed genotype",
            archive=archive, unchanged_scientific_inputs=True,
            implementation_sha256=base.shared.fingerprint_files([
                Path(__file__), base.PROJECT / "tests/test_phasefree_alleles.py",
                base.PROJECT / "run_phasefree_allele_recovery.sh"])))
    prepare_chromosome(out, config)


def prepare_chromosome(out, config):
    chrom = 4
    directory = out / "chr4"
    experiment_hash = base.adapter.digest(out / "experiment.json")
    if base.verified(directory / "prepared.json", experiment_hash):
        return
    if directory.exists():
        raise ValueError("Inspect incomplete exact-allele preparation before retrying")
    base.shared.verify_hashes(config["input_sha256"])
    cohort = json.loads((base.SCREEN / "cohort.json").read_text())
    with base.MODEL.open("rb") as stream:
        keys = [key for key in pickle.load(stream)["marker_keys"] if key[0] == str(chrom)]
    ref = base.REFERENCE / "beagle_references" / base.shared.REFERENCE_PATTERN.format(chrom=chrom)
    map_path = base.REFERENCE / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map"
    for path in (ref, Path(str(ref) + ".tbi"), map_path):
        if base.adapter.digest(path) != config["reference_sha256"][str(path)]:
            raise ValueError(f"Frozen reference/map changed: {path}")
    directory.mkdir()
    span = cohort["chromosome_spans"][str(chrom)]
    observed = {
        name: {(row[0], int(row[1]), row[2], row[3])
               for line in Path(path).read_text().splitlines()
               if (row := line.split("\t"))[0] == str(chrom) and span[0] <= int(row[1]) < span[1]}
        for name, path in config["masks"].items()}
    for name, mask in observed.items():
        if len({key[:2] for key in mask}) != len(mask):
            raise ValueError(f"Ambiguous observed alleles within a single array mask: {name}")
    wanted = set(keys).union(*observed.values())
    regions = directory / "positions.tsv"
    regions.write_text("".join(f"{chrom}\t{pos}\n" for pos in sorted({key[1] for key in wanted})))
    query = subprocess.check_output([
        "bcftools", "query", "-R", str(regions), "-s", ",".join(config["donor_ids"]),
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n", str(ref)], text=True)
    found, donor_values, seen = [], [], set()
    for line in query.splitlines():
        fields = line.split("\t")
        key = (fields[0], int(fields[1]), fields[2], fields[3])
        if key not in wanted:
            continue
        if key in seen:
            raise ValueError("Duplicate exact dense-reference variant")
        seen.add(key)
        values = dict(field.split("=") for field in fields[4:])
        found.append(key)
        donor_values.append([allele for sample in config["donor_ids"]
                             for allele in base.adapter.phased_gt(values[sample])])
    if not set(keys) <= seen:
        raise ValueError("Missing frozen model alleles")
    donors = np.asarray(donor_values, dtype=np.int8).T
    if not np.isin(donors, [0, 1]).all():
        raise ValueError("Incomplete dense founder haplotypes")
    positions, index = np.asarray([key[1] for key in found]), {key: i for i, key in enumerate(found)}
    model_indices = [index[key] for key in keys]
    names, truth = [sample["sample_id"] for sample in config["samples"]], []
    for sample in config["samples"]:
        batch = config["batches"][sample["batch"]]
        if sample["track"] == "pure":
            offset = 2 * config["donor_ids"].index(sample["sample_id"])
            values = donors[offset:offset + 2]
        else:
            plan = json.loads(Path(batch["plan"]).read_text())
            bp = base.load_bp(base.PEDIGREE / f"chromosomes/pdb_{sample['batch']}_chr{chrom}.bp", chrom)
            values = base.shared.pedigree_haplotypes(
                plan, bp, sample["sample_id"], config["donor_ids"], donors, positions)
        original = base.adapter.load_haplotypes(batch["vcf"], [sample["sample_id"]], chrom, keys)
        np.testing.assert_array_equal(original, values[:, model_indices])
        truth.append(values)
    truth = np.asarray(truth)
    truth_path = directory / "truth.vcf"
    base.write_truth_vcf(truth_path, chrom, found, names, truth)
    base.shared.checked(["bcftools", "view", "-Oz", "-o", directory / "truth.vcf.gz", truth_path],
                        directory, "compress_truth")
    base.shared.checked(["bcftools", "index", "-f", "-t", directory / "truth.vcf.gz"], directory, "index_truth")
    truth_path.unlink()
    pd.DataFrame(found, columns=["chrom", "pos", "ref", "alt"]).to_csv(
        directory / "truth_markers.tsv", sep="\t", index=False)
    np.save(directory / "truth_haplotypes.npy", truth)
    original_ids = subprocess.check_output(["bcftools", "query", "-l", str(ref)], text=True).splitlines()
    keep = [sample for sample in original_ids if sample not in set(config["excluded_reference_ids"])]
    if len(keep) != config["n_dense_reference_after_exclusion"] or set(keep) & set(config["donor_ids"]):
        raise ValueError("Incorrect recovery reference exclusions")
    keep_path, disjoint = directory / "reference.keep.txt", directory / "reference.disjoint.vcf.gz"
    keep_path.write_text("".join(sample + "\n" for sample in keep))
    base.shared.checked(["bcftools", "view", "-S", keep_path, "-Oz", "--threads", "2",
                         "-o", disjoint, ref], directory, "prepare_reference")
    actual = subprocess.check_output(["bcftools", "query", "-l", str(disjoint)], text=True).splitlines()
    if actual != keep:
        raise ValueError("Recovery imputation-reference identities changed")
    inputs = base.shared.fingerprint_files([
        ref, map_path, base.PREVIOUS / "beagle.jar", out / "allele_recovery.json",
        directory / "truth.vcf.gz", directory / "truth_markers.tsv", directory / "truth_haplotypes.npy"])
    base.adapter.write_json(directory / "prepared.json", dict(
        status="complete", experiment_sha256=experiment_hash, exact_allele_recovery=True,
        n_model_markers=len(keys), n_reconstructed_variants=len(found),
        reference_cache_sha256=base.adapter.digest(disjoint), input_sha256=inputs,
        output_sha256=base.shared.fingerprint_files([
            keep_path, directory / "truth.vcf.gz", directory / "truth.vcf.gz.tbi",
            directory / "truth_markers.tsv", directory / "truth_haplotypes.npy"])))
    print("chr4 exact-allele preparation preserves all model and available observed variants.", flush=True)


def case(out, info):
    config = base.configuration(out)
    annex(out)
    experiment_hash = base.adapter.digest(out / "experiment.json")
    directory, work = out / "chr4", out / "chr4" / info["case_id"]
    if base.verified(work / "completion.json", experiment_hash):
        return np.load(work / "dosages.npy")
    prepared = json.loads((directory / "prepared.json").read_text())
    disjoint = directory / "reference.disjoint.vcf.gz"
    if base.adapter.digest(disjoint) != prepared["reference_cache_sha256"]:
        raise ValueError("Changed recovered imputation reference")
    if work.exists():
        attempt = 1
        archive = work.with_name(f"{work.name}.incomplete_{attempt}")
        while archive.exists():
            attempt += 1
            archive = work.with_name(f"{work.name}.incomplete_{attempt}")
        work.rename(archive)
    work.mkdir()
    found = [(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
             for row in base.adapter.records(directory / "truth_markers.tsv")]
    seen, index = set(found), {key: i for i, key in enumerate(found)}
    truth = np.load(directory / "truth_haplotypes.npy")
    names = [sample["sample_id"] for sample in config["samples"]]
    sample_index = names.index(info["sample_id"])
    span = json.loads((base.SCREEN / "cohort.json").read_text())["chromosome_spans"]["4"]
    mask = {(row[0], int(row[1]), row[2], row[3])
            for line in Path(config["masks"][info["mask"]]).read_text().splitlines()
            if (row := line.split("\t"))[0] == "4" and span[0] <= int(row[1]) < span[1]}
    with base.MODEL.open("rb") as stream:
        keys = [key for key in pickle.load(stream)["marker_keys"] if key[0] == "4"]
    array_indices = [i for i, key in enumerate(found) if key in mask]
    if not array_indices:
        raise ValueError("No available observed recovery genotypes")
    array_keys = [found[i] for i in array_indices]
    target = work / "array.unphased.vcf"
    base.write_truth_vcf(target, 4, array_keys, [info["sample_id"]],
                         truth[sample_index:sample_index + 1, :, array_indices], unphased=True)
    prefix = work / "beagle"
    seconds = base.shared.checked([
        "java", "-Xmx16g", "-jar", base.PREVIOUS / "beagle.jar", f"gt={target}", f"ref={disjoint}",
        f"map={base.REFERENCE / 'plink.GRCh38.map/plink.chr4.GRCh38.map'}", f"out={prefix}",
        "nthreads=6", f"seed={config['seed']}", "gp=true"], work, "impute")
    imputed, retained = Path(str(prefix) + ".vcf.gz"), work / "retained.vcf.gz"
    base.shared.checked(["bcftools", "index", "-t", imputed], work, "index_imputed")
    base.shared.checked(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz", "-o",
                         retained, imputed, directory / "truth.vcf.gz"], work, "retain")
    base.shared.checked(["bcftools", "index", "-t", retained], work, "index_retained")
    np.testing.assert_array_equal(load_exact(retained, info["sample_id"], array_keys),
                                  truth[sample_index][:, array_indices].sum(axis=0))
    values = load_exact(retained, info["sample_id"], keys)
    np.save(work / "dosages.npy", values)
    model_truth = truth[sample_index][:, [index[key] for key in keys]].sum(axis=0)
    base.adapter.write_json(work / "completion.json", dict(
        status="complete", experiment_sha256=experiment_hash, seconds=seconds,
        exact_allele_recovery=True, observed_genotypes_preserved=True,
        n_observed=len(array_indices), n_unavailable_observed=len(mask - seen),
        n_model_markers=len(keys), n_genotype_errors=int(np.count_nonzero(values != model_truth)),
        input_sha256=base.shared.fingerprint_files([
            directory / "prepared.json", out / "allele_recovery.json", base.PREVIOUS / "beagle.jar", target]),
        output_sha256=base.shared.fingerprint_files([retained, work / "dosages.npy"])))
    imputed.unlink()
    Path(str(imputed) + ".tbi").unlink()
    print(f"chr4 {info['case_id']}: exact observed/model alleles verified in {seconds:.1f}s", flush=True)
    return values


def chromosome(out):
    config = base.configuration(out)
    if base.verified(out / "chr4/completion.json", base.adapter.digest(out / "experiment.json")):
        return
    prepare(out)
    values = [case(out, info) for info in config["cases"]]
    directory = out / "chr4"
    found = [(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
             for row in base.adapter.records(directory / "truth_markers.tsv")]
    index = {key: i for i, key in enumerate(found)}
    with base.MODEL.open("rb") as stream:
        keys = [key for key in pickle.load(stream)["marker_keys"] if key[0] == "4"]
    truth = np.load(directory / "truth_haplotypes.npy")
    np.savez_compressed(directory / "dosages.npz",
                        original=truth[:, :, [index[key] for key in keys]].sum(axis=1),
                        array_imputed=np.asarray(values))
    base.adapter.write_json(directory / "marker_keys.json", keys)
    prepared = json.loads((directory / "prepared.json").read_text())
    base.adapter.write_json(directory / "completion.json", dict(
        status="complete", experiment_sha256=base.adapter.digest(out / "experiment.json"),
        exact_allele_recovery=True, n_cases=28, original_haplotypes_reproduced=True,
        n_reference_after=config["n_dense_reference_after_exclusion"], input_sha256=prepared["input_sha256"],
        output_sha256=base.shared.fingerprint_files([
            directory / "dosages.npz", directory / "marker_keys.json",
            *[directory / info["case_id"] / "completion.json" for info in config["cases"]]])))
    (directory / "reference.disjoint.vcf.gz").unlink()


def run(out):
    annex(out)
    lock = out / "scheduler.lock"
    lock.mkdir()
    state = dict(status="running", active_chromosomes=[4, 6],
                 completed_chromosomes=[1, 2, 3, 5, 22], exact_allele_recovery=True)
    base.adapter.write_json(out / "durable_scheduler.json", state)
    try:
        scheduler.run_batch([
            ("chr4_alleles", [sys.executable, str(Path(__file__)), "--stage", "chromosome", "--output", str(out)]),
            ("chr6", [sys.executable, str(Path(base.__file__)), "--stage", "chromosome",
                      "--chromosome", "6", "--output", str(out), "--retry-incomplete"])], out, state)
    finally:
        lock.rmdir()
    scheduler.schedule(out, recover_preparation=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "smoke", "chromosome", "run"), required=True)
    parser.add_argument("--output", type=Path, default=Path("/validation"))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.output)
    elif args.stage == "smoke":
        config = base.configuration(args.output)
        prepare(args.output)
        chosen = next(info for info in config["cases"] if info["mask"] == "frederik")
        case(args.output, chosen)
        base.adapter.write_json(args.output / "allele_smoke_passed.json", dict(
            status="complete", case_id=chosen["case_id"],
            recovery_sha256=base.adapter.digest(args.output / "allele_recovery.json")))
    elif args.stage == "chromosome":
        chromosome(args.output)
    else:
        run(args.output)
