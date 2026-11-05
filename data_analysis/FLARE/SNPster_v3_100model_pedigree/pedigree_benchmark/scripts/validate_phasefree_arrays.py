#!/usr/bin/env python3
"""Frozen-model, seen-donor seven-source stress tests under two real array masks."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import subprocess

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as shared
from realized_ancestry import load_bp
from screen_dosage_models import load_dosages
from screen_phasefree_stack import features
from screen_v1_decoders import predict_decoder
from train_matched_global import select_training, truth_vector


PROJECT = Path(__file__).resolve().parents[2]
MODEL = PROJECT / (
    "results_v3_phasefree_stack_01/candidates/phasefree_stack_boost_clean/model.pkl")
MODEL_SHA = "95f3e3c9abaf76fb357d4ef51b2e29bf212a033002c02a04b825ec01b0231686"
SCREEN = Path("/results")
REFERENCE = Path("/imputer")
PREVIOUS = PROJECT / "results_v3_matched_global_01"
PEDIGREE = PROJECT / "results_v3_pedigree_pilot_10/pedigree"
MASKS = ("lonni", "frederik")
SOURCES = list(adapter.SOURCES)


def pure_controls(rows, eligible):
    selected = []
    for source in SOURCES:
        choices = [row for row in rows if row["sample_id"] in eligible
                   and ("EAS_SEA" if row["source_group"] in ("EAS", "SEA")
                        else row["source_group"]) == source]
        if not choices:
            raise ValueError(f"No eligible dense pure-validation donor for {source}")
        row = min(choices, key=lambda row: hashlib.sha256(
            f"{adapter.SEED}:array-pure:{row['sample_id']}".encode()).digest())
        selected.append(dict(
            sample_id=row["sample_id"], family_id=f"pure:{row['sample_id']}",
            ancestors=[row["sample_id"]], truth=[float(value == source) for value in SOURCES],
            track="pure", batch="pure"))
    return selected


def configuration(out):
    config = json.loads((out / "experiment.json").read_text())
    shared.verify_hashes(config["implementation_sha256"])
    if adapter.digest(MODEL) != config["model_sha256"] or config["model_sha256"] != MODEL_SHA:
        raise ValueError("Frozen model changed")
    return config


def verified(path, experiment_hash):
    if not path.exists():
        return False
    record = json.loads(path.read_text())
    if record["status"] != "complete" or record["experiment_sha256"] != experiment_hash:
        raise ValueError(f"Incompatible checkpoint: {path}")
    shared.verify_hashes(record["output_sha256"])
    shared.verify_hashes(record["input_sha256"])
    return True


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Preparation requires fresh output")
    out.mkdir(parents=True, exist_ok=True)
    if adapter.digest(MODEL) != MODEL_SHA:
        raise ValueError("Frozen candidate changed")
    cohort = json.loads((SCREEN / "cohort.json").read_text())
    old = json.loads((PREVIOUS / "pilot.json").read_text())
    dense = set(subprocess.check_output([
        "bcftools", "query", "-l",
        str(REFERENCE / "beagle_references" / shared.REFERENCE_PATTERN.format(chrom=22))
    ], text=True).splitlines())
    eligible_file = PEDIGREE / "inputs/eligible_heldout_founders.tsv"
    eligible = {row["sample"] for row in adapter.records(eligible_file)} & dense
    reference_ids = {row["sample_id"] for row in adapter.records(
        cohort["reference"]["labels_tsv"])}
    train_ids = {row["sample_id"] for row in adapter.records(cohort["pure_train"]["labels_tsv"])}
    if eligible & (reference_ids | train_ids):
        raise ValueError("Eligible stress-test donors overlap fitted identities")
    pure_batches = cohort["evaluation"]["pure"]["batches"]
    if len(pure_batches) != 1:
        raise ValueError("Expected the frozen single pure-validation batch")
    pure_batch = pure_batches[0]
    samples = pure_controls(adapter.records(pure_batch["samples_tsv"]), eligible)
    batches = {"pure": dict(vcf=pure_batch["vcf"])}
    inputs = [MODEL, SCREEN / "cohort.json", PREVIOUS / "pilot.json", eligible_file,
              Path(cohort["reference"]["labels_tsv"]), Path(cohort["pure_train"]["labels_tsv"]),
              Path(pure_batch["samples_tsv"]), Path(pure_batch["vcf"]), PREVIOUS / "beagle.jar"]
    prior_cases = {row["sample_id"] for row in old["samples"] if row["role"] == "evaluation"}
    candidates = []
    for batch in cohort["evaluation"]["pedigree"]["batches"]:
        name = batch["batch"]
        plan_path = PEDIGREE / f"plans/plans/{name}/pedigrees.json"
        plan = json.loads(plan_path.read_text())
        batches[name] = dict(vcf=batch["vcf"], plan=str(plan_path))
        inputs += [plan_path, Path(batch["samples_tsv"]), Path(batch["vcf"])]
        for row in adapter.records(batch["samples_tsv"]):
            donors = shared.ancestor_donors(plan, row["sample_id"])
            if donors <= eligible and row["sample_id"] not in prior_cases:
                candidates.append(dict(
                    sample_id=row["sample_id"], family_id=row["family_id"], batch=name,
                    ancestors=sorted(donors), truth=truth_vector(row).tolist(), track="balanced"))
    samples += select_training(candidates)
    ids = [row["sample_id"] for row in samples]
    if len(set(ids)) != 14 or not np.allclose(np.sum([row["truth"] for row in samples], axis=0), 2):
        raise ValueError("Expected exactly source-balanced seven pure plus seven mixed controls")
    donor_ids = sorted({donor for row in samples for donor in row["ancestors"]})
    excluded = set(old["excluded_reference_ids"])
    if not set(donor_ids) <= excluded or set(donor_ids) & (reference_ids | train_ids):
        raise ValueError("Stress-test donors are not excluded from every fitted reference")
    masks = {}
    for name in MASKS:
        path = Path(f"/data/flare2_test/source_model/real_user_validation/{name}/{name}.autosomes.vcf.gz")
        target = out / f"{name}.observed_keys.tsv"
        with target.open("w") as stream:
            subprocess.run(["bcftools", "query", "-i", "INFO/IMP!=1",
                            "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(path)],
                           check=True, stdout=stream)
        masks[name] = str(target)
        inputs += [path, target]
    used_batches = {row["batch"] for row in samples}
    batches = {name: value for name, value in batches.items() if name in used_batches}
    for chrom in range(1, 23):
        for name in batches:
            if name != "pure":
                inputs.append(PEDIGREE / f"chromosomes/pdb_{name}_chr{chrom}.bp")
    implementation = [Path(__file__), Path(adapter.__file__), Path(shared.__file__),
                      Path(__file__).with_name("realized_ancestry.py"),
                      Path(__file__).with_name("screen_dosage_models.py"),
                      Path(__file__).with_name("screen_phasefree_stack.py"),
                      Path(__file__).with_name("screen_fine_profiles.py"),
                      Path(__file__).with_name("screen_v1_decoders.py"),
                      Path(__file__).with_name("train_matched_global.py")]
    cases = [dict(case_id=f"{mask}__{sample['sample_id']}", mask=mask, **sample)
             for mask in MASKS for sample in samples]
    config = dict(
        status="prepared", model_sha256=MODEL_SHA, seed=adapter.SEED,
        model_frozen=True, new_model_training=False, independent_final_test=False,
        seen_donor_stress_test=True, production_changed=False, pgs_qc_changed=False,
        input_contract="Complete pre-QC imputed hard-GT counts; no DS substitution or missing padding.",
        selection="Hash-selected one pure per source; seven donor/family-disjoint balanced pairs selected by fixed MILP, not predictions.",
        n_unique_individuals=14, n_array_conditions=28, n_singleton_imputations=616,
        samples=samples, cases=cases, masks=masks, batches=batches, donor_ids=donor_ids,
        excluded_reference_ids=sorted(excluded), n_dense_reference_after_exclusion=len(dense - excluded),
        minor_2_5_support=0, minor_2_5_certification=False,
        minor_scope="This bounded stage tests pure and 50:50 transfer, not low-level ancestry sensitivity.",
        diagnostic_mean_tv_max=.10,
        source_reference_exclusion_policy_unchanged=True,
        input_sha256=shared.fingerprint_files(inputs),
        implementation_sha256=shared.fingerprint_files(implementation),
        reference_sha256={path: value for path, value in old["input_sha256"].items()
                          if path.startswith("/imputer/")},
    )
    pd.DataFrame([dict(
        case_id=row["case_id"], sample_id=row["sample_id"], track=row["track"],
        mask=row["mask"], family_id=row["family_id"], ancestors=",".join(row["ancestors"]),
        **dict(zip(("truth_" + source for source in SOURCES), row["truth"])))
        for row in cases]).to_csv(out / "cases.tsv", sep="\t", index=False)
    config["input_sha256"][str(out / "cases.tsv")] = adapter.digest(out / "cases.tsv")
    adapter.write_json(out / "experiment.json", config)
    print(f"Frozen design: {len(samples)} people, {len(cases)} array conditions, "
          f"{len(donor_ids)} seen donors; no minor-ancestry certification.", flush=True)


def write_truth_vcf(path, chrom, keys, names, haplotypes, unphased=False):
    if haplotypes.shape != (len(names), 2, len(keys)) or not np.isin(haplotypes, [0, 1]).all():
        raise ValueError("Invalid truth/sample/marker axes")
    with path.open("w") as stream:
        stream.write(f'##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                     + "\t".join(names) + "\n")
        for i, key in enumerate(keys):
            if unphased:
                calls = "\t".join("/".join(map(str, sorted(x[:, i]))) for x in haplotypes)
            else:
                calls = "\t".join(f"{x[0, i]}|{x[1, i]}" for x in haplotypes)
            stream.write(f"{chrom}\t{key[1]}\t.\t{key[2]}\t{key[3]}\t.\tPASS\t.\tGT\t{calls}\n")


def prepare_chromosome(out, chrom):
    config = configuration(out)
    experiment_hash = adapter.digest(out / "experiment.json")
    directory = out / f"chr{chrom}"
    if verified(directory / "prepared.json", experiment_hash):
        return
    if directory.exists():
        raise ValueError(f"Incomplete chromosome preparation requires inspection: {directory}")
    shared.verify_hashes(config["input_sha256"])
    cohort = json.loads((SCREEN / "cohort.json").read_text())
    with MODEL.open("rb") as stream:
        keys = [key for key in pickle.load(stream)["marker_keys"] if key[0] == str(chrom)]
    ref = REFERENCE / "beagle_references" / shared.REFERENCE_PATTERN.format(chrom=chrom)
    map_path = REFERENCE / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map"
    for path in (ref, Path(str(ref) + ".tbi"), map_path):
        if adapter.digest(path) != config["reference_sha256"][str(path)]:
            raise ValueError(f"Frozen imputation reference/map changed: {path}")
    directory.mkdir(exist_ok=True)
    span = cohort["chromosome_spans"][str(chrom)]
    observed = {}
    for name, path in config["masks"].items():
        observed[name] = {(row[0], int(row[1]), row[2], row[3])
                          for line in Path(path).read_text().splitlines()
                          if (row := line.split("\t"))[0] == str(chrom)
                          and span[0] <= int(row[1]) < span[1]}
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
            raise ValueError("Duplicate exact dense-reference marker")
        seen.add(key)
        values = dict(field.split("=") for field in fields[4:])
        found.append(key)
        donor_values.append([allele for sample in config["donor_ids"]
                             for allele in adapter.phased_gt(values[sample])])
    if not set(keys) <= seen or len({key[1] for key in found}) != len(found):
        raise ValueError("Missing frozen markers or ambiguous observed array positions")
    donors = np.asarray(donor_values, dtype=np.int8).T
    if not np.isin(donors, [0, 1]).all():
        raise ValueError("Dense founder panel has missing alleles")
    positions = np.asarray([key[1] for key in found])
    index = {key: i for i, key in enumerate(found)}
    model_indices = [index[key] for key in keys]
    names = [row["sample_id"] for row in config["samples"]]
    truth = []
    for sample in config["samples"]:
        batch = config["batches"][sample["batch"]]
        if sample["track"] == "pure":
            offset = 2 * config["donor_ids"].index(sample["sample_id"])
            x = donors[offset:offset + 2]
        else:
            plan = json.loads(Path(batch["plan"]).read_text())
            bp = load_bp(PEDIGREE / f"chromosomes/pdb_{sample['batch']}_chr{chrom}.bp", chrom)
            x = shared.pedigree_haplotypes(
                plan, bp, sample["sample_id"], config["donor_ids"], donors, positions)
        original = adapter.load_haplotypes(batch["vcf"], [sample["sample_id"]], chrom, keys)
        if not np.array_equal(original, x[:, model_indices]):
            raise ValueError("Dense reconstruction does not reproduce original canonical haplotypes")
        truth.append(x)
    truth = np.asarray(truth)
    truth_path = directory / "truth.vcf"
    write_truth_vcf(truth_path, chrom, found, names, truth)
    shared.checked(["bcftools", "view", "-Oz", "-o", directory / "truth.vcf.gz", truth_path],
                   directory, "compress_truth")
    shared.checked(["bcftools", "index", "-f", "-t", directory / "truth.vcf.gz"],
                   directory, "index_truth")
    truth_path.unlink()
    pd.DataFrame(found, columns=["chrom", "pos", "ref", "alt"]).to_csv(
        directory / "truth_markers.tsv", sep="\t", index=False)
    np.save(directory / "truth_haplotypes.npy", truth)
    original_ids = subprocess.check_output(["bcftools", "query", "-l", str(ref)], text=True).splitlines()
    keep = [sample for sample in original_ids if sample not in set(config["excluded_reference_ids"])]
    if len(keep) != config["n_dense_reference_after_exclusion"] or set(keep) & set(config["donor_ids"]):
        raise ValueError("Incorrect imputation-reference exclusion")
    keep_path, disjoint = directory / "reference.keep.txt", directory / "reference.disjoint.vcf.gz"
    keep_path.write_text("".join(sample + "\n" for sample in keep))
    if not disjoint.exists():
        shared.checked(["bcftools", "view", "-S", keep_path, "-Oz", "--threads", "2",
                        "-o", disjoint, ref], directory, "prepare_reference")
    actual = subprocess.check_output(["bcftools", "query", "-l", str(disjoint)], text=True).splitlines()
    if actual != keep:
        raise ValueError("Prepared imputation-reference identity axis mismatch")
    inputs = shared.fingerprint_files(
        [ref, map_path, PREVIOUS / "beagle.jar", directory / "truth.vcf.gz",
         directory / "truth_markers.tsv", directory / "truth_haplotypes.npy"])
    adapter.write_json(directory / "prepared.json", dict(
        status="complete", experiment_sha256=experiment_hash,
        reference_cache_sha256=adapter.digest(disjoint),
        input_sha256=inputs,
        output_sha256=shared.fingerprint_files([keep_path, directory / "truth.vcf.gz",
                                               directory / "truth.vcf.gz.tbi",
                                               directory / "truth_markers.tsv",
                                               directory / "truth_haplotypes.npy"])))


def chromosome(out, chrom, retry_incomplete=False):
    config = configuration(out)
    shared.verify_hashes(config["input_sha256"])
    experiment_hash = adapter.digest(out / "experiment.json")
    directory = out / f"chr{chrom}"
    if verified(directory / "completion.json", experiment_hash):
        print(f"chr{chrom}: verified completed checkpoint", flush=True)
        return
    prepare_chromosome(out, chrom)
    prepared = json.loads((directory / "prepared.json").read_text())
    inputs = prepared["input_sha256"]
    disjoint = directory / "reference.disjoint.vcf.gz"
    if adapter.digest(disjoint) != prepared["reference_cache_sha256"]:
        raise ValueError("Incomplete/changed prepared reference cache")
    cohort = json.loads((SCREEN / "cohort.json").read_text())
    span = cohort["chromosome_spans"][str(chrom)]
    observed = {
        name: {(row[0], int(row[1]), row[2], row[3])
               for line in Path(path).read_text().splitlines()
               if (row := line.split("\t"))[0] == str(chrom)
               and span[0] <= int(row[1]) < span[1]}
        for name, path in config["masks"].items()
    }
    found = [(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
             for row in adapter.records(directory / "truth_markers.tsv")]
    seen = set(found)
    truth = np.load(directory / "truth_haplotypes.npy")
    names = [row["sample_id"] for row in config["samples"]]
    with MODEL.open("rb") as stream:
        keys = [key for key in pickle.load(stream)["marker_keys"] if key[0] == str(chrom)]
    index = {key: i for i, key in enumerate(found)}
    model_indices = [index[key] for key in keys]
    keep = (directory / "reference.keep.txt").read_text().splitlines()
    map_path = REFERENCE / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map"
    recovered = []
    for case in config["cases"]:
        work = directory / case["case_id"]
        checkpoint = work / "completion.json"
        if verified(checkpoint, experiment_hash):
            recovered.append(np.load(work / "dosages.npy"))
            continue
        if work.exists():
            if not retry_incomplete:
                raise ValueError(f"Incomplete case requires explicit --resume: {work}")
            attempt = 1
            archive = work.with_name(f"{work.name}.incomplete_{attempt}")
            while archive.exists():
                attempt += 1
                archive = work.with_name(f"{work.name}.incomplete_{attempt}")
            work.rename(archive)
            print(f"Preserved incomplete attempt: {archive}", flush=True)
        work.mkdir()
        sample_index = names.index(case["sample_id"])
        array_indices = [i for i, key in enumerate(found) if key in observed[case["mask"]]]
        if not array_indices:
            raise ValueError("No reconstructable observed sites for an array/chromosome")
        array_keys = [found[i] for i in array_indices]
        target = work / "array.unphased.vcf"
        write_truth_vcf(target, chrom, array_keys, [case["sample_id"]],
                        truth[sample_index:sample_index + 1, :, array_indices], unphased=True)
        prefix = work / "beagle"
        seconds = shared.checked([
            "java", "-Xmx16g", "-jar", PREVIOUS / "beagle.jar", f"gt={target}",
            f"ref={disjoint}", f"map={map_path}", f"out={prefix}", "nthreads=6",
            f"seed={config['seed']}", "gp=true"], work, "impute")
        imputed = Path(str(prefix) + ".vcf.gz")
        shared.checked(["bcftools", "index", "-t", imputed], work, "index_imputed")
        retained = work / "retained.vcf.gz"
        shared.checked(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz",
                        "-o", retained, imputed, directory / "truth.vcf.gz"], work, "retain")
        shared.checked(["bcftools", "index", "-t", retained], work, "index_retained")
        array_regions = work / "array.positions.tsv"
        array_regions.write_text("".join(f"{chrom}\t{key[1]}\n" for key in array_keys))
        array_doses, _ = load_dosages(retained, [case["sample_id"]], array_keys, array_regions)
        if not np.array_equal(array_doses[0], truth[sample_index][:, array_indices].sum(axis=0)):
            raise ValueError("Imputation changed a supplied observed genotype")
        model_regions = work / "model.positions.tsv"
        model_regions.write_text("".join(f"{chrom}\t{key[1]}\n" for key in keys))
        values, _ = load_dosages(retained, [case["sample_id"]], keys, model_regions)
        np.save(work / "dosages.npy", values[0])
        stats = dict(
            status="complete", experiment_sha256=experiment_hash, seconds=seconds,
            observed_genotypes_preserved=True, n_observed=len(array_indices),
            n_unavailable_observed=len(observed[case["mask"]] - seen),
            n_model_markers=len(keys), n_genotype_errors=int(np.count_nonzero(
                values[0] != truth[sample_index][:, model_indices].sum(axis=0))),
            input_sha256=shared.fingerprint_files(
                [directory / "prepared.json", PREVIOUS / "beagle.jar", target]),
            output_sha256=shared.fingerprint_files([retained, work / "dosages.npy"]))
        adapter.write_json(checkpoint, stats)
        recovered.append(values[0])
        imputed.unlink()
        Path(str(imputed) + ".tbi").unlink()
        print(f"chr{chrom} {case['case_id']}: verified in {seconds:.1f}s", flush=True)
    np.savez_compressed(directory / "dosages.npz",
                        original=truth[:, :, model_indices].sum(axis=1),
                        array_imputed=np.asarray(recovered))
    adapter.write_json(directory / "marker_keys.json", keys)
    outputs = [directory / "dosages.npz", directory / "marker_keys.json",
               *[directory / row["case_id"] / "completion.json" for row in config["cases"]]]
    adapter.write_json(directory / "completion.json", dict(
        status="complete", experiment_sha256=experiment_hash,
        n_cases=len(recovered), n_reference_after=len(keep),
        original_haplotypes_reproduced=True, input_sha256=inputs,
        output_sha256=shared.fingerprint_files(outputs)))
    disjoint.unlink()
    print(f"chr{chrom}: all {len(recovered)} cases complete", flush=True)


def summarize(out):
    config = configuration(out)
    shared.verify_hashes(config["input_sha256"])
    with MODEL.open("rb") as stream:
        bundle = pickle.load(stream)
    keys = bundle["marker_keys"]
    offsets = {key: i for i, key in enumerate(keys)}
    original = np.empty((14, len(keys)), dtype=np.int8)
    imputed = np.empty((28, len(keys)), dtype=np.int8)
    seen = set()
    for chrom in range(1, 23):
        directory = out / f"chr{chrom}"
        if not verified(directory / "completion.json", adapter.digest(out / "experiment.json")):
            raise ValueError(f"Missing chromosome {chrom}")
        chrom_keys = [tuple(key) for key in json.loads((directory / "marker_keys.json").read_text())]
        if seen & set(chrom_keys):
            raise ValueError("Duplicate chromosome marker coverage")
        seen.update(chrom_keys)
        indices = [offsets[key] for key in chrom_keys]
        with np.load(directory / "dosages.npz") as data:
            original[:, indices] = data["original"]
            imputed[:, indices] = data["array_imputed"]
    if seen != set(keys):
        raise ValueError("Incomplete frozen marker axis")
    with threadpool_limits(limits=6):
        clean_predictions = predict_decoder(bundle["decoder"], features(bundle["feature_models"], original))
        array_predictions = predict_decoder(bundle["decoder"], features(bundle["feature_models"], imputed))
    records = []
    by_id = {row["sample_id"]: i for i, row in enumerate(config["samples"])}
    for i, case in enumerate(config["cases"]):
        truth = np.asarray(case["truth"])
        for stage, values in (("clean", clean_predictions[by_id[case["sample_id"]]]),
                              ("array_imputed", array_predictions[i])):
            records.append(dict(
                case_id=case["case_id"], sample_id=case["sample_id"], mask=case["mask"],
                track=case["track"], stage=stage,
                tv=float(np.abs(values - truth).sum() / 2),
                off_target=float(values[truth == 0].sum()),
                **dict(zip(("truth_" + source for source in SOURCES), truth)),
                **dict(zip(SOURCES, values))))
    frame = pd.DataFrame(records)
    frame.to_csv(out / "predictions.tsv", sep="\t", index=False)
    groups = [dict(mask=mask, track=track, stage=stage, n=len(group),
                   mean_tv=float(group.tv.mean()), p95_tv=float(group.tv.quantile(.95)),
                   max_tv=float(group.tv.max()), mean_off_target=float(group.off_target.mean()))
              for (mask, track, stage), group in frame.groupby(["mask", "track", "stage"])]
    pd.DataFrame(groups).to_csv(out / "metrics.tsv", sep="\t", index=False)
    adapter.write_json(out / "completion.json", dict(
        status="complete_seen_donor_stress_test", model_sha256=adapter.digest(MODEL),
        experiment_sha256=adapter.digest(out / "experiment.json"),
        n_unique_individuals=14, n_array_conditions=28, metrics=groups,
        diagnostic_mean_tv_max=config["diagnostic_mean_tv_max"],
        array_mean_tv_diagnostics=[
            dict(mask=row["mask"], track=row["track"],
                 passed=row["mean_tv"] <= config["diagnostic_mean_tv_max"])
            for row in groups if row["stage"] == "array_imputed"],
        independent_final_test=False, minor_2_5_certification=False,
        production_promotion=False, pgs_qc_changed=False,
        output_sha256=shared.fingerprint_files([out / "predictions.tsv", out / "metrics.tsv"])))
    print(frame.groupby(["mask", "track", "stage"]).tv.mean().to_string(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "chromosome", "summarize"), required=True)
    parser.add_argument("--output", type=Path, default=Path("/validation"))
    parser.add_argument("--chromosome", type=int, choices=range(1, 23))
    parser.add_argument("--retry-incomplete", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.output)
    elif args.stage == "chromosome":
        if args.chromosome is None:
            parser.error("--chromosome is required for chromosome stage")
        chromosome(args.output, args.chromosome, args.retry_incomplete)
    else:
        summarize(args.output)
