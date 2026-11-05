#!/usr/bin/env python3
"""Independently reread frozen array-test genotypes and reproduce final scores."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import pickle
import subprocess

import numpy as np
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from screen_phasefree_stack import features
from screen_v1_decoders import predict_decoder


PROJECT = Path(__file__).resolve().parents[2]
EXPECTED_MODEL_SHA = "95f3e3c9abaf76fb357d4ef51b2e29bf212a033002c02a04b825ec01b0231686"


def dosage(value):
    alleles = value.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in ("0", "1") for allele in alleles):
        raise ValueError(f"Incomplete/nonbiallelic verification genotype: {value}")
    return sum(map(int, alleles))


def one_sample_vcf(path, sample):
    header, seen = False, set()
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()) as stream:
        for line in stream:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if fields[9:] != [sample]:
                    raise ValueError(f"Verification sample-axis mismatch: {path}")
                header = True
                continue
            if not header or len(fields) != 10:
                raise ValueError(f"Invalid single-sample VCF: {path}")
            values = dict(zip(fields[8].split(":"), fields[9].split(":"), strict=True))
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key in seen:
                raise ValueError(f"Duplicate verification VCF marker: {path}")
            seen.add(key)
            yield key, dosage(values["GT"])
    if not header:
        raise ValueError(f"Missing verification VCF header: {path}")


def verify(out):
    config_path, completion_path = out / "experiment.json", out / "completion.json"
    config, completion = json.loads(config_path.read_text()), json.loads(completion_path.read_text())
    if (completion["status"] != "complete_seen_donor_stress_test"
            or completion["experiment_sha256"] != adapter.digest(config_path)):
        raise ValueError("Stress-test completion/configuration mismatch")
    checked = {}

    def hashes(values):
        for path, expected in values.items():
            if path in checked:
                if checked[path] != expected:
                    raise ValueError(f"Inconsistent expected digest: {path}")
                continue
            if adapter.digest(path) != expected:
                raise ValueError(f"Verification digest mismatch: {path}")
            checked[path] = expected

    hashes(config["input_sha256"])
    hashes(config["implementation_sha256"])
    hashes(completion["output_sha256"])
    model_path = PROJECT / (
        "results_v3_phasefree_stack_01/candidates/phasefree_stack_boost_clean/model.pkl")
    if adapter.digest(model_path) != EXPECTED_MODEL_SHA or config["model_sha256"] != EXPECTED_MODEL_SHA:
        raise ValueError("Frozen model digest mismatch")
    with model_path.open("rb") as stream:
        model = pickle.load(stream)
    keys = model["marker_keys"]
    offsets = {key: i for i, key in enumerate(keys)}
    if len(offsets) != 26066 or len(config["samples"]) != 14 or len(config["cases"]) != 28:
        raise ValueError("Frozen verification sample/marker counts changed")
    names = [row["sample_id"] for row in config["samples"]]
    if len(set(names)) != 14 or len({row["case_id"] for row in config["cases"]}) != 28:
        raise ValueError("Duplicate verification individuals/conditions")
    cohort = json.loads(Path("/results/cohort.json").read_text())
    truth_rows = {
        row["sample_id"]: row
        for track in ("pure", "pedigree") for batch in cohort["evaluation"][track]["batches"]
        for row in adapter.records(batch["samples_tsv"])
    }
    fitted = {
        row["sample_id"] for group in ("reference", "pure_train")
        for row in adapter.records(cohort[group]["labels_tsv"])
    }
    if set(config["donor_ids"]) & fitted:
        raise ValueError("Verification found fitted/evaluation donor overlap")
    for sample in config["samples"]:
        row = truth_rows[sample["sample_id"]]
        expected = [float(row["realized_EAS"]) + float(row["realized_SEA"])
                    if source == "EAS_SEA" else float(row["realized_" + source])
                    for source in adapter.SOURCES]
        np.testing.assert_allclose(expected, sample["truth"], atol=1e-12, rtol=0)
    clean, imputed = np.empty((14, len(keys)), dtype=np.int8), np.empty((28, len(keys)), dtype=np.int8)
    all_keys, observed_count = set(), 0
    for chrom in range(1, 23):
        directory = out / f"chr{chrom}"
        record = json.loads((directory / "completion.json").read_text())
        if record["status"] != "complete" or record["experiment_sha256"] != adapter.digest(config_path):
            raise ValueError(f"Chromosome {chrom} completion mismatch")
        hashes(record["input_sha256"])
        hashes(record["output_sha256"])
        prepared = json.loads((directory / "prepared.json").read_text())
        if prepared["status"] != "complete" or prepared["experiment_sha256"] != adapter.digest(config_path):
            raise ValueError("Chromosome preparation verification mismatch")
        hashes(prepared["input_sha256"])
        hashes(prepared["output_sha256"])
        keep = (directory / "reference.keep.txt").read_text().splitlines()
        if (len(set(keep)) != len(keep) or len(keep) != 3731
                or set(keep) & set(config["excluded_reference_ids"])):
            raise ValueError("Verification reference exclusion mismatch")
        chrom_keys = [tuple(key) for key in json.loads((directory / "marker_keys.json").read_text())]
        if len(set(chrom_keys)) != len(chrom_keys) or all_keys & set(chrom_keys):
            raise ValueError("Duplicate chromosome verification marker")
        all_keys.update(chrom_keys)
        positions = {key: i for i, key in enumerate(chrom_keys)}
        indices = [offsets[key] for key in chrom_keys]
        original = np.empty((14, len(indices)), dtype=np.int8)
        for batch, info in config["batches"].items():
            people = [row["sample_id"] for row in config["samples"] if row["batch"] == batch]
            text = subprocess.check_output([
                "bcftools", "query", "-R", str(directory / "positions.tsv"),
                "-s", ",".join(people), "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n",
                info["vcf"]], text=True)
            found = set()
            for line in text.splitlines():
                fields = line.split("\t")
                key = (fields[0], int(fields[1]), fields[2], fields[3])
                if key not in positions:
                    continue
                if key in found:
                    raise ValueError("Duplicate original verification genotype")
                found.add(key)
                calls = dict(value.split("=") for value in fields[4:])
                if set(calls) != set(people):
                    raise ValueError("Original verification sample axis mismatch")
                for person in people:
                    original[names.index(person), positions[key]] = dosage(calls[person])
            if found != set(chrom_keys):
                raise ValueError("Original verification marker axis incomplete")
        values = []
        for case in config["cases"]:
            work = directory / case["case_id"]
            checkpoint = json.loads((work / "completion.json").read_text())
            if (checkpoint["status"] != "complete"
                    or checkpoint["experiment_sha256"] != adapter.digest(config_path)):
                raise ValueError("Case verification checkpoint mismatch")
            hashes(checkpoint["input_sha256"])
            hashes(checkpoint["output_sha256"])
            observed = dict(one_sample_vcf(work / "array.unphased.vcf", case["sample_id"]))
            found_observed, found_model = set(), {}
            for key, value in one_sample_vcf(work / "retained.vcf.gz", case["sample_id"]):
                if key in observed:
                    if key in found_observed or value != observed[key]:
                        raise ValueError("Supplied observed genotype changed")
                    found_observed.add(key)
                if key in positions:
                    if key in found_model:
                        raise ValueError("Duplicate imputed verification marker")
                    found_model[key] = value
            if found_observed != set(observed) or set(found_model) != set(chrom_keys):
                raise ValueError("Incomplete observed/model verification coverage")
            observed_count += len(observed)
            row = np.asarray([found_model[key] for key in chrom_keys], dtype=np.int8)
            np.testing.assert_array_equal(row, np.load(work / "dosages.npy"))
            if np.count_nonzero(row != original[names.index(case["sample_id"])]) != checkpoint["n_genotype_errors"]:
                raise ValueError("Incorrect saved genotype-error count")
            values.append(row)
        values = np.asarray(values)
        with np.load(directory / "dosages.npz") as saved:
            np.testing.assert_array_equal(original, saved["original"])
            np.testing.assert_array_equal(values, saved["array_imputed"])
        clean[:, indices], imputed[:, indices] = original, values
        print(f"chr{chrom}: independent original/imputed and observed-genotype checks passed", flush=True)
    if all_keys != set(keys):
        raise ValueError("Whole-genome verification marker axis incomplete")
    with threadpool_limits(limits=6):
        clean_predictions = predict_decoder(model["decoder"], features(model["feature_models"], clean))
        array_predictions = predict_decoder(model["decoder"], features(model["feature_models"], imputed))
    with (out / "predictions.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    indexed = {(row["case_id"], row["stage"]): row for row in rows}
    if len(rows) != 56 or len(indexed) != 56:
        raise ValueError("Unexpected verification prediction row count")
    scores, off_targets = {}, {}
    for i, case in enumerate(config["cases"]):
        truth = np.asarray(case["truth"])
        for stage, expected in (
            ("clean", clean_predictions[names.index(case["sample_id"])]),
            ("array_imputed", array_predictions[i]),
        ):
            row = indexed[case["case_id"], stage]
            if (row["sample_id"] != case["sample_id"] or row["mask"] != case["mask"]
                    or row["track"] != case["track"]):
                raise ValueError("Verification prediction metadata mismatch")
            np.testing.assert_allclose([float(row[source]) for source in adapter.SOURCES],
                                       expected, atol=1e-12, rtol=0)
            tv = float(np.abs(expected - truth).sum() / 2)
            off_target = float(expected[truth == 0].sum())
            if (abs(float(row["tv"]) - tv) > 1e-12
                    or abs(float(row["off_target"]) - off_target) > 1e-12):
                raise ValueError("Verification TV/off-target mismatch")
            scores.setdefault((case["mask"], case["track"], stage), []).append(tv)
            off_targets.setdefault((case["mask"], case["track"], stage), []).append(off_target)
    with (out / "metrics.tsv").open() as stream:
        metrics = list(csv.DictReader(stream, delimiter="\t"))
    if len(metrics) != len(scores):
        raise ValueError("Verification metric group count mismatch")
    for row in metrics:
        values = scores[row["mask"], row["track"], row["stage"]]
        if int(row["n"]) != len(values):
            raise ValueError("Verification metric support mismatch")
        np.testing.assert_allclose(
            [float(row["mean_tv"]), float(row["p95_tv"]), float(row["max_tv"]),
             float(row["mean_off_target"])],
            [np.mean(values), np.quantile(values, .95), max(values),
             np.mean(off_targets[row["mask"], row["track"], row["stage"]])],
            atol=1e-12, rtol=0)
    for row in completion["array_mean_tv_diagnostics"]:
        expected = np.mean(scores[row["mask"], row["track"], "array_imputed"]) <= .10
        if row["passed"] != expected:
            raise ValueError("Verification array-diagnostic mismatch")
    if (adapter.digest(model_path) != EXPECTED_MODEL_SHA
            or completion["model_sha256"] != EXPECTED_MODEL_SHA
            or completion["independent_final_test"] or completion["minor_2_5_certification"]
            or completion["production_promotion"]):
        raise ValueError("Frozen-model/interpretation constraints violated")
    adapter.write_json(out / "verification.json", dict(
        status="complete_verified_seen_donor_stress_test",
        experiment_sha256=adapter.digest(config_path),
        completion_sha256=adapter.digest(completion_path), verifier_sha256=adapter.digest(__file__),
        model_sha256=EXPECTED_MODEL_SHA, n_unique_hashes=len(checked),
        n_independent_original_genotype_cells=int(clean.size),
        n_independent_imputed_genotype_cells=int(imputed.size),
        n_observed_genotypes_verified=observed_count,
        prediction_absolute_tolerance=1e-12,
        independent_final_test=False, minor_2_5_certification=False, production_promotion=False))
    print("Independent frozen-model array stress-test verification complete.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("/validation"))
    verify(parser.parse_args().output)
