#!/usr/bin/env python3
"""Eight-source native ancestry models on verified, unchanged production outputs."""
import argparse
from collections import defaultdict
import copy
import gzip
import json
from pathlib import Path
import pickle
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from audit_production_ancestry import SOURCES, digest, write_json
import model_family_gnomix as native
from production_matched_preprocessing import (
    PROJECT, ROOT, compress, nucleotide_calls, positions, read, read_call_matrix, reference, write_vcf,
)
from sgdp_genetic_qc import query_keys, run


def validate_vector(value):
    value = np.asarray(value, dtype=float)
    if (value.shape != (8,) or not np.isfinite(value).all() or (value < 0).any()
            or abs(float(value.sum()) - 1) > 1e-10):
        raise ValueError("Invalid eight-source probability vector")
    return value


def constrained_frequencies(frequencies, dosage):
    frequencies, dosage = np.asarray(frequencies), np.asarray(dosage)
    if (frequencies.ndim != 2 or frequencies.shape[1] != 8
            or dosage.shape != (len(frequencies),) or len(dosage) < 100
            or not np.isfinite(frequencies).all() or not np.isin(dosage, [0, 1, 2]).all()):
        raise ValueError("Require at least100 called, aligned genotype/frequency SNPs")
    a, b = 2 * frequencies, dosage
    gram, cross = a.T @ a / len(b), a.T @ b / len(b)
    result = minimize(lambda x: float(x @ gram @ x - 2 * cross @ x), np.ones(8) / 8,
                      jac=lambda x: 2 * (gram @ x - cross), method="SLSQP",
                      bounds=[(0, 1)] * 8,
                      constraints={"type": "eq", "fun": lambda x: x.sum() - 1,
                                   "jac": lambda x: np.ones(8)},
                      options={"ftol": 1e-12, "maxiter": 1000})
    if not result.success:
        raise ValueError(f"Constrained frequency estimation failed: {result.message}")
    return validate_vector(result.x)


def panels(cohort_path, chromosome, masks, imputed_folders, merged_files, output, truth_folder=None):
    output.mkdir()
    cohort, protocol = read(cohort_path), read(cohort_path)["protocol"]
    references = cohort["reference"]
    ids = [row["sample_id"] for row in references]
    (output / "reference.ids").write_text("\n".join(ids) + "\n")
    policy = protocol["ancestry_panel"]
    run(["bcftools", "view", "-S", output / "reference.ids", "-r", str(chromosome),
         "-m2", "-M2", "-v", "snps", "-q", f"{policy['minimum_maf']}:minor",
         "-c", f"{policy['minimum_mac']}:minor", "-Oz", "-o", output / "filtered.vcf.gz",
         reference(chromosome)], output / "reference_quality")
    subprocess.run(["bcftools", "index", str(output / "filtered.vcf.gz")], check=True)
    keys = query_keys(output / "filtered.vcf.gz", chromosome)
    span = cohort["chromosome_spans"][str(chromosome)]
    map_path = ROOT / "plink.GRCh38.map" / f"plink.chr{chromosome}.GRCh38.map"
    map_points = native.map_points(map_path, chromosome)
    common = set(keys)
    for folder in imputed_folders:
        common &= set(query_keys(folder / f"imputed_chr{chromosome}.risk.vcf.gz", chromosome))
    chosen, previous = [], -policy["minimum_spacing_bp"]
    position_counts = defaultdict(int)
    for key in keys:
        position_counts[key[1]] += 1
    for key in keys:
        if (key in common and position_counts[key[1]] == 1
                and span[0] <= key[1] < span[1] and map_points[0, 0] <= key[1] <= map_points[-1, 0]
                and key[1] - previous >= policy["minimum_spacing_bp"]):
            chosen.append(key)
            previous = key[1]
    if len(chosen) < 100:
        raise ValueError("Insufficient reference-quality, shared imputed ancestry markers")
    region = positions(output / "panel_positions.tsv", chosen)
    run(["bcftools", "view", "-R", region, "-Oz", "-o", output / "reference.vcf.gz",
         output / "filtered.vcf.gz"], output / "thin_reference")
    subprocess.run(["bcftools", "index", str(output / "reference.vcf.gz")], check=True)
    ref = read_call_matrix(output / "reference.vcf.gz", chosen, ids, output, phased=True)
    reference_haps = ref.reshape(2 * len(ids), len(chosen))
    labels = np.asarray([SOURCES.index(row["source"]) for row in references])
    frequencies = np.stack([ref[labels == index].mean(axis=(0, 1)) for index in range(8)], axis=1)
    np.savez_compressed(output / "reference.npz", haplotypes=reference_haps, labels=labels,
                        frequencies=frequencies)
    (output / "full.panel").write_text("".join(f"{row['sample_id']}\t{row['source']}\n" for row in references))
    by_fine = defaultdict(list)
    for row in references:
        by_fine[row["fine"]].append(row)
    capped = [row for fine in sorted(by_fine)
              for row in sorted(by_fine[fine], key=lambda row: row["sample_id"])[:policy["fine_cap"]]]
    (output / "finecap.ids").write_text("".join(row["sample_id"] + "\n" for row in capped))
    (output / "finecap.panel").write_text("".join(f"{row['sample_id']}\t{row['source']}\n" for row in capped))
    run(["bcftools", "view", "-S", output / "finecap.ids", "-Oz", "-o", output / "finecap.vcf.gz",
         output / "reference.vcf.gz"], output / "cap_reference")
    subprocess.run(["bcftools", "index", str(output / "finecap.vcf.gz")], check=True)
    evaluation = [row for row in cohort["people"] if row["role"] == "evaluation"]
    training = [row for row in cohort["people"] if row["role"] == "training"]
    training_truth, training_ancestry = None, None
    if cohort.get("truth_design") == "synthetic_reference_haplotype_development":
        from dense_split_preprocessing import synthesise
        if truth_folder is None or not (truth_folder / "cases.json").is_file():
            raise ValueError("Synthetic training requires its declared verified case inputs")
        training_founders = [row["sample_id"] for row in cohort["training_founders"]]
        founder_calls = read_call_matrix(reference(chromosome), chosen, training_founders, output, phased=True)
        training_truth, training_ancestry = synthesise(
            cohort, chromosome, chosen, training_founders, founder_calls, training)
    panel_records = {}
    for mask, folder, merged in zip(masks, imputed_folders, merged_files):
        imputed = folder / f"imputed_chr{chromosome}.risk.vcf.gz"
        for role, people in (("evaluation", evaluation), ("training", training)):
            production_ids = [f"IMPID_{row['production_ids'][mask]}" for row in people]
            calls = read_call_matrix(imputed, chosen, production_ids, output, phased=True)
            path = write_vcf(output / f"{mask}.{role}.vcf", chosen,
                             [row["sample_id"] for row in people], calls, phased=True)
            arrays = dict(haplotypes=calls.reshape(2 * len(people), len(chosen)))
            if role == "training" and training_truth is not None:
                from dense_split_preprocessing import aligned_training_labels
                aligned, informative = aligned_training_labels(
                    training_truth, calls, training_ancestry, [key[1] for key in chosen])
                arrays["labels"] = aligned.reshape(2 * len(people), len(chosen))
                write_json(output / f"{mask}.training_alignment.json", dict(
                    policy="Training-only homolog alignment at concordant heterozygotes; nearest informative physical position through noninformative sites; no genotype changes. Modal SNP ancestry per upstream Gnomix window.",
                    people=[row["sample_id"] for row in training], informative_markers=informative,
                    evaluation_truth_used=False))
            np.savez_compressed(output / f"{mask}.{role}.npz", **arrays)
        raw_positions = {int(value) for value in subprocess.check_output(
            ["bcftools", "query", "-f", "%POS\\n", str(merged)], text=True).splitlines()}
        observed_keys = [key for key in keys if key[1] in raw_positions and position_counts[key[1]] == 1]
        observed_region = positions(output / f"{mask}.observed_positions.tsv", observed_keys)
        ref_observed = read_call_matrix(output / "filtered.vcf.gz", observed_keys, ids, output, phased=True)
        observed_freq = np.stack([ref_observed[labels == index].mean(axis=(0, 1))
                                  for index in range(8)], axis=1)
        observed_dosage = np.full((len(evaluation), len(observed_keys)), -1, dtype=np.int8)
        for person_index, person in enumerate(evaluation):
            actual = nucleotide_calls(merged, observed_region, f"IMPID_{person['production_ids'][mask]}")
            for index, key in enumerate(observed_keys):
                genotype = actual.get(key[1])
                if genotype is not None:
                    if any(base not in (key[2], key[3]) for base in genotype):
                        raise ValueError("Observed phase-independent alleles differ from exact reference keys")
                    observed_dosage[person_index, index] = genotype.count(key[3])
        np.savez_compressed(output / f"{mask}.observed.npz", frequencies=observed_freq,
                            dosage=observed_dosage)
        panel_records[mask] = dict(observed_reference_quality_markers=len(observed_keys))
    write_json(output / "panel.json", dict(
        chromosome=chromosome, masks=masks, keys=chosen, span=span, map=str(map_path),
        evaluation=[row["sample_id"] for row in evaluation],
        training=[row["sample_id"] for row in training],
        training_sources=[row.get("source") for row in training],
        reference_people=len(ids), finecap_people=len(capped), marker_count=len(chosen),
        stage_native_observed_panels=panel_records, source_order=list(SOURCES),
        raw_control_uses_only_actual_observed_calls=True))
    (output / "filtered.vcf.gz").unlink()
    Path(str(output / "filtered.vcf.gz") + ".csi").unlink()


def physical_flare(path, source_order, keys, span, expected_ids):
    wanted = {key: index for index, key in enumerate(keys)}
    positions_array = np.asarray([key[1] for key in keys], dtype=float)
    weights = np.diff(np.r_[span[0], (positions_array[:-1] + positions_array[1:]) / 2, span[1]])
    ids = None
    values = np.zeros((len(expected_ids), 8))
    seen = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                ids = line.rstrip().split("\t")[9:]
                if ids != expected_ids:
                    raise ValueError("FLARE changed expected sample identity/order")
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key not in wanted or key in seen or ids is None:
                raise ValueError("FLARE changed exact marker coverage")
            seen.add(key)
            names = fields[8].split(":")
            first, second = names.index("ANP1"), names.index("ANP2")
            for index, call in enumerate(fields[9:]):
                parts = call.split(":")
                a, b = [np.asarray(list(map(float, parts[column].split(",")))) for column in (first, second)]
                if (len(a) != 8 or len(b) != 8 or not np.isfinite(a).all() or not np.isfinite(b).all()
                        or (a < 0).any() or (b < 0).any()
                        or abs(a.sum() - 1) > .005 or abs(b.sum() - 1) > .005):
                    raise ValueError("Invalid FLARE posterior tensor")
                # FLARE writes rounded local probabilities; restore their probability masses explicitly.
                posterior = .5 * (a / a.sum() + b / b.sum())
                values[index] += posterior[[source_order.index(source) for source in SOURCES]] * weights[wanted[key]]
    if seen != set(keys):
        raise ValueError("Incomplete FLARE local posterior coverage")
    values /= span[1] - span[0]
    return np.stack([validate_vector(row) for row in values])


def fit_gnomix(folder, metadata, policy, mask, output):
    gnomix_class, windowed_class = native.load_upstream()
    keys = [tuple(key) for key in metadata["keys"]]
    configuration = dict(window_size_cM=policy["window_size_cM"], architecture="fast", smooth_size=1)
    layout = native.window_layout(keys, native.map_points(metadata["map"], metadata["chromosome"]),
                                  metadata["span"], configuration)
    data = np.load(folder / "reference.npz")
    x_path, y_path = output / f"{mask}.reference_X.npy", output / f"{mask}.reference_y.npy"
    np.save(x_path, data["haplotypes"])
    labels = np.repeat(data["labels"], 2)
    np.save(y_path, np.repeat(labels[:, None], len(keys), axis=1).astype(np.int8))
    model = gnomix_class(C=len(keys), M=layout["M"], A=8, S=1, mode="fast",
                         snp_pos=layout["positions"], snp_ref=np.asarray([key[2] for key in keys]),
                         snp_alt=np.asarray([key[3] for key in keys]), population_order=np.asarray(SOURCES),
                         missing_encoding=2, n_jobs=1, calibrate=False, context_ratio=.5,
                         mode_filter=False, seed=20261022, verbose=False)
    model.base.base_multithread = False
    for estimator in model.base.models:
        estimator.set_params(random_state=20261022)
    dataset = windowed_class([str(x_path)], [str(y_path)], len(keys), layout["M"])
    model.base.train_windowed(dataset, verbose=False)
    training_data = np.load(folder / f"{mask}.training.npz")
    training = training_data["haplotypes"]
    if "labels" in training_data:
        marker_labels = training_data["labels"]
        if marker_labels.shape != training.shape or not np.isin(marker_labels, range(8)).all():
            raise ValueError("Invalid native-eight-source local training ancestry")
        smooth_labels = native.window_labels(marker_labels, layout)
    else:
        training_labels = np.repeat([SOURCES.index(source) for source in metadata["training_sources"]], 2)
        smooth_labels = np.repeat(training_labels[:, None], layout["W"], axis=1)
    model.smooth.train(model.base.predict_proba(training), smooth_labels)
    with (output / f"{mask}.model.pkl").open("wb") as handle:
        pickle.dump(dict(model=model, layout=layout, source_order=SOURCES), handle)
    with (output / f"{mask}.model.pkl").open("rb") as handle:
        restored = pickle.load(handle)["model"]
    evaluation = np.load(folder / f"{mask}.evaluation.npz")["haplotypes"]
    before = model.smooth.predict_proba(model.base.predict_proba(evaluation))
    after = restored.smooth.predict_proba(restored.base.predict_proba(evaluation))
    np.testing.assert_allclose(before, after, rtol=0, atol=1e-10)
    if list(map(int, model.smooth.model.classes_)) != list(range(8)):
        raise ValueError("Gnomix changed eight-source class axis")
    if (after.shape != (2 * len(metadata["evaluation"]), layout["W"], 8)
            or not np.isfinite(after).all() or (after < 0).any()
            or not np.allclose(after.sum(axis=2), 1, rtol=0, atol=1e-10)):
        raise ValueError("Gnomix returned malformed local posterior probabilities")
    paired = after.reshape(-1, 2, layout["W"], 8).mean(axis=1)
    values = np.einsum("nwa,w->na", paired, layout["physical_weights"]) / (metadata["span"][1] - metadata["span"][0])
    np.savez_compressed(output / f"{mask}.local_posteriors.npz", probabilities=after,
                        weights=layout["physical_weights"])
    x_path.unlink()
    y_path.unlink()
    return np.stack([validate_vector(value) for value in values])


def evaluate(cohort_path, folder, candidate, output):
    output.mkdir()
    cohort, metadata = read(cohort_path), read(folder / "panel.json")
    policy = cohort["protocol"]["candidates"][candidate]
    keys = [tuple(key) for key in metadata["keys"]]
    predictions = []
    for mask in metadata["masks"]:
        start = time.monotonic()
        if policy["family"] == "flare":
            jar = PROJECT / "results_v3_global_precision_01/precision.jar"
            prefix = output / mask
            refname = "reference.vcf.gz" if policy["reference"] == "full" else "finecap.vcf.gz"
            run(["java", "-Xmx6g", "-jar", jar, f"gt={folder / f'{mask}.evaluation.vcf.gz'}",
                 f"ref={folder / refname}", f"ref-panel={folder / (policy['reference'] + '.panel')}",
                 f"map={metadata['map']}", f"out={prefix}", "nthreads=1", "seed=20261022",
                 f"array={str(policy['array']).lower()}", f"em={str(policy['em']).lower()}",
                 "update-p=false", "probs=true", "min-maf=0", "min-mac=0"], output / f"{mask}.flare")
            from screen_model_families import read_model
            source_order = read_model(str(prefix) + ".model")["ancestries"]
            values = physical_flare(str(prefix) + ".anc.vcf.gz", source_order, keys, metadata["span"],
                                    metadata["evaluation"])
        elif policy["family"] == "gnomix":
            values = fit_gnomix(folder, metadata, policy, mask, output)
        else:
            if policy["input"] == "production_phased_imputed":
                dosage = np.load(folder / f"{mask}.evaluation.npz")["haplotypes"].reshape(-1, 2, len(keys)).sum(axis=1)
                frequencies = np.load(folder / "reference.npz")["frequencies"]
            else:
                data = np.load(folder / f"{mask}.observed.npz")
                dosage, frequencies = data["dosage"], data["frequencies"]
            values = np.stack([constrained_frequencies(frequencies[row >= 0], row[row >= 0]) for row in dosage])
        for identity, value in zip(metadata["evaluation"], values):
            validate_vector(value)
            predictions.append(dict(candidate=candidate, mask=mask, chromosome=metadata["chromosome"],
                                    sample_id=identity, span_bp=metadata["span"][1] - metadata["span"][0],
                                    seconds=time.monotonic() - start,
                                    **{source: float(value[index]) for index, source in enumerate(SOURCES)}))
        print(f"{candidate}/{mask}/chr{metadata['chromosome']}: eight-source predictions complete.", flush=True)
    pd.DataFrame(predictions).to_csv(output / "predictions.tsv", sep="\t", index=False)
    write_json(output / "complete.json", dict(
        candidate=candidate, chromosome=metadata["chromosome"], rows=len(predictions),
        predictions_sha256=digest(output / "predictions.tsv"), panel_sha256=digest(folder / "panel.json"),
        cohort_sha256=digest(cohort_path), source_order=SOURCES, production_changed=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("panels", "evaluate"))
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--masks", nargs="+")
    parser.add_argument("--folders", type=Path, nargs="+")
    parser.add_argument("--merged", type=Path, nargs="+")
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--candidate")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.stage == "panels":
            panels(args.cohort, args.chromosome, args.masks, args.folders, args.merged, args.output, args.truth)
        else:
            evaluate(args.cohort, args.folder, args.candidate, args.output)
