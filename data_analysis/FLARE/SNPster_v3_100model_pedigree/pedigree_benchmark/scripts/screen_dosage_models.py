#!/usr/bin/env python3
"""Small phase-independent global ancestry screen on physically thinned dosages."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import re
import shutil
import subprocess
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as provenance
from audit_real_transfer import frequency_fits, thin_sites
from screen_model_families import score, verify
from screen_v1_decoders import predict_decoder, training_weights
from train_matched_global import training_table


PROJECT = Path(__file__).resolve().parents[2]
SOURCES = list(adapter.SOURCES)
MODES = ("frequency_simplex_100kb", "dosage_pca64_ridge_clean",
         "dosage_pca64_ridge_matched", "dosage_pca64_trees_matched")
SPACING = 100_000
COMPONENTS = 64


def genotype_dosage(text):
    if (len(text) != 3 or text[1] not in ("|", "/")
            or text[0] not in ("0", "1") or text[2] not in ("0", "1")):
        raise ValueError(f"Complete biallelic diploid genotype required, got {text!r}")
    return (text[0] == "1") + (text[2] == "1")


def load_dosages(path, ids, keys, regions, check_phase=False):
    if len(set(ids)) != len(ids) or not ids:
        raise ValueError("Empty/duplicate dosage sample axis")
    wanted = {key[1]: (index, key) for index, key in enumerate(keys)}
    if len(wanted) != len(keys):
        raise ValueError("Dosage marker positions must be unique")
    values = np.empty((len(ids), len(keys)), dtype=np.int8)
    altered = np.empty_like(values) if check_phase else None
    seen, reorder = set(), None
    command = ["bcftools", "query", "-H", "-R", str(regions), "-s", ",".join(ids),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                selected = [re.sub(r"^\[\d+\]", "", value).removesuffix(":GT")
                            for value in fields[4:]]
                if len(selected) != len(ids) or set(selected) != set(ids):
                    raise ValueError("bcftools dosage sample header mismatch")
                reorder = [selected.index(sample) for sample in ids]
                continue
            if reorder is None or len(fields) != 4 + len(ids):
                raise ValueError("Missing sample header or malformed genotype row")
            key = (fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3])
            if key[1] not in wanted:
                raise ValueError("Position query returned an unrequested dosage site")
            index, expected = wanted[key[1]]
            if key != expected or index in seen:
                raise ValueError(f"Duplicate or incompatible dosage marker: {key} vs {expected}")
            seen.add(index)
            genotypes = [fields[4 + i] for i in reorder]
            values[:, index] = [genotype_dosage(gt) for gt in genotypes]
            if check_phase:
                # Reverse alternate orientations and erase phase, across every held-out genotype.
                changed = [gt[::-1].replace("|", "/") if (i + index) % 2
                           else gt.replace("|", "/") for i, gt in enumerate(genotypes)]
                altered[:, index] = [genotype_dosage(gt) for gt in changed]
        errors = process.stderr.read()
        if process.wait():
            raise ValueError(f"Dosage query failed for {path}: {errors.strip()}")
    finally:
        if process.poll() is None:
            process.terminate()
        process.wait()
        process.stdout.close()
        process.stderr.close()
    if len(seen) != len(keys):
        raise ValueError(f"Missing {len(keys) - len(seen)} required dosage markers")
    if check_phase and not np.array_equal(values, altered):
        raise ValueError("Phase changes altered diploid dosage")
    return values, altered


def sample_batches(cohort, matched, chrom):
    groups = [("reference", cohort["reference"]["vcf"],
               adapter.sample_records(cohort["reference"]["labels_tsv"])),
              ("clean_training", cohort["pure_train"]["vcf"],
               adapter.sample_records(cohort["pure_train"]["labels_tsv"]))]
    groups += [("clean_training", batch["vcf"], adapter.sample_records(batch["samples_tsv"]))
               for batch in cohort["pedigree_train"]["batches"]]
    groups += [(name, batch["vcf"], adapter.sample_records(batch["samples_tsv"]))
               for name, group in cohort["evaluation"].items() for batch in group["batches"]]
    plan = json.loads((matched / "pilot.json").read_text())
    groups += [("imputed_training" if row["role"] == "training" else "array_imputed",
                str(matched / f"chr{chrom}/{row['sample_id']}/truth_panel.vcf.gz"),
                [{"sample_id": row["sample_id"]}]) for row in plan["samples"]]
    return groups


def configuration(out):
    result = json.loads((out / "experiment.json").read_text())
    provenance.verify_hashes(result["implementation_sha256"])
    if result["models"] != list(MODES) or result["spacing_bp"] != SPACING:
        raise ValueError("Frozen dosage experiment configuration changed")
    return result


def complete(path, out):
    if not path.exists():
        return False
    audit = json.loads(path.read_text())
    if (audit["status"] != "complete"
            or audit["experiment_sha256"] != adapter.digest(out / "experiment.json")):
        raise ValueError(f"Incomplete/incompatible dosage checkpoint: {path}")
    provenance.verify_hashes(audit["input_sha256"])
    provenance.verify_hashes(audit["output_sha256"])
    print(f"Verified {path}", flush=True)
    return True


def prepare(screen, matched, out):
    if (out / "experiment.json").exists():
        raise ValueError("Dosage preparation requires a new output")
    cohort = verify(screen)
    counts = adapter.validate_splits(cohort)
    previous = json.loads((matched / "completion.json").read_text())
    checked = json.loads((matched / "verification/final_verification.json").read_text())
    if (previous["status"] != "complete" or checked["status"] != "complete_verified"
            or checked["completion_sha256"] != adapter.digest(matched / "completion.json")):
        raise ValueError("Previously verified matched imputation results changed")
    provenance.verify_hashes(previous["input_sha256"])
    provenance.verify_hashes(previous["output_sha256"])
    inputs = adapter.input_hashes(cohort)
    for chrom in range(1, 23):
        path = matched / f"chr{chrom}/audit.json"
        audit = json.loads(path.read_text())
        if audit["status"] != "complete":
            raise ValueError("Matched imputation is incomplete")
        inputs[str(path)] = adapter.digest(path)
        for path, expected in audit["output_sha256"].items():
            if path.endswith("truth_panel.vcf.gz") or path.endswith("truth_panel.vcf.gz.tbi"):
                inputs[path] = expected
    paths = [Path(__file__), Path(adapter.__file__),
             PROJECT / "pedigree_benchmark/scripts/audit_real_transfer.py",
             PROJECT / "pedigree_benchmark/scripts/pilot_matched_imputation.py",
             PROJECT / "pedigree_benchmark/scripts/train_matched_global.py",
             PROJECT / "pedigree_benchmark/scripts/screen_v1_decoders.py",
             PROJECT / "pedigree_benchmark/scripts/screen_model_families.py",
             PROJECT / "config/unified_selection.json", PROJECT / "config/real_benchmarks.tsv"]
    inputs.update(provenance.fingerprint_files([
        screen / "cohort.json", matched / "completion.json", matched / "pilot.json",
        matched / "experiment.json", screen / "baseline/metrics.json"]))
    provenance.verify_hashes(inputs)
    settings = dict(
        status="configured", models=list(MODES), spacing_bp=SPACING,
        thinning="First frozen common marker in each100kb physical bin; not r-squared LD pruning",
        n_original_common_markers=197260, seed=adapter.SEED, sample_counts=counts,
        n_components=COMPONENTS, reference_only_preprocessing=True,
        pca=dict(svd_solver="randomized", whiten=True, random_state=adapter.SEED,
                 iterated_power=3),
        scale="Reference allele count mean and sqrt(2*p*(1-p))",
        ridge=dict(alpha=.01, fit_intercept=True),
        trees=dict(n_estimators=200, min_samples_leaf=5, max_features=1., random_state=adapter.SEED),
        matched_training_weight=.25, clean_training_weights="Seven-source/family-balanced, half each",
        no_evaluation_truth_in_fit=True, missing_genotypes="Explicit failure; no missing-call padding",
        phase_invariance="All held-out GTs are reloaded with phase erased and alternating allele orientation reversed",
        implementation_sha256=provenance.fingerprint_files(paths), input_sha256=inputs,
        local_ancestry_models_unchanged=True, production_promotion=False,
        limitations=[
            "Physically thinned markers, not formal LD pruning.",
            "Previously inspected development cohorts, not independent certification.",
            "Only seven actually-imputed training genomes and five related/shared-donor validation controls.",
            "Phase invariance does not remove genotype-call errors, reference bias or label limitations.",
            "Hard GT allele counts, not posterior-mean FORMAT/DS dosages.",
            "One array mask; production GP/DS-QC missingness and other arrays are not covered."
        ])
    adapter.write_json(out / "experiment.json", settings)
    shutil.copyfile(screen / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(screen / "baseline/metrics.json", out / "baseline/metrics.json")
    print("Prepared four fixed phase-independent genotype models; no imputation needed", flush=True)


def extract(screen, matched, out, chrom):
    settings = configuration(out)
    directory = out / f"chr{chrom}"
    if complete(directory / "audit.json", out):
        return
    if directory.exists():
        raise ValueError("Incomplete dosage directory exists; inspect before retry")
    directory.mkdir(parents=True)
    cohort = json.loads((screen / "cohort.json").read_text())
    all_keys = []
    for line in Path(cohort["markers_tsv"]).read_text().splitlines():
        chm, pos, ref, alt = line.split("\t")
        if chm.removeprefix("chr") == str(chrom):
            all_keys.append((str(chrom), int(pos), ref, alt))
    if not all_keys or all_keys != sorted(all_keys, key=lambda key: key[1]):
        raise ValueError("Missing or unordered frozen marker identities")
    mask = thin_sites(all_keys, np.ones(len(all_keys), dtype=bool), SPACING)
    keys = [key for key, selected in zip(all_keys, mask) if selected]
    regions = directory / "positions.tsv"
    regions.write_text("".join(f"{chrom}\t{key[1]}\n" for key in keys))
    batches = sample_batches(cohort, matched, chrom)
    paths = [Path(path) for _, path, _ in batches]
    paths += [Path(cohort["markers_tsv"]), screen / "cohort.json", matched / "pilot.json"]
    expected = {str(path): settings["input_sha256"][str(path)] for path in paths}
    provenance.verify_hashes(expected)
    axes, matrix, changed, phase_rows = [], [], [], []
    for group, path, rows in batches:
        ids = [row["sample_id"] for row in rows]
        evaluation = group in ("pure", "pedigree", "real", "array_imputed")
        x, altered = load_dosages(path, ids, keys, regions, evaluation)
        first = len(axes)
        axes.extend(dict(group=group, sample_id=sample) for sample in ids)
        matrix.append(x)
        if evaluation:
            changed.append(altered)
            phase_rows.extend(range(first, first + len(ids)))
    identities = pd.DataFrame(axes)
    if len(identities) != 4541 or identities.duplicated(["group", "sample_id"]).any():
        raise ValueError("Unexpected dosage sample-axis coverage")
    np.savez_compressed(directory / "dosages.npz", X=np.concatenate(matrix),
                        phase_X=np.concatenate(changed), phase_rows=np.array(phase_rows))
    identities.to_csv(directory / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(keys, columns=["chrom", "pos", "ref", "alt"]).to_csv(
        directory / "markers.tsv", sep="\t", index=False)
    provenance.verify_hashes(expected)
    audit = dict(status="complete", chromosome=chrom, n_original_markers=len(all_keys),
                 n_markers=len(keys), n_samples=len(identities), phase_checked_individuals=len(phase_rows),
                 phase_checked_genotypes=len(phase_rows) * len(keys),
                 experiment_sha256=adapter.digest(out / "experiment.json"), input_sha256=expected,
                 output_sha256=provenance.fingerprint_files([
                     directory / name for name in ("dosages.npz", "samples.tsv", "markers.tsv")]))
    adapter.write_json(directory / "audit.json", audit)
    print(f"chr{chrom}: {len(keys)} physically thinned SNPs, 4541 samples; phase invariance checked",
          flush=True)


def standardize(x, mean, scale):
    x = np.asarray(x)
    if (x.ndim != 2 or x.shape[1] != len(mean) or len(mean) != len(scale)
            or not np.all(np.isfinite(x)) or np.any(x < 0) or np.any(x > 2)
            or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(scale)) or np.any(scale <= 0)):
        raise ValueError("Invalid complete dosage matrix or reference scaling")
    return ((x - mean) / scale).astype(np.float32)


def train_heads(scores, axes, cohort):
    pure, ped = training_table(cohort)
    metadata = pd.concat([pure, ped], ignore_index=True)
    clean_index = axes.index[axes.group.eq("clean_training")]
    imputed_index = axes.index[axes.group.eq("imputed_training")]
    clean_ids, imputed_ids = axes.loc[clean_index, "sample_id"], axes.loc[imputed_index, "sample_id"]
    if set(clean_ids) != set(metadata.sample_id) or len(clean_ids) != 799:
        raise ValueError("Clean PCA training identities changed")
    if len(imputed_ids) != 7 or not set(imputed_ids) <= set(ped.sample_id):
        raise ValueError("Imputed PCA training contains a nontraining individual")
    truth = ["truth_" + source for source in SOURCES]
    y = metadata.set_index("sample_id").loc[clean_ids, truth].to_numpy(float)
    extra_y = ped.set_index("sample_id").loc[imputed_ids, truth].to_numpy(float)
    weights = pd.Series(training_weights(pure, ped), index=metadata.sample_id).loc[clean_ids].to_numpy()
    clean_x, extra_x = scores[clean_index].astype(float), scores[imputed_index].astype(float)
    a, b = np.r_[clean_x, extra_x], np.r_[y, extra_y]
    w = np.r_[.75 * weights, np.repeat(.25 / 7, 7)]
    return {
        MODES[1]: Ridge(alpha=.01, fit_intercept=True).fit(clean_x, y, sample_weight=weights),
        MODES[2]: Ridge(alpha=.01, fit_intercept=True).fit(a, b, sample_weight=w),
        MODES[3]: ExtraTreesRegressor(
            n_estimators=200, min_samples_leaf=5, max_features=1.,
            random_state=adapter.SEED, n_jobs=1).fit(a, b, sample_weight=w),
    }


def predict(bundle, x):
    if bundle["kind"] == "frequency":
        values, _ = frequency_fits(bundle["frequencies"], x)
        return values
    x = standardize(x, bundle["mean"], bundle["scale"])
    scores = bundle["pca"].transform(x).astype(float)
    return predict_decoder(bundle["head"], scores)


def fit(screen, matched, out):
    settings = configuration(out)
    if (out / "completion.json").exists():
        raise ValueError("Completed dosage results must not be overwritten")
    arrays, phase_arrays, markers, axes, phase_rows = [], [], [], None, None
    total_original = 0
    for chrom in range(1, 23):
        directory = out / f"chr{chrom}"
        if not complete(directory / "audit.json", out):
            raise ValueError(f"Missing dosage chromosome{chrom}")
        audit = json.loads((directory / "audit.json").read_text())
        total_original += audit["n_original_markers"]
        current = pd.read_csv(directory / "samples.tsv", sep="\t")
        with np.load(directory / "dosages.npz") as saved:
            if axes is not None and (not axes.equals(current)
                                      or not np.array_equal(phase_rows, saved["phase_rows"])):
                raise ValueError("Dosage sample/phase axes differ between chromosomes")
            arrays.append(saved["X"])
            phase_arrays.append(saved["phase_X"])
            phase_rows = saved["phase_rows"]
        axes = current
        marker_frame = pd.read_csv(directory / "markers.tsv", sep="\t", dtype={"chrom": str})
        markers.extend(list(marker_frame.itertuples(index=False, name=None)))
    if total_original != 197260 or len(markers) != 26066 or len(set(markers)) != len(markers):
        raise ValueError("Unexpected physical thinning/marker coverage")
    x, phase_x = np.concatenate(arrays, axis=1), np.concatenate(phase_arrays, axis=1)
    if not np.array_equal(x[phase_rows], phase_x) or len(phase_rows) != 288:
        raise ValueError("Whole-genome phase changes altered dosage inputs")
    cohort = json.loads((screen / "cohort.json").read_text())
    reference_rows = adapter.sample_records(cohort["reference"]["labels_tsv"], True)
    ref_index = axes.index[axes.group.eq("reference")]
    if axes.loc[ref_index, "sample_id"].tolist() != [row["sample_id"] for row in reference_rows]:
        raise ValueError("Reference PCA identity axis differs from labels")
    reference = x[ref_index]
    mean = reference.mean(axis=0, dtype=np.float64)
    p = mean / 2
    scale = np.sqrt(2 * p * (1 - p))
    frequencies = np.column_stack([
        reference[[i for i, row in enumerate(reference_rows) if row["source_group"] == source]]
        .mean(axis=0, dtype=np.float64) / 2 for source in SOURCES])
    pca = PCA(n_components=COMPONENTS, svd_solver="randomized", whiten=True,
              random_state=adapter.SEED, iterated_power=3)
    with threadpool_limits(limits=6):
        pca.fit(standardize(reference, mean, scale))
        scores = pca.transform(standardize(x, mean, scale))
        models = train_heads(scores, axes, cohort)
    scores_frame = axes.copy()
    for index in range(COMPONENTS):
        scores_frame[f"PC{index + 1}"] = scores[:, index]
    scores_frame.to_csv(out / "pca_scores.tsv", sep="\t", index=False)
    results, outputs = {}, [out / "pca_scores.tsv"]
    start = time.monotonic()
    for mode in MODES:
        directory = out / "candidates" / mode
        directory.mkdir(parents=True)
        bundle = dict(
            kind="frequency" if mode == MODES[0] else "pca", marker_keys=markers,
            sources=SOURCES, experiment_sha256=adapter.digest(out / "experiment.json"),
            frequency_informative_range=.03, phase_independent=True,
            reference_only_preprocessing=True)
        if mode == MODES[0]:
            bundle["frequencies"] = frequencies
        else:
            bundle.update(mean=mean, scale=scale, pca=pca, head=models[mode])
        path = directory / "model.pkl"
        with path.with_suffix(".pending").open("wb") as stream:
            pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
        path.with_suffix(".pending").replace(path)
        with path.open("rb") as stream:
            restored = pickle.load(stream)
        with threadpool_limits(limits=6):
            before = predict(bundle, x[phase_rows])
            after = predict(restored, x[phase_rows])
            altered = predict(restored, phase_x)
        reload_delta = float(np.max(np.abs(before - after)))
        phase_delta = float(np.max(np.abs(after - altered)))
        if (reload_delta > 1e-10 or phase_delta > 1e-10
                or not np.all(np.isfinite(after)) or np.any(after < -1e-10)
                or not np.allclose(after.sum(axis=1), 1, atol=1e-8)):
            raise ValueError("Invalid probabilities, saved model mismatch or phase dependence")
        predictions = axes.loc[phase_rows].rename(columns={"group": "dataset"}).copy()
        for index, source in enumerate(SOURCES):
            predictions["pred_" + source] = after[:, index]
        array = predictions[predictions.dataset.eq("array_imputed")].copy()
        array["TV_distance"] = np.abs(array[["pred_" + s for s in SOURCES]].to_numpy()
                                      - [.5, .5, 0, 0, 0, 0, 0]).sum(axis=1) / 2
        if len(array) != 5:
            raise ValueError("Incomplete imputed validation axis")
        array.to_csv(directory / "array_imputed_predictions.tsv", sep="\t", index=False)
        predictions[~predictions.dataset.eq("array_imputed")].to_csv(
            directory / "predictions.tsv", sep="\t", index=False)
        score(out, mode)
        metrics = json.loads((directory / "metrics.json").read_text())
        results[mode] = dict(
            pure_mean_tv=metrics["pure"]["mean_tv"], pedigree_mean_tv=metrics["pedigree"]["mean_tv"],
            array_imputed_mean_tv=float(array.TV_distance.mean()), real=metrics["real"],
            unchanged_gates_pass=metrics["all_gates_pass"],
            failed_gates=[key for key, value in metrics["gates"].items() if not value],
            max_roundtrip_prediction_error=reload_delta, max_phase_prediction_error=phase_delta)
        outputs += [directory / name for name in
                    ("model.pkl", "array_imputed_predictions.tsv", "predictions.tsv", "metrics.json")]
    previous = json.loads((matched / "completion.json").read_text())
    pd.DataFrame([dict(candidate=key, **{k: result[k] for k in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv", "unchanged_gates_pass")})
        for key, result in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    provenance.verify_hashes(settings["input_sha256"])
    summary = dict(
        status="complete", n_candidates=4, n_markers=len(markers), n_chromosomes=22,
        n_reference=3447, n_clean_training=799, n_imputed_training=7,
        evaluation_n=dict(pure=159, pedigree=123, real=1, array_imputed=5),
        pca_reference_explained_variance=float(pca.explained_variance_ratio_.sum()),
        results=results, previous_haplotype_model_results=previous["results"],
        reference_only_preprocessing=True, phase_invariance_checked_on_all_heldout_genotypes=True,
        production_promotion=False, independent_final_test=False,
        limitations=settings["limitations"], final_prediction_seconds=time.monotonic() - start,
        input_sha256=provenance.fingerprint_files([out / "experiment.json"]),
        output_sha256=provenance.fingerprint_files(outputs + [out / "comparison.tsv"]))
    adapter.write_json(out / "completion.json", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "extract", "fit"))
    parser.add_argument("--screen", type=Path, default=Path("/results"))
    parser.add_argument("--matched", type=Path, default=Path("/pilot"))
    parser.add_argument("--outdir", type=Path, default=Path("/dosage"))
    parser.add_argument("--chromosome", type=int, choices=range(1, 23))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.screen, args.matched, args.outdir)
    elif args.stage == "fit":
        fit(args.screen, args.matched, args.outdir)
    elif args.chromosome is None:
        parser.error("--chromosome is required for extraction")
    else:
        extract(args.screen, args.matched, args.outdir, args.chromosome)


if __name__ == "__main__":
    main()
