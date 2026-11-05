#!/usr/bin/env python3
"""Matched-imputation training for a small global decoder, not local ancestry."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import shutil
import subprocess

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
import pilot_matched_imputation as pilot
from screen_model_families import score, verify
from screen_v1_decoders import predict_decoder, training_weights


PROJECT = Path(__file__).resolve().parents[2]
FEATURES = [f"{stage}_{source}" for stage in ("base", "smooth") for source in adapter.SOURCES]
ARMS = ("global_ridge_clean", "global_ridge_matched")
N_SELECTED = 7


def truth_vector(row):
    return np.array([float(row["realized_EAS"]) + float(row["realized_SEA"])
                     if source == "EAS_SEA" else float(row["realized_" + source])
                     for source in adapter.SOURCES])


def select_training(candidates):
    eligible = [row for row in candidates if len(row["ancestors"]) == 2
                and np.count_nonzero(np.isclose(row["truth"], .5, atol=1e-12, rtol=0)) == 2
                and np.count_nonzero(np.asarray(row["truth"]) > 1e-12) == 2]
    eligible.sort(key=lambda row: row["sample_id"])
    if not eligible:
        raise ValueError("No dense-founder-available two-source training descendants")
    matrix = [np.ones(len(eligible))]
    lower, upper = [N_SELECTED], [N_SELECTED]
    for index in range(len(adapter.SOURCES)):
        matrix.append([float(row["truth"][index] > .1) for row in eligible])
        lower.append(2)
        upper.append(2)
    for key in ("family_id", "ancestors"):
        identities = sorted({identity for row in eligible
                             for identity in (row[key] if key == "ancestors" else [row[key]])})
        for identity in identities:
            matrix.append([float(identity in row[key] if key == "ancestors"
                                 else identity == row[key]) for row in eligible])
            lower.append(0)
            upper.append(1)
    cost = [int.from_bytes(hashlib.sha256(
        f"{adapter.SEED}:{row['sample_id']}".encode()).digest()[:8], "little") / 2**64
        for row in eligible]
    result = milp(np.asarray(cost), integrality=np.ones(len(eligible)),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(
                      np.asarray(matrix), np.asarray(lower), np.asarray(upper)))
    if not result.success or result.x is None:
        raise ValueError(f"Balanced donor-disjoint training selection is infeasible: {result.message}")
    selected = [row for row, value in zip(eligible, result.x) if value > .5]
    if (len(selected) != N_SELECTED
            or len({row["family_id"] for row in selected}) != N_SELECTED
            or len({donor for row in selected for donor in row["ancestors"]}) != 2 * N_SELECTED
            or not np.allclose(np.sum([row["truth"] for row in selected], axis=0), 1)):
        raise ValueError("Selected cohort does not satisfy exact balance and disjointness")
    return selected


def frozen_config(out):
    config = json.loads((out / "experiment.json").read_text())
    pilot.verify_hashes(config["implementation_sha256"])
    if config["pilot_sha256"] != adapter.digest(out / "pilot.json"):
        raise ValueError("Matched cohort configuration changed")
    return config


def checkpoint(path, out):
    if not path.exists():
        return False
    record = json.loads(path.read_text())
    if (record["status"] != "complete"
            or record["experiment_sha256"] != adapter.digest(out / "experiment.json")):
        raise ValueError(f"Incomplete/incompatible checkpoint: {path}")
    pilot.verify_hashes(record["input_sha256"])
    pilot.verify_hashes(record["output_sha256"])
    print(f"Verified {path}", flush=True)
    return True


def prepare(screen, previous, out, reference):
    if (out / "experiment.json").exists() or (out / "pilot.json").exists():
        raise ValueError("Preparation requires a fresh experiment")
    cohort = verify(screen)
    adapter.validate_splits(cohort)
    old = json.loads((previous / "pilot.json").read_text())
    training_root = PROJECT / "results_v3_pedigree_training_screen_01/cohort/pedigree"
    development_root = Path(old["plan_root"])
    dense = set(subprocess.check_output([
        "bcftools", "query", "-l",
        str(reference / "beagle_references" / pilot.REFERENCE_PATTERN.format(chrom=22))
    ], text=True).splitlines())
    candidates, batches, inputs = [], {}, []
    for role, root, groups in (
            ("training", training_root, cohort["pedigree_train"]["batches"]),
            ("evaluation", development_root, cohort["evaluation"]["pedigree"]["batches"])):
        for batch in groups:
            name = f"{role}_{batch['batch']}"
            plan_path = root / f"plans/plans/{batch['batch']}/pedigrees.json"
            plan = json.loads(plan_path.read_text())
            batches[name] = dict(plan=str(plan_path), vcf=batch["vcf"],
                                 original_batch=batch["batch"], original_root=str(root))
            inputs.extend([plan_path, Path(batch["vcf"]), Path(batch["samples_tsv"])])
            if role == "training":
                for row in adapter.records(batch["samples_tsv"]):
                    donors = pilot.ancestor_donors(plan, row["sample_id"])
                    if donors <= dense:
                        candidates.append(dict(sample_id=row["sample_id"], family_id=row["family_id"],
                                               batch=name, ancestors=sorted(donors),
                                               truth=truth_vector(row).tolist(), role="training"))
    selected = select_training(candidates)
    samples = list(selected)
    for row in old["samples"]:
        samples.append({**row, "batch": "evaluation_" + row["batch"], "role": "evaluation",
                        "truth": [.5, .5, 0, 0, 0, 0, 0]})
    pure_train = {row["sample_id"] for row in adapter.records(cohort["pure_train"]["labels_tsv"])}
    train_donors = {donor for row in selected for donor in row["ancestors"]}
    heldout = set(old["excluded_reference_ids"])
    ancestry_reference = {row["sample_id"] for row in adapter.records(cohort["reference"]["labels_tsv"])}
    if (not train_donors <= pure_train or train_donors & (heldout | ancestry_reference)
            or pure_train & heldout):
        raise ValueError("Training founder/reference/evaluation identity separation failed")
    excluded = heldout | pure_train
    used_batches = {row["batch"] for row in samples}
    batches = {name: batch for name, batch in batches.items() if name in used_batches}
    bp_dir = out / "pedigree/chromosomes"
    bp_dir.mkdir(parents=True)
    for chrom in range(1, 23):
        for name, batch in batches.items():
            source = Path(batch["original_root"]) / (
                f"chromosomes/pdb_{batch['original_batch']}_chr{chrom}.bp")
            target = bp_dir / f"pdb_{name}_chr{chrom}.bp"
            target.symlink_to(source)
            inputs.extend([source, target])
    for name, expected in (
            ("real_observed_keys.tsv", old["prepared_output_sha256"]["/pilot/real_observed_keys.tsv"]),
            ("beagle.jar", old["input_sha256"]["/pilot/beagle.jar"])):
        source = previous / name
        if adapter.digest(source) != expected:
            raise ValueError(f"Original matched pilot artifact changed: {source}")
        shutil.copyfile(source, out / name)
    (out / "excluded_reference_ids.txt").write_text("".join(s + "\n" for s in sorted(excluded)))
    pd.DataFrame(samples).to_csv(out / "selected_samples.tsv", sep="\t", index=False)
    inputs += [Path(pilot.__file__), Path(adapter.__file__), screen / "cohort.json",
               Path(cohort["markers_tsv"]), previous / "pilot.json", out / "beagle.jar"]
    for chrom in range(1, 23):
        ref = reference / "beagle_references" / pilot.REFERENCE_PATTERN.format(chrom=chrom)
        inputs += [ref, Path(str(ref) + ".tbi"),
                   reference / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map"]
        model = screen / f"candidates/{pilot.CANDIDATE}/models/chr{chrom}.pkl"
        inputs += [model, model.with_suffix(".json")]
    plan = dict(
        status="prepared", candidate=pilot.CANDIDATE, seed=adapter.SEED,
        reference=str(reference), beagle=str(out / "beagle.jar"),
        samples=samples, donor_ids=sorted({d for row in samples for d in row["ancestors"]}),
        batches=batches, plan_root=str(out / "pedigree"),
        excluded_reference_ids=sorted(excluded), input_sha256=pilot.fingerprint_files(inputs),
        prepared_output_sha256=pilot.fingerprint_files([
            out / "real_observed_keys.tsv", out / "excluded_reference_ids.txt",
            out / "selected_samples.tsv"]),
        singleton_imputation=True, training_and_evaluation_reference_policy_identical=True,
        production_changed=False, no_new_families=True)
    adapter.write_json(out / "pilot.json", plan)
    implementation = pilot.fingerprint_files([
        Path(__file__), Path(adapter.__file__), Path(pilot.__file__),
        PROJECT / "pedigree_benchmark/scripts/screen_v1_decoders.py",
        PROJECT / "pedigree_benchmark/scripts/screen_model_families.py"])
    experiment = dict(
        status="configured", model="global Ridge, alpha=.001, fit_intercept=False",
        arms=list(ARMS), global_only=True, local_models_frozen=True,
        features=FEATURES, imputed_training_weight=.25,
        matched_group_weights=dict(clean_pure=.375, clean_pedigree=.375, imputed_training=.25),
        pure_weights="Seven-source balanced", pedigree_weights="Family balanced",
        imputed_weights="Equal per individual; each of seven sources totals exactly one genome",
        n_imputed_training=N_SELECTED, n_array_validation=5,
        selected_training_ids=[row["sample_id"] for row in selected],
        pure_training_n=307, clean_pedigree_training_n=492,
        no_evaluation_truth_in_fit=True, no_haplotype_labels_after_imputation=True,
        implementation_sha256=implementation, pilot_sha256=adapter.digest(out / "pilot.json"),
        no_production_promotion=True,
        limitations=[
            "Previously inspected development cohorts, not independent final certification.",
            "Identity separation does not establish kinship separation.",
            "Seven imputed training genomes cannot establish broad population robustness.",
            "Observed array mask comes from one real individual.",
            "Direct complete Beagle genotypes; does not reproduce production GP/DS QC missingness.",
            "Expanded exclusions change the old imputation reference; validation is rerun under the same new policy.",
            "Two global heads only: clean-only control and matched-data primary; no local ancestry model is retrained."
        ])
    adapter.write_json(out / "experiment.json", experiment)
    shutil.copyfile(screen / "cohort.json", out / "cohort.json")
    (out / "baseline").mkdir()
    shutil.copyfile(screen / "baseline/metrics.json", out / "baseline/metrics.json")
    print(f"Prepared {N_SELECTED} balanced donor-disjoint training pedigrees and five validation controls",
          flush=True)
    print(pd.DataFrame(selected)[["sample_id", "ancestors", "truth"]].to_string(index=False), flush=True)


def features(screen, out, chrom):
    frozen_config(out)
    directory = out / f"features/chr{chrom}"
    if checkpoint(directory / "audit.json", out):
        return
    if directory.exists():
        raise ValueError("Incomplete feature directory exists; inspect before retry")
    audit_path = out / f"chr{chrom}/audit.json"
    imputation_audit = json.loads(audit_path.read_text())
    if imputation_audit["status"] != "complete":
        raise ValueError("Imputation must complete before feature extraction")
    pilot.verify_hashes(imputation_audit["output_sha256"])
    pilot.verify_hashes(imputation_audit["input_sha256"])
    directory.mkdir(parents=True)
    cohort = json.loads((screen / "cohort.json").read_text())
    plan = json.loads((out / "pilot.json").read_text())
    adapter.load_upstream()
    bundle, _ = adapter.saved_chromosome(cohort, adapter.config_for(pilot.CANDIDATE), chrom,
                                         screen / "candidates" / pilot.CANDIDATE)
    model, layout, keys = bundle["model"], bundle["layout"], bundle["keys"]
    rows, inputs = [], [audit_path]

    def append(group, ids, b):
        smoothed = adapter.ordered_probabilities(model, b)
        base, _ = adapter.aggregate_probabilities(b, layout)
        smooth, _ = adapter.aggregate_probabilities(smoothed, layout)
        combined = np.concatenate([base, smooth], axis=1) / len(keys)
        for sample, values in zip(ids, combined):
            rows.append(dict(group=group, sample_id=sample, chromosome=chrom,
                             n_markers=len(keys), **dict(zip(FEATURES, values))))

    with threadpool_limits(limits=1):
        cache = Path(f"/robust/profiles/chr{chrom}")
        prior = json.loads((cache / "audit.json").read_text())
        if prior["status"] != "complete":
            raise ValueError("Original clean-training profile checkpoint is incomplete")
        pilot.verify_hashes(prior["input_sha256"])
        pilot.verify_hashes(prior["output_sha256"])
        ids = json.loads((cache / "samples.json").read_text())
        if len(ids) != 799:
            raise ValueError("Original training sample axis changed")
        with np.load(cache / "original.npz") as data:
            append("clean_training", ids, data["probabilities"])
        inputs += [cache / "audit.json", cache / "original.npz", cache / "samples.json"]
        for track, group in cohort["evaluation"].items():
            for batch in group["batches"]:
                inputs += [Path(batch["vcf"]), Path(batch["samples_tsv"])]
                for part in adapter.chunks(adapter.sample_records(batch["samples_tsv"]), 64):
                    ids = [row["sample_id"] for row in part]
                    x = adapter.load_haplotypes(batch["vcf"], ids, chrom, keys)
                    append(track, ids, model.base.predict_proba(x))
        for sample in plan["samples"]:
            path = out / f"chr{chrom}/{sample['sample_id']}/truth_panel.vcf.gz"
            inputs.append(path)
            x = adapter.load_haplotypes(path, [sample["sample_id"]], chrom, keys)
            append("imputed_training" if sample["role"] == "training" else "array_imputed",
                   [sample["sample_id"]], model.base.predict_proba(x))
    frame = pd.DataFrame(rows)
    if len(frame) != 1094 or frame.duplicated(["group", "sample_id"]).any():
        raise ValueError("Unexpected feature sample-axis coverage")
    for stage in ("base", "smooth"):
        block = frame[[f"{stage}_{source}" for source in adapter.SOURCES]].to_numpy()
        if (not np.isfinite(block).all() or np.any(block < 0)
                or not np.allclose(block.sum(axis=1), 1, atol=2e-6)):
            raise ValueError("Invalid global feature probabilities")
    frame.to_csv(directory / "features.tsv", sep="\t", index=False)
    record = dict(status="complete", chromosome=chrom, n_markers=len(keys), n_rows=len(frame),
                  experiment_sha256=adapter.digest(out / "experiment.json"),
                  input_sha256=pilot.fingerprint_files(inputs),
                  output_sha256=pilot.fingerprint_files([directory / "features.tsv"]))
    adapter.write_json(directory / "audit.json", record)
    print(f"chr{chrom}: all {len(frame)} global feature rows checkpointed", flush=True)


def training_table(cohort):
    pure = pd.read_csv(cohort["pure_train"]["labels_tsv"], sep="\t")
    for source in adapter.SOURCES:
        pure[f"truth_{source}"] = pure.source_group.eq(source).astype(float)
    pure["family_id"] = pure.sample_id
    ped = pd.concat([pd.read_csv(b["samples_tsv"], sep="\t")
                     for b in cohort["pedigree_train"]["batches"]], ignore_index=True)
    values = np.stack([truth_vector(row) for row in ped.to_dict("records")])
    for index, source in enumerate(adapter.SOURCES):
        ped[f"truth_{source}"] = values[:, index]
    return pure, ped


def fit_global_heads(clean, imputed, pure, ped):
    truth = ["truth_" + source for source in adapter.SOURCES]
    metadata = pd.concat([pure, ped], ignore_index=True)
    joined = metadata.merge(clean, on="sample_id", how="outer", validate="one_to_one", indicator=True)
    if len(joined) != 799 or not joined._merge.eq("both").all():
        raise ValueError("Clean training metadata/features do not match exactly")
    if len(imputed) != N_SELECTED or not set(imputed.sample_id) <= set(ped.sample_id):
        raise ValueError("Only selected training descendants may enter the imputed fit")
    extra = imputed.merge(ped[["sample_id", *truth]], on="sample_id", validate="one_to_one")
    weights = pd.Series(training_weights(pure, ped), index=metadata.sample_id)
    clean_weights = weights.loc[joined.sample_id].to_numpy()
    x, y = joined[FEATURES].to_numpy(float), joined[truth].to_numpy(float)
    models = {}
    for arm in ARMS:
        if arm == "global_ridge_clean":
            a, b, w = x, y, clean_weights
        else:
            a = np.concatenate([x, extra[FEATURES].to_numpy(float)])
            b = np.concatenate([y, extra[truth].to_numpy(float)])
            w = np.r_[.75 * clean_weights, np.repeat(.25 / len(extra), len(extra))]
        model = Ridge(alpha=.001, fit_intercept=False).fit(a, b, sample_weight=w)
        models[arm] = model
    return models


def summarize(screen, out):
    config = frozen_config(out)
    plan = json.loads((out / "pilot.json").read_text())
    cohort = json.loads((screen / "cohort.json").read_text())
    pilot.verify_hashes(plan["input_sha256"])
    frames = []
    for chrom in range(1, 23):
        path = out / f"features/chr{chrom}"
        if not checkpoint(path / "audit.json", out):
            raise ValueError(f"Missing global features for chromosome{chrom}")
        frame = pd.read_csv(path / "features.tsv", sep="\t")
        frame["chromosome"] = chrom
        frames.append(frame)
    frame = pd.concat(frames, ignore_index=True)
    if frame.duplicated(["group", "sample_id", "chromosome"]).any():
        raise ValueError("Duplicate chromosome features")
    for _, group in frame.groupby(["group", "sample_id"]):
        if set(group.chromosome) != set(range(1, 23)) or group.n_markers.sum() != 197260:
            raise ValueError("Incomplete feature marker/chromosome coverage")
    frame[FEATURES] = frame[FEATURES].mul(frame.n_markers, axis=0)
    total = frame.groupby(["group", "sample_id"])[[*FEATURES, "n_markers"]].sum()
    total[FEATURES] = total[FEATURES].div(total.n_markers, axis=0)
    total = total.reset_index()
    total.to_csv(out / "global_features.tsv", sep="\t", index=False)
    pure, ped = training_table(cohort)
    clean = total[total.group.eq("clean_training")][["sample_id", *FEATURES]]
    imputed = total[total.group.eq("imputed_training")][["sample_id", *FEATURES]]
    if set(imputed.sample_id) != set(config["selected_training_ids"]):
        raise ValueError("Imputed training identities differ from the prepared subset")
    models = fit_global_heads(clean, imputed, pure, ped)
    evaluation = total[total.group.isin(["pure", "pedigree", "real", "array_imputed"])].copy()
    if len(evaluation) != 288:
        raise ValueError("Evaluation features have incorrect coverage")
    frozen = evaluation[["group", "sample_id", *["smooth_" + s for s in adapter.SOURCES]]].rename(
        columns={"group": "dataset", **{f"smooth_{s}": f"pred_{s}" for s in adapter.SOURCES}})
    original = pd.read_csv(screen / f"candidates/{pilot.CANDIDATE}/predictions.tsv", sep="\t")
    check = frozen[~frozen.dataset.eq("array_imputed")].merge(
        original, on=["dataset", "sample_id"], validate="one_to_one", suffixes=("_new", "_old"))
    max_error = float(np.max(np.abs(
        check[[f"pred_{s}_new" for s in adapter.SOURCES]].to_numpy()
        - check[[f"pred_{s}_old" for s in adapter.SOURCES]].to_numpy())))
    if len(check) != 283 or max_error > 1e-10:
        raise ValueError("Original frozen validation predictions did not reproduce")
    outputs, results = [], {}
    for arm in ("frozen_crf", *ARMS):
        directory = out / "candidates" / arm
        directory.mkdir(parents=True, exist_ok=True)
        if arm == "frozen_crf":
            predictions = frozen.copy()
            delta = 0.
        else:
            bundle = dict(model=models[arm], sources=list(adapter.SOURCES), feature_columns=FEATURES,
                          global_only=True, local_models_frozen=True, experiment=arm,
                          experiment_sha256=adapter.digest(out / "experiment.json"))
            model_path = directory / "global_decoder.pkl"
            with model_path.with_suffix(".pending").open("wb") as stream:
                pickle.dump(bundle, stream, protocol=pickle.HIGHEST_PROTOCOL)
            model_path.with_suffix(".pending").replace(model_path)
            with model_path.open("rb") as stream:
                reloaded = pickle.load(stream)
            x = evaluation[FEATURES].to_numpy(float)
            before = predict_decoder(models[arm], x)
            after = predict_decoder(reloaded["model"], x)
            delta = float(np.max(np.abs(before - after)))
            if delta > 1e-10:
                raise ValueError("Serialized decoder changed evaluation predictions")
            predictions = evaluation[["group", "sample_id"]].rename(columns={"group": "dataset"})
            for index, source in enumerate(adapter.SOURCES):
                predictions["pred_" + source] = after[:, index]
            outputs.append(model_path)
        array = predictions[predictions.dataset.eq("array_imputed")].copy()
        array["TV_distance"] = np.abs(array[["pred_" + s for s in adapter.SOURCES]].to_numpy()
                                      - [.5, .5, 0, 0, 0, 0, 0]).sum(axis=1) / 2
        array.to_csv(directory / "array_imputed_predictions.tsv", sep="\t", index=False)
        predictions[~predictions.dataset.eq("array_imputed")].to_csv(
            directory / "predictions.tsv", sep="\t", index=False)
        score(out, arm)
        metrics = json.loads((directory / "metrics.json").read_text())
        results[arm] = dict(
            pure_mean_tv=metrics["pure"]["mean_tv"], pedigree_mean_tv=metrics["pedigree"]["mean_tv"],
            array_imputed_mean_tv=float(array.TV_distance.mean()), real=metrics["real"],
            unchanged_gates_pass=metrics["all_gates_pass"],
            failed_gates=[key for key, value in metrics["gates"].items() if not value],
            max_decoder_roundtrip_error=delta)
        outputs += [directory / name for name in
                    ("array_imputed_predictions.tsv", "predictions.tsv", "metrics.json")]
    pd.DataFrame([dict(candidate=key, **{k: value[k] for k in (
        "pure_mean_tv", "pedigree_mean_tv", "array_imputed_mean_tv", "unchanged_gates_pass")})
        for key, value in results.items()]).to_csv(out / "comparison.tsv", sep="\t", index=False)
    summary = dict(
        status="complete", n_chromosomes=22, n_imputed_training=7, n_array_validation=5,
        n_markers=197260, results=results, local_ancestry_models_retrained=False,
        max_original_clean_validation_error=max_error, production_promotion=False,
        independent_final_test=False, limitations=config["limitations"],
        input_sha256=pilot.fingerprint_files([out / "experiment.json", out / "pilot.json"]),
        output_sha256=pilot.fingerprint_files(outputs + [
            out / "global_features.tsv", out / "comparison.tsv"]))
    adapter.write_json(out / "completion.json", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("prepare", "impute", "features", "summarize"))
    parser.add_argument("--screen", type=Path, default=Path("/results"))
    parser.add_argument("--previous", type=Path, default=Path("/previous"))
    parser.add_argument("--outdir", type=Path, default=Path("/pilot"))
    parser.add_argument("--reference", type=Path, default=Path("/imputer"))
    parser.add_argument("--chromosome", type=int, choices=range(1, 23))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.screen, args.previous, args.outdir, args.reference)
    elif args.stage == "summarize":
        summarize(args.screen, args.outdir)
    elif args.chromosome is None:
        parser.error("--chromosome is required")
    elif args.stage == "impute":
        frozen_config(args.outdir)
        pilot.chromosome(args.screen, args.outdir, args.chromosome)
    else:
        features(args.screen, args.outdir, args.chromosome)


if __name__ == "__main__":
    main()
