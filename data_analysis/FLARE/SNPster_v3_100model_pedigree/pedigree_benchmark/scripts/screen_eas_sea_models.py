#!/usr/bin/env python3
"""Four fixed pedigree-trained decoders with population-purged diagnostics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import joblib
import numpy as np
import pandas as pd

from audit_representations import collapsed_eas_sea_tv
from screen_global_precision import global_frame
from screen_precision_retraining import common, replace_features, verify_inputs
from screen_training_decoders import assemble_training, fit_matched, unified_metrics
from screen_v1_decoders import (
    SOURCES, SEED, digest, family_donors, feature_frame, fractions, predict_decoder,
    reference_metrics, report, scored, write_json,
)
from fit_evaluate_candidate import evaluate_constraints

REPRESENTATIONS = {"diag_k30": 1, "diag_k36": 3}
MODES = tuple(f"{alias}_{kind}" for alias in REPRESENTATIONS for kind in ("ridge", "extratrees"))


def population_split(pure, training, donors, population):
    if (pure.sample_id.duplicated().any() or pure.fine.isna().any()
            or set(training.family_id) != set(donors)):
        raise ValueError("Population purge requires complete metadata and family donors")
    excluded = set(pure.loc[pure.fine.eq(population), "sample_id"])
    known = set(pure.sample_id)
    if not set().union(*donors.values()).issubset(known):
        raise ValueError("Pedigree founder missing from pure-training metadata")
    families = {family for family, ids in donors.items() if ids & excluded}
    clean_pure = pure[~pure.sample_id.isin(excluded)].copy()
    clean_training = training[~training.family_id.isin(families)].copy()
    remaining = set().union(*(ids for family, ids in donors.items() if family not in families))
    if (not len(clean_training) or set(clean_pure.source_group) != set(SOURCES)
            or remaining & excluded or clean_pure.fine.eq(population).any()):
        raise ValueError("Population purge leaves invalid training coverage or donor overlap")
    return clean_pure, clean_training, {
        "heldout_population": population, "excluded_donors": sorted(excluded),
        "excluded_families": sorted(families),
        "training_pure_donors": clean_pure.sample_id.tolist(),
        "training_founder_donors": sorted(remaining),
        "training_families": sorted(set(clean_training.family_id)),
        "training_populations": sorted(set(clean_pure.fine)),
        "n_pure_train": len(clean_pure), "n_pedigree_train": len(clean_training),
    }


def source_diagnostics(frame):
    subset = frame[frame.source_group.isin(["EAS", "SEA"])]
    if set(subset.source_group) != {"EAS", "SEA"}:
        raise ValueError("EAS/SEA diagnostics require both sources")
    return {
        "n": len(subset),
        "eas_mean_tv": float(subset.loc[subset.source_group.eq("EAS"), "TV_distance"].mean()),
        "sea_mean_tv": float(subset.loc[subset.source_group.eq("SEA"), "TV_distance"].mean()),
        "balanced_dominant_accuracy": float(subset.groupby("source_group").dominant_correct.mean().mean()),
    }


def gates(unified, real_pass, cfg):
    if (unified["n_reference"] < cfg["min_reference_samples"]
            or unified["n_pedigree"] < cfg["min_pedigree_samples"]
            or not set(cfg["required_scenarios"]).issubset(
                {r["scenario"] for r in unified["scenarios"]})):
        raise ValueError("Insufficient coverage for unchanged evaluation gates")
    result = {k: unified[k] <= v for k, v in cfg["max_gates"].items()}
    result.update({k: unified[k] >= v for k, v in cfg["min_gates"].items()})
    evaluable = unified["minor_2_5_n"] >= cfg["minor_gate_min_support"]
    if evaluable:
        result["minor_2_5_recall"] = unified["minor_2_5_recall"] >= cfg["minor_recall_min"]
    if cfg["require_real_hard_gate"]:
        result["real_world_regression"] = real_pass
    return result, evaluable


def export_features(directory, cols):
    global_path = directory / "run.global.anc.gz"
    frame = global_frame(global_path)
    if set(frame.columns) != set(cols):
        raise ValueError("Training features differ from saved representation columns")
    x = frame[cols].to_numpy(float)
    x /= x.sum(axis=1, keepdims=True)
    fractions(x, "Training features")
    pd.DataFrame(x, index=frame.index, columns=cols).to_csv(
        directory / "features.tsv", sep="\t", index_label="sample_id")
    write_json(directory / "output_audit.json", {
        "n_samples": len(frame),
        "output_sha256": {str(p): digest(p) for p in [directory / "features.tsv", global_path]},
    })
    return len(frame)


def recover_outputs(root, reuse, out, tasks, representations):
    previous = json.loads((reuse / "input_audit.json").read_text())
    tracked = [reuse / "input_audit.json", reuse / "original_screen_eas_sea_models.py",
               reuse / "original_launcher.sh"]
    for path, expected in previous["input_sha256"].items():
        actual = Path(path)
        if actual == Path(__file__):
            actual = reuse / "original_screen_eas_sea_models.py"
        elif actual == root / "run_eas_sea_model_screen.sh":
            actual = reuse / "original_launcher.sh"
        elif actual.is_relative_to("/results"):
            actual = reuse / actual.relative_to("/results")
        if digest(actual) != expected:
            raise ValueError(f"Recovery input changed: {path}")
    old_tasks = pd.read_csv(reuse / "tasks.tsv", sep="\t")
    pd.testing.assert_frame_equal(pd.DataFrame(tasks), old_tasks)
    recovered, skipped = [], []
    for row in tasks:
        old = reuse / "inference" / row["task"]
        log = old / "run.log"
        if not log.exists() or "End Time" not in log.read_text():
            skipped.append({"task": row["task"], "reason": "No completed FLARE log"})
            continue
        command = json.loads((old / "command.json").read_text())
        rep = representations[row["alias"]]
        cfg = json.loads(Path(rep["config"]).read_text())
        if command != inference_command(row["vcf"], rep, cfg, Path("/results/inference") / row["task"]):
            raise ValueError("Recovery inference parameters differ")
        new = out / "inference" / row["task"]
        new.mkdir(parents=True, exist_ok=False)
        for name in ("run.global.anc.gz", "run.log", "command.json"):
            original = old / name
            expected = digest(original)
            shutil.copyfile(original, new / name)
            if digest(new / name) != expected:
                raise ValueError("Copied inference output checksum differs")
            tracked.append(original)
        n = export_features(new, rep["feature_columns"])
        recovered.append(row["task"])
        print(f"Recovered {row['task']}: {n} descendants without rerunning FLARE", flush=True)
    write_json(out / "recovery_audit.json", {
        "source_directory": str(reuse), "recovered_tasks": recovered,
        "rerun_tasks": skipped, "all_previous_input_hashes_verified": True,
        "reused_source_sha256": {str(p): digest(p) for p in tracked},
    })
    return set(recovered), tracked


def prepare(root, source, out, reuse=None):
    data, training, cohort, paths = common(root, source)
    *_, provenance = data
    manifest = pd.read_csv(root / "config/frozen_pilot_10.tsv", sep="\t")
    representations, tasks, baselines = {}, [], []
    for alias, index in REPRESENTATIONS.items():
        directory = root / f"results_v3_representation_audit_01/rep_{index:02d}"
        audit_path = directory / "input_audit.json"
        audit = json.loads(audit_path.read_text())
        for kind in ("common_input_sha256", "representation_input_sha256"):
            for path, expected in audit[kind].items():
                if digest(path) != expected:
                    raise ValueError(f"Audited representation input changed: {path}")
            provenance["input_sha256"].update(audit[kind])
        row = manifest.iloc[index]
        if row.representation != audit["representation"]:
            raise ValueError("Audited representation identity differs from manifest")
        saved = joblib.load(Path(row.candidate_dir) / "model.joblib")
        cols = saved["feature_columns"]
        representations[alias] = {
            "representation": row.representation, "model": row.model_path,
            "config": row.model_config, "feature_columns": cols,
        }
        vcfs = sorted((root / "results_v3_pedigree_training_screen_01/cohort/pedigree/panel_genomes").glob(
            "*.autosomes.vcf.gz"))
        if len(vcfs) != 11:
            raise ValueError("Expected all eleven reusable training pedigree batches")
        for vcf in vcfs:
            batch = vcf.name.removesuffix(".autosomes.vcf.gz")
            tasks.append({"task": f"{alias}__{batch}", "alias": alias, "batch": batch, "vcf": str(vcf)})
        real_path = root / f"results_v3_sweep_100/features/dev/{row.representation}__lonni_real.features.tsv"
        feature_frame([real_path], cols)
        metrics_path = directory / "metrics.json"
        m = json.loads(metrics_path.read_text())
        baselines.append({"experiment": f"{alias}_frozen", "mean_tv": m["pedigree"]["mean_tv"],
                          "reference_mean_tv": m["pure_validation"]["mean_tv"],
                          "eas_reference_mean_tv": m["pure_validation_by_source"]["EAS"]["mean_tv"],
                          "sea_reference_mean_tv": m["pure_validation_by_source"]["SEA"]["mean_tv"]})
        paths.extend([audit_path, real_path, metrics_path, directory / "mixture_diagnostics.tsv"])
    v1_metrics = root / "results_v3_representation_audit_01/rep_00/metrics.json"
    m = json.loads(v1_metrics.read_text())
    baselines.append({"experiment": "v1_frozen", "mean_tv": m["pedigree"]["mean_tv"],
                      "reference_mean_tv": m["pure_validation"]["mean_tv"],
                      "eas_reference_mean_tv": m["pure_validation_by_source"]["EAS"]["mean_tv"],
                      "sea_reference_mean_tv": m["pure_validation_by_source"]["SEA"]["mean_tv"]})
    pd.DataFrame(baselines).to_csv(out / "frozen_baselines.tsv", sep="\t", index=False)
    pd.DataFrame(tasks).to_csv(out / "tasks.tsv", sep="\t", index=False)
    write_json(out / "representations.json", representations)
    recovered = set()
    if reuse is not None:
        recovered, recovery_paths = recover_outputs(root, reuse, out, tasks, representations)
        paths.extend(recovery_paths)
        paths.append(out / "recovery_audit.json")
    pd.DataFrame([row for row in tasks if row["task"] not in recovered],
                 columns=["task", "alias", "batch", "vcf"]).to_csv(
                     out / "pending_tasks.tsv", sep="\t", index=False)
    jar = Path("/opt/flare/flare.jar")
    evidence = root / "results_v3_representation_audit_01/collisions/global_precision_evidence.json"
    if digest(jar) != json.loads(evidence.read_text())["flare_jar_sha256"]:
        raise ValueError("Installed FLARE jar differs from audited stock jar")
    paths.extend([
        *(Path(row["vcf"]) for row in tasks), jar, evidence, v1_metrics,
        out / "tasks.tsv", out / "pending_tasks.tsv", out / "representations.json", out / "frozen_baselines.tsv",
        source / "source_reference.autosomes.vcf.gz",
        Path("/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map"),
        root / "config/unified_selection.json", root / "run_eas_sea_model_screen.sh",
        *sorted(Path(__file__).parent.glob("*.py")),
    ])
    provenance["input_sha256"].update({str(p): digest(p) for p in paths})
    provenance.update({
        "training_n_descendants": len(training), "training_n_families": cohort["n_families"],
        "training_weight": {"pure": .75, "pedigree": .25},
        "uses_old_mosaic_training_examples": False, "feature_precision": "legacy_three_decimal",
        "evaluation": "fixed donor-disjoint development pilot",
        "reference_evaluation": "validation donors disjoint from pure/pedigree decoder training",
        "population_holdout": "Decoder training only; frozen FLARE representation and reference not rebuilt.",
    })
    write_json(out / "input_audit.json", provenance)
    print("Prepared 22 training-only inference tasks, four fixed decoders; existing evaluation features reused",
          flush=True)
    print(f"FLARE tasks pending: {len(tasks) - len(recovered)}; recovered: {len(recovered)}", flush=True)


def inference_command(vcf, rep, cfg, directory):
    return [
        "java", "-Xmx32g", "-jar", "/opt/flare/flare.jar", f"gt={vcf}",
        "ref=/data/flare2_test/source_model/source_reference.autosomes.vcf.gz",
        "ref-panel=/data/flare2_test/source_model/source_reference.ref.panel",
        "map=/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map",
        f"model={rep['model']}", "update-p=false", "em=false", "probs=false", "array=false",
        f"min-maf={cfg['min_maf']}", f"min-mac={cfg['min_mac']}", "nthreads=6",
        f"seed={cfg['flare_seed']}", f"out={directory / 'run'}",
    ]


def infer(root, out, task):
    verify_inputs(out)
    tasks = pd.read_csv(out / "pending_tasks.tsv", sep="\t")
    selected = tasks[tasks.task.eq(task)]
    if len(selected) != 1:
        raise ValueError("Unknown/nonunique pending training inference task")
    row = selected.iloc[0]
    rep = json.loads((out / "representations.json").read_text())[row.alias]
    cfg = json.loads(Path(rep["config"]).read_text())
    directory = out / "inference" / task
    directory.mkdir(parents=True, exist_ok=False)
    command = inference_command(row.vcf, rep, cfg, directory)
    write_json(directory / "command.json", command)
    subprocess.run(command, check=True)
    n = export_features(directory, rep["feature_columns"])
    print(f"{task}: training features complete ({n} descendants)", flush=True)


def load_representation(root, source, out, alias):
    data, training, *_ = common(root, source)
    _, _, evaluation, pure, reference, real, *_ = data
    rep = json.loads((out / "representations.json").read_text())[alias]
    cols, representation = rep["feature_columns"], rep["representation"]
    historical = root / "results_v3_sweep_100/features/dev"
    pure = replace_features(pure, feature_frame([historical / f"{representation}__source_train.features.tsv"], cols), cols)
    reference = replace_features(reference, feature_frame([historical / f"{representation}__source_validation.features.tsv"], cols), cols)
    real = replace_features(real, feature_frame([historical / f"{representation}__lonni_real.features.tsv"], cols), cols)
    paths = sorted((root / "results_v3_pedigree_pilot_10/pedigree/features").glob(
        f"{representation}__*.features.tsv"))
    evaluation = replace_features(evaluation, feature_frame(paths, cols), cols)
    tasks = pd.read_csv(out / "tasks.tsv", sep="\t")
    paths = []
    for task in tasks.loc[tasks.alias.eq(alias), "task"]:
        directory = out / "inference" / task
        audit = json.loads((directory / "output_audit.json").read_text())
        for path, expected in audit["output_sha256"].items():
            if digest(path) != expected:
                raise ValueError(f"New training feature changed: {path}")
        paths.append(directory / "features.tsv")
    if len(paths) != 11:
        raise ValueError("Incomplete training representation inference")
    training = replace_features(training, feature_frame(paths, cols), cols)
    plan = pd.read_csv(root / "results_v3_pedigree_training_screen_01/cohort/pedigree/plans/plans/expected_ancestry.tsv",
                       sep="\t")
    return cols, pure, training, evaluation, reference, real, family_donors(plan)


def population_holdouts(kind, cols, pure, training, reference, donors, eligible, out):
    subset = reference[reference.source_group.isin(["EAS", "SEA"])]
    frames, folds = [], []
    directory = out / "population_fold_models"
    directory.mkdir()
    for index, population in enumerate(sorted(set(subset.fine))):
        clean_pure, clean_training, fold = population_split(pure, training, donors, population)
        x, y, weights = assemble_training(clean_pure, clean_training, cols, True)
        model = fit_matched(kind, x, y, weights)
        validation = subset[subset.fine.eq(population)]
        if set(validation.sample_id) & (set(clean_pure.sample_id) | set(fold["training_founder_donors"])):
            raise ValueError("Population-held-out validation overlaps training donors")
        path = directory / f"fold_{index:02d}.joblib"
        joblib.dump({"decoder": model, "feature_columns": cols, "sources": SOURCES,
                     "fold": fold, "development_screen_only": True}, path)
        pred = predict_decoder(model, validation[cols].to_numpy(float))
        np.testing.assert_allclose(pred, predict_decoder(joblib.load(path)["decoder"],
                                                        validation[cols].to_numpy(float)),
                                   atol=1e-12, rtol=0)
        frame = reference_metrics(validation, pred)
        frame["fine"] = validation.fine.to_numpy()
        frame["pedigree_eligible"] = frame.sample_id.isin(eligible)
        frame["fold"] = index
        frames.append(frame)
        folds.append({**fold, "validation_donors": validation.sample_id.tolist()})
    result = pd.concat(frames, ignore_index=True)
    if result.sample_id.duplicated().any() or set(result.sample_id) != set(subset.sample_id):
        raise ValueError("Population-held-out predictions have incomplete coverage")
    result.to_csv(out / "population_holdout_predictions.tsv", sep="\t", index=False)
    write_json(out / "population_folds.json", folds)
    return {
        "all_eas_sea": source_diagnostics(result),
        "pedigree_eligible": source_diagnostics(result[result.pedigree_eligible]),
        "n_folds": len(folds),
        "scope": "Decoder-only holdout; representation/reference are fixed, not unseen-population training.",
    }


def experiment(root, source, out, mode):
    verify_inputs(out)
    alias, kind = mode.rsplit("_", 1)
    cols, pure, training, evaluation, reference, real, donors = load_representation(root, source, out, alias)
    x, y, weights = assemble_training(pure, training, cols, True)
    model = fit_matched(kind, x, y, weights)
    directory = out / mode
    directory.mkdir(exist_ok=False)
    path = directory / "development_model.joblib"
    joblib.dump({
        "decoder": model, "calibrator": None, "feature_columns": cols, "sources": SOURCES,
        "experiment": mode, "development_screen_only": True,
        "uses_old_mosaic_training_examples": False, "input_audit_sha256": digest(out / "input_audit.json"),
    }, path)
    decoded = {}
    for dataset, frame in (("development", evaluation), ("reference", reference), ("real", real)):
        pred = predict_decoder(model, frame[cols].to_numpy(float))
        np.testing.assert_allclose(pred, predict_decoder(joblib.load(path)["decoder"],
                                                        frame[cols].to_numpy(float)), atol=1e-12, rtol=0)
        decoded[dataset] = pred
    ped = scored(evaluation, decoded["development"])
    ped["eas_sea_collapsed_tv"] = collapsed_eas_sea_tv(
        decoded["development"], evaluation[[f"realized_{s}" for s in SOURCES]].to_numpy(float))
    ped.to_csv(directory / "development_predictions.tsv", sep="\t", index=False)
    ref = reference_metrics(reference, decoded["reference"])
    ref["fine"] = reference.fine.to_numpy()
    eligible = set(pd.read_csv(root / "results_v3_pedigree_pilot_10/pedigree/inputs/eligible_heldout_founders.tsv",
                              sep="\t")["sample"])
    ref["pedigree_eligible"] = ref.sample_id.isin(eligible)
    ref.to_csv(directory / "reference_predictions.tsv", sep="\t", index=False)
    population_metrics = population_holdouts(kind, cols, pure, training, reference, donors, eligible, directory)
    rules = pd.read_csv(root / "config/real_benchmarks.tsv", sep="\t")
    rules = rules[rules.dataset_id.eq("lonni_real")]
    if len(rules) != 1 or len(real) != 1:
        raise ValueError("Expected one known real-user regression case")
    prediction = dict(zip(SOURCES, decoded["real"][0]))
    passed, details = evaluate_constraints(prediction, json.loads(rules.iloc[0].constraints_json))
    cfg = json.loads((root / "config/unified_selection.json").read_text())
    unified = unified_metrics(ped, ref, cfg)
    statuses, minor_evaluable = gates(unified, bool(passed), cfg)
    write_json(directory / "metrics.json", {
        "experiment": mode, "pedigree": report(ped), "reference": report(ref),
        "pure_by_source": {str(s): report(g) for s, g in ref.groupby("source_group")},
        "eas_sea": source_diagnostics(ref),
        "eligible_eas_sea": source_diagnostics(ref[ref.pedigree_eligible]),
        "population_holdout": population_metrics, "unified": unified,
        "real_regression": {"pass": bool(passed), "prediction": prediction, "details": details},
        "gate_results": statuses, "minor_gate_evaluable": minor_evaluable,
        "mean_eas_sea_collapsed_pedigree_tv": float(ped.eas_sea_collapsed_tv.mean()),
        "input_audit_sha256": digest(out / "input_audit.json"),
        "parameters": model.get_params(), "seed": SEED,
        "promoted": False, "development_screen_only": True,
    })
    print(f"{mode}: complete; pedigree mean TV={report(ped)['mean_tv']:.6f}", flush=True)


def summarize(out):
    verify_inputs(out)
    rows = []
    for mode in MODES:
        m = json.loads((out / mode / "metrics.json").read_text())
        if m["input_audit_sha256"] != digest(out / "input_audit.json"):
            raise ValueError("Input versions differ between models")
        rows.append({
            "experiment": mode, "mean_tv": m["pedigree"]["mean_tv"],
            "p95_tv": m["pedigree"]["p95_tv"], "reference_mean_tv": m["reference"]["mean_tv"],
            "eas_reference_mean_tv": m["eas_sea"]["eas_mean_tv"],
            "sea_reference_mean_tv": m["eas_sea"]["sea_mean_tv"],
            "eas_sea_balanced_accuracy": m["eas_sea"]["balanced_dominant_accuracy"],
            "population_holdout_eas_tv": m["population_holdout"]["all_eas_sea"]["eas_mean_tv"],
            "population_holdout_sea_tv": m["population_holdout"]["all_eas_sea"]["sea_mean_tv"],
            "population_holdout_balanced_accuracy": m["population_holdout"]["all_eas_sea"]["balanced_dominant_accuracy"],
            "eligible_population_holdout_balanced_accuracy": m["population_holdout"]["pedigree_eligible"]["balanced_dominant_accuracy"],
            "balanced_mean_tv": m["unified"]["balanced_mean_tv"],
            "balanced_p95_tv": m["unified"]["balanced_p95_tv"],
            "balanced_off_target_mass": m["unified"]["balanced_off_target_mass"],
            "real_regression_pass": m["real_regression"]["pass"],
            "minor_gate_evaluable": m["minor_gate_evaluable"],
            "passes_evaluable_gates": all(m["gate_results"].values()),
            "gate_failures": ",".join(k for k, status in m["gate_results"].items() if not status),
            "development_screen_only": True, "promoted": False,
        })
    for kind in ("ridge", "extratrees"):
        configs = [json.loads((out / f"{alias}_{kind}/metrics.json").read_text()) for alias in REPRESENTATIONS]
        if configs[0]["parameters"] != configs[1]["parameters"] or configs[0]["seed"] != configs[1]["seed"]:
            raise ValueError("Fixed learner configurations differ between representations")
    pd.DataFrame(rows).to_csv(out / "eas_sea_model_screen.tsv", sep="\t", index=False)
    write_json(out / "completion.json", {
        "training_inference_tasks": 22, "main_models": 4,
        "population_fold_models": sum(json.loads((out / mode / "metrics.json").read_text())[
            "population_holdout"]["n_folds"] for mode in MODES),
        "input_hashes_verified": True, "saved_predictions_reproduced": True,
        "promoted": False, "independent_final_validation": False,
    })
    print("Four-model EAS/SEA screen complete; no model promoted", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("prepare", "infer", "summarize", *MODES))
    parser.add_argument("--task")
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reuse-dir", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.root, args.source, args.outdir, args.reuse_dir)
    elif args.stage == "infer":
        infer(args.root, args.outdir, args.task)
    elif args.stage == "summarize":
        summarize(args.outdir)
    else:
        experiment(args.root, args.source, args.outdir, args.stage)


if __name__ == "__main__":
    main()
