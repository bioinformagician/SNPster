#!/usr/bin/env python3
"""Matched rounded versus precise feature training on audited, fixed cohorts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

import joblib
import numpy as np
import pandas as pd

from screen_global_precision import global_frame, REP
from screen_training_decoders import (
    assemble_training, fit_matched, load_training, unified_metrics,
)
from screen_v1_decoders import (
    SOURCES, SEED, digest, feature_frame, fractions, load_inputs,
    predict_decoder, reference_metrics, report, scored, write_json,
)
from fit_evaluate_candidate import evaluate_constraints

KINDS = ("ridge", "stochastic", "extratrees")
ARMS = ("rounded", "precision")
MODES = tuple(f"{kind}_{arm}" for kind in KINDS for arm in ARMS)


def paired_features(frame):
    precise = frame.to_numpy(float)
    if (not np.isfinite(precise).all() or (precise < 0).any()
            or not np.allclose(precise.sum(axis=1), 1., atol=1e-12, rtol=0)):
        raise ValueError("Precise globals must be finite nonnegative unit-sum vectors")
    rounded = np.rint(1000 * precise) / 1000
    if np.any(rounded.sum(axis=1) <= 0):
        raise ValueError("Quantization produced a zero-sum feature row")
    return {"rounded": fractions(rounded / rounded.sum(axis=1, keepdims=True), "rounded"),
            "precision": fractions(precise / precise.sum(axis=1, keepdims=True), "precision")}


def replace_features(frame, features, cols):
    if (frame.sample_id.duplicated().any() or features.sample_id.duplicated().any()
            or set(frame.sample_id) != set(features.sample_id)):
        raise ValueError("Feature replacement requires exact unique sample coverage")
    result = frame.copy()
    result[cols] = features.set_index("sample_id").loc[frame.sample_id, cols].to_numpy(float)
    fractions(result[cols], "replacement")
    return result


def verify_inputs(out):
    audit = json.loads((out / "input_audit.json").read_text())
    for path, expected in audit["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Retraining input changed: {path}")
    return audit


def common(root, source):
    data = load_inputs(root, source)
    saved, cols, evaluation, pure, reference, real, donors, *_ = data
    training_dir = root / "results_v3_pedigree_training_screen_01/cohort"
    training, cohort_audit, paths = load_training(
        training_dir, cols, pure, evaluation, set().union(*donors.values()), root)
    return data, training, cohort_audit, paths


def prepare(root, source, out):
    data, training, cohort_audit, paths = common(root, source)
    _, cols, evaluation, pure, reference, real, *_, provenance = data
    precision_dir = root / "results_v3_global_precision_01"
    summary_path = precision_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    if (summary["stock_matches_recompiled"] is not True
            or summary["rounding_precision_vectors_reproduces_stock"] is not True
            or len(set(summary["local_ancestry_hashes"].values())) != 1):
        raise ValueError("Verified export-only control is required")
    old_audit = json.loads((precision_dir / "input_audit.json").read_text())
    jar = precision_dir / "precision.jar"
    if digest(jar) != old_audit["input_sha256"][str(jar).replace(str(precision_dir), "/results")]:
        raise ValueError("Precision jar differs from verified experiment")
    tasks = []
    historical = root / "results_v3_sweep_100/features/dev"
    for dataset in ("source_train", "source_validation", "lonni_real"):
        vcf = (source / f"{dataset}.autosomes.vcf.gz" if dataset != "lonni_real" else
               source / "real_user_validation/lonni/lonni.autosomes.vcf.gz")
        tasks.append({"task": dataset, "dataset": dataset, "vcf": str(vcf),
                      "historical_features": str(historical / f"{REP}__{dataset}.features.tsv")})
    for dataset, directory in (
        ("training", root / "results_v3_pedigree_training_screen_01/cohort"),
        ("evaluation", root / "results_v3_pedigree_pilot_10"),
    ):
        vcfs = sorted((directory / "pedigree/panel_genomes").glob("*.autosomes.vcf.gz"))
        old = sorted((directory / "pedigree/features").glob(f"{REP}__*.features.tsv"))
        if not vcfs or len(vcfs) != len(old):
            raise ValueError(f"{dataset}: incomplete reusable inference inputs")
        for vcf in vcfs:
            batch = vcf.name.removesuffix(".autosomes.vcf.gz")
            tasks.append({"task": f"{dataset}_{batch}", "dataset": dataset, "vcf": str(vcf),
                          "historical_features": str(directory / f"pedigree/features/{REP}__{batch}.features.tsv")})
    frame = pd.DataFrame(tasks)
    if frame.task.duplicated().any() or len(frame) != 17:
        raise ValueError("Expected exactly 17 distinct inference tasks")
    frame.to_csv(out / "tasks.tsv", sep="\t", index=False)
    paths.extend(Path(row[field]) for row in tasks for field in ("vcf", "historical_features"))
    paths.extend([
        jar, summary_path, precision_dir / "input_audit.json", out / "tasks.tsv",
        source / "source_reference.autosomes.vcf.gz", source / "source_reference.ref.panel",
        Path("/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map"),
        root / "config/unified_selection.json",
        root / "run_precision_retraining.sh",
        *sorted(Path(__file__).parent.glob("*.py")),
    ])
    provenance["input_sha256"].update({str(p): digest(p) for p in paths})
    provenance.update({
        "precision_jar": str(jar), "feature_columns": cols,
        "training_n_descendants": len(training), "training_n_families": cohort_audit["n_families"],
        "n_pure_train": len(pure), "n_reference": len(reference), "n_real": len(real),
        "n_descendants": len(evaluation), "n_families": evaluation.family_id.nunique(),
        "evaluation": "fixed training pool versus donor-disjoint, previously inspected development pilot",
        "reference_evaluation": "all validation donors disjoint from decoder training donors",
        "feature_precision_control": "Same inferred posteriors; three-decimal np.rint then normalization versus double export.",
        "uses_old_mosaic_training_examples": False,
        "training_weight": {"pure": .75, "pedigree": .25},
    })
    write_json(out / "input_audit.json", provenance)
    print("Prepared 17 inference tasks; fixed 307 pure + 492 pedigree training rows, "
          "159 pure validation + 123 development descendants + one real regression", flush=True)


def infer(root, out, task):
    audit = verify_inputs(out)
    tasks = pd.read_csv(out / "tasks.tsv", sep="\t")
    selected = tasks[tasks.task.eq(task)]
    if len(selected) != 1:
        raise ValueError(f"Unknown/nonunique inference task: {task}")
    row = selected.iloc[0]
    directory = out / "inference" / task
    directory.mkdir(parents=True, exist_ok=False)
    cfg = json.loads((root / f"results_v3_sweep_100/models/{REP}.config.json").read_text())
    command = [
        "java", "-Xmx32g", "-jar", audit["precision_jar"], f"gt={row.vcf}",
        "ref=/data/flare2_test/source_model/source_reference.autosomes.vcf.gz",
        "ref-panel=/data/flare2_test/source_model/source_reference.ref.panel",
        "map=/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map",
        f"model={root / f'results_v3_sweep_100/models/{REP}.final.model'}",
        "update-p=false", "em=false", "probs=false", "array=false",
        f"min-maf={cfg['min_maf']}", f"min-mac={cfg['min_mac']}",
        "nthreads=6", f"seed={cfg['flare_seed']}", f"out={directory / 'run'}",
    ]
    write_json(directory / "command.json", command)
    subprocess.run(command, check=True)
    global_path = directory / "run.global.anc.gz"
    frame = global_frame(global_path)
    cols = audit["feature_columns"]
    if set(frame.columns) != set(cols):
        raise ValueError("Precise feature columns differ from frozen representation")
    frame = frame[cols]
    old = pd.read_csv(row.historical_features, sep="\t").set_index("sample_id")
    if old.index.has_duplicates or set(old.index) != set(frame.index):
        raise ValueError("Precise versus historical inference sample IDs differ")
    paired = paired_features(frame)
    for arm, x in paired.items():
        pd.DataFrame(x, index=frame.index, columns=cols).to_csv(
            directory / f"{arm}.features.tsv", sep="\t", index_label="sample_id")
    difference = np.abs(paired["rounded"] - old.loc[frame.index, cols].to_numpy(float))
    write_json(directory / "metrics.json", {
        "task": task, "dataset": row.dataset, "n_samples": len(frame),
        "historical_max_abs_difference": float(difference.max()),
        "historical_n_different_rows": int((difference.max(axis=1) > 1e-12).sum()),
        "matched_rounded_versus_precision": True,
        "rounding_definition": "Nearest three-decimal values of same precise posteriors, normalized.",
        "output_sha256": {str(p): digest(p) for p in
                          [global_path, directory / "rounded.features.tsv", directory / "precision.features.tsv"]},
    })
    print(f"{task}: paired features complete; historical max difference={difference.max():.8g}",
          flush=True)


def load_paired(root, source, out, arm):
    if arm not in ARMS:
        raise ValueError("Unknown precision arm")
    data, training, *_ = common(root, source)
    _, cols, evaluation, pure, reference, real, *_ = data
    frames = {}
    tasks = pd.read_csv(out / "tasks.tsv", sep="\t")
    for dataset, group in tasks.groupby("dataset"):
        paths = []
        for task in group.task:
            directory = out / "inference" / task
            audit = json.loads((directory / "metrics.json").read_text())
            for path, expected in audit["output_sha256"].items():
                if digest(path) != expected:
                    raise ValueError(f"Paired features changed: {path}")
            paths.append(directory / f"{arm}.features.tsv")
        frames[dataset] = feature_frame(paths, cols)
    ordered = [
        replace_features(frame, frames[dataset], cols)
        for dataset, frame in (("source_train", pure), ("training", training),
                               ("evaluation", evaluation), ("source_validation", reference),
                               ("lonni_real", real))
    ]
    return cols, ordered


def experiment(root, source, out, mode):
    audit = verify_inputs(out)
    kind, arm = mode.split("_")
    cols, (pure, training, evaluation, reference, real) = load_paired(root, source, out, arm)
    x, y, weights = assemble_training(pure, training, cols, True)
    model = fit_matched(kind, x, y, weights)
    directory = out / mode
    directory.mkdir(exist_ok=False)
    model_path = directory / "development_model.joblib"
    joblib.dump({
        "decoder": model, "calibrator": None, "sources": SOURCES, "feature_columns": cols,
        "experiment": mode, "feature_precision": arm, "development_screen_only": True,
        "uses_old_mosaic_training_examples": False, "input_audit_sha256": digest(out / "input_audit.json"),
    }, model_path)
    reloaded = joblib.load(model_path)["decoder"]
    scored_sets = {}
    for dataset, frame in (("development", evaluation), ("reference", reference), ("real", real)):
        pred = predict_decoder(model, frame[cols].to_numpy(float))
        np.testing.assert_allclose(pred, predict_decoder(reloaded, frame[cols].to_numpy(float)),
                                   atol=1e-12, rtol=0)
        if dataset == "real":
            rules = pd.read_csv(root / "config/real_benchmarks.tsv", sep="\t")
            rules = rules[rules.dataset_id.eq("lonni_real")]
            if len(rules) != 1 or len(frame) != 1:
                raise ValueError("Expected one known real-user regression case")
            passed, details = evaluate_constraints(dict(zip(SOURCES, pred[0])),
                                                   json.loads(rules.iloc[0].constraints_json))
            real_metrics = {"pass": bool(passed), "details": details,
                            "prediction": dict(zip(SOURCES, pred[0]))}
        else:
            scored_sets[dataset] = (scored(frame, pred) if dataset == "development"
                                     else reference_metrics(frame, pred))
            scored_sets[dataset].to_csv(directory / f"{dataset}_predictions.tsv", sep="\t", index=False)
    cfg = json.loads((root / "config/unified_selection.json").read_text())
    unified = unified_metrics(scored_sets["development"], scored_sets["reference"], cfg)
    parameters = ({"kind": model["kind"], "alpha": model["alpha"]}
                  if isinstance(model, dict) else model.get_params())
    write_json(directory / "experiment.json", {
        "experiment": mode, "parameters": parameters, "seed": SEED,
        "training_weight": audit["training_weight"], "calibrator": None,
        "uses_old_mosaic_training_examples": False, "development_screen_only": True, "promoted": False,
    })
    write_json(directory / "metrics.json", {
        "experiment": mode, "pedigree": report(scored_sets["development"]),
        "reference": report(scored_sets["reference"]), "unified": unified,
        "real_regression": real_metrics, "saved_model_predictions_reproduced": True,
        "input_audit_sha256": digest(out / "input_audit.json"),
        "promoted": False, "development_screen_only": True,
    })
    print(f"{mode}: complete, pedigree mean TV={report(scored_sets['development'])['mean_tv']:.6f}",
          flush=True)


def summarize(root, out):
    verify_inputs(out)
    tasks = pd.read_csv(out / "tasks.tsv", sep="\t")
    rows = []
    cfg = json.loads((root / "config/unified_selection.json").read_text())
    for mode in MODES:
        m = json.loads((out / mode / "metrics.json").read_text())
        if m["input_audit_sha256"] != digest(out / "input_audit.json"):
            raise ValueError("Model experiment input audits differ")
        u = m["unified"]
        if u["n_reference"] < cfg["min_reference_samples"] or u["n_pedigree"] < cfg["min_pedigree_samples"]:
            raise ValueError("Insufficient samples for existing evaluation gates")
        if not set(cfg["required_scenarios"]).issubset({r["scenario"] for r in u["scenarios"]}):
            raise ValueError("Missing required evaluation scenarios")
        gates = {k: u[k] <= v for k, v in cfg["max_gates"].items()}
        gates.update({k: u[k] >= v for k, v in cfg["min_gates"].items()})
        minor_evaluable = u["minor_2_5_n"] >= cfg["minor_gate_min_support"]
        if minor_evaluable:
            gates["minor_2_5_recall"] = u["minor_2_5_recall"] >= cfg["minor_recall_min"]
        if cfg["require_real_hard_gate"]:
            gates["real_world_regression"] = m["real_regression"]["pass"]
        rows.append({
            "experiment": mode, "mean_tv": m["pedigree"]["mean_tv"],
            "p95_tv": m["pedigree"]["p95_tv"], "reference_mean_tv": m["reference"]["mean_tv"],
            "balanced_mean_tv": m["unified"]["balanced_mean_tv"],
            "balanced_p95_tv": m["unified"]["balanced_p95_tv"],
            "balanced_off_target_mass": m["unified"]["balanced_off_target_mass"],
            "minor_2_5_n": m["pedigree"]["minor_2_5_n"],
            "minor_2_5_recall": m["pedigree"]["minor_2_5_recall"],
            "real_regression_pass": m["real_regression"]["pass"],
            "max_gates_pass": all(m["unified"][k] <= v for k, v in cfg["max_gates"].items()),
            "minor_gate_evaluable": minor_evaluable,
            "passes_evaluable_gates": all(gates.values()),
            "gate_failures": ",".join(k for k, passed in gates.items() if not passed),
            "development_screen_only": True, "promoted": False,
        })
    for kind in KINDS:
        configs = [json.loads((out / f"{kind}_{arm}/experiment.json").read_text()) for arm in ARMS]
        for field in ("parameters", "seed", "training_weight", "calibrator"):
            if configs[0][field] != configs[1][field]:
                raise ValueError(f"{kind}: matched {field} differs")
        old, new = [pd.read_csv(out / f"{kind}_{arm}/development_predictions.tsv", sep="\t")
                    for arm in ARMS]
        if set(old.sample_id) != set(new.sample_id):
            raise ValueError("Matched predictions have different sample IDs")
        delta = new.groupby("family_id").TV_distance.mean() - old.groupby("family_id").TV_distance.mean()
        delta.rename("precision_minus_rounded_mean_tv").to_csv(
            out / f"{kind}.paired_family_deltas.tsv", sep="\t")
    pd.DataFrame(rows).to_csv(out / "matched_precision_screen.tsv", sep="\t", index=False)
    metrics = [json.loads((out / "inference" / task / "metrics.json").read_text()) for task in tasks.task]
    pd.DataFrame([{k: v for k, v in metric.items() if k != "output_sha256"} for metric in metrics]).to_csv(
        out / "feature_precision_audit.tsv", sep="\t", index=False)
    write_json(out / "completion.json", {
        "inference_tasks": len(tasks), "models": len(MODES), "input_hashes_verified": True,
        "all_saved_predictions_reproduced": True, "promoted": False,
        "historical_different_rows": sum(m["historical_n_different_rows"] for m in metrics),
        "limitations": "Previously inspected development cohort; no independent final validation. "
                        "Matched rounded arm is explicit quantization of the same posterior summaries.",
    })
    print("Matched precision screen complete; no model promoted", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("prepare", "infer", "summarize", *MODES))
    parser.add_argument("--task")
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.root, args.source, args.outdir)
    elif args.stage == "infer":
        infer(args.root, args.outdir, args.task)
    elif args.stage == "summarize":
        summarize(args.root, args.outdir)
    else:
        experiment(args.root, args.source, args.outdir, args.stage)


if __name__ == "__main__":
    main()
