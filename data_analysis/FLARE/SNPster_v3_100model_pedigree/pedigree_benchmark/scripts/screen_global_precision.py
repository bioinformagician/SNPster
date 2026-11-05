#!/usr/bin/env python3
"""Matched export-only precision control; never modifies the installed FLARE jar."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

import joblib
import numpy as np
import pandas as pd

from screen_v1_decoders import SOURCES, digest, fractions, frozen_prediction, write_json

SOURCE_BLOB = "6ff16b37b23ccbd25f42b82cc503c75f84130183"
CLASS = "admix/EstimatedGlobalAncProportions.class"
REP = "v1_frozen_k18"
ARMS = ("stock", "recompiled", "precision")
VALIDATION = (
    "HGDP00720", "HG02087", "HG01597", "HG02142", "HG02178",
    "NA18605", "NA07031", "HGDP00662", "HGDP00545",
)
TRAINING = ("GRC11052377", "GRC11053498")


def precision_source(source):
    replacements = (
        ("String[] df3 = df3();", "String[] df3 = null;"),
        ("float sum = 0f;", "double sum = 0.0;"),
        ("float factor = 1000f/sum;", "double factor = 1.0/sum;"),
        ("sb.append(df3[ (int) Math.rint(factor*fa.get(j)) ]);",
         "sb.append(Double.toString(factor*fa.get(j)));"),
    )
    for old, new in replacements:
        if source.count(old) != 1:
            raise ValueError(f"Verified FLARE exporter differs: {old}")
        source = source.replace(old, new)
    return source


def local_hash(path):
    h = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for line in stream:
            if not line.startswith(b"##"):
                h.update(line)
    return h.hexdigest()


def global_frame(path):
    frame = pd.read_csv(path, sep=r"\s+")
    cols = [c for c in frame if c.startswith("anc_")]
    if not cols or frame.SAMPLE.duplicated().any():
        raise ValueError("Invalid FLARE global schema")
    frame = frame.rename(columns={"SAMPLE": "sample_id"}).set_index("sample_id")
    x = frame[cols].to_numpy(float)
    if not np.isfinite(x).all() or (x < 0).any() or (x.sum(axis=1) <= 0).any():
        raise ValueError("Invalid global proportions")
    return frame[cols]


def prepare(root, source, out):
    original = out / "input/EstimatedGlobalAncProportions.java"
    data = original.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if blob != SOURCE_BLOB:
        raise ValueError("Downloaded source is not the verified FLARE source blob")
    jar = Path("/opt/flare/flare.jar")
    evidence = root / "results_v3_representation_audit_01/collisions/global_precision_evidence.json"
    if digest(jar) != json.loads(evidence.read_text())["flare_jar_sha256"]:
        raise ValueError("Installed jar changed since collision audit")
    implementation = out / "input/inference_screen.py"
    implementation.write_bytes(Path(__file__).read_bytes())
    tracked = [original, jar, evidence, implementation]
    for arm in ("recompiled", "precision"):
        build = out / f"build/{arm}"
        build.mkdir(parents=True, exist_ok=False)
        java = build / "EstimatedGlobalAncProportions.java"
        java.write_text(data.decode() if arm == "recompiled" else precision_source(data.decode()))
        subprocess.run(["/jdk/bin/javac", "--release", "17", "-cp", str(jar),
                        "-d", str(build), str(java)], check=True)
        with zipfile.ZipFile(jar) as src, zipfile.ZipFile(out / f"{arm}.jar", "w") as dst:
            for info in src.infolist():
                dst.writestr(info, (build / CLASS).read_bytes() if info.filename == CLASS
                             else src.read(info.filename))
        with zipfile.ZipFile(jar) as src, zipfile.ZipFile(out / f"{arm}.jar") as dst:
            changed = [name for name in src.namelist() if src.read(name) != dst.read(name)]
        if changed != [CLASS]:
            raise ValueError(f"Unexpected jar changes: {changed}")
        tracked.extend([java, build / CLASS, out / f"{arm}.jar"])
    harness = out / "input/ExportControl.java"
    harness.write_text(
        'import admix.EstimatedGlobalAncProportions;\n'
        'import java.io.PrintWriter;\n'
        'public class ExportControl {\n'
        ' public static void main(String[] args) {\n'
        '  var p = new EstimatedGlobalAncProportions(new String[]{"control"}, '
        'new String[]{"anc_0", "anc_1"});\n'
        '  p.add(0, 1, new double[]{0.9998, 0.0002});\n'
        '  p.add(1, 1, new double[]{0.9996, 0.0004});\n'
        '  p.writeGlobalAncestry(new PrintWriter(System.out, true));\n'
        ' }\n'
        '}\n'
    )
    subprocess.run(["/jdk/bin/javac", "--release", "17", "-cp", str(jar),
                    "-d", str(out / "input"), str(harness)], check=True)
    for arm in ARMS:
        arm_jar = jar if arm == "stock" else out / f"{arm}.jar"
        result = subprocess.run(["java", "-cp", f"{out / 'input'}:{arm_jar}", "ExportControl"],
                                check=True, capture_output=True, text=True)
        path = out / f"input/{arm}.control.tsv"
        path.write_text(result.stdout)
        actual = global_frame(path).to_numpy()
        np.testing.assert_allclose(actual, [[.9997, .0003]] if arm == "precision" else [[1., 0.]],
                                   atol=1e-15, rtol=0)
    ref = source / "source_reference.autosomes.vcf.gz"
    reference_ids = set(subprocess.run(["bcftools", "query", "-l", str(ref)],
                                      check=True, capture_output=True, text=True).stdout.splitlines())
    rows = []
    subsets = []
    for role, ids in (("train", TRAINING), ("validation", VALIDATION)):
        vcf = source / f"source_{role}.autosomes.vcf.gz"
        metadata = source / f"source_{role}_metadata.tsv"
        meta = pd.read_csv(metadata, sep="\t").set_index("sample").loc[list(ids)].reset_index()
        if reference_ids & set(ids):
            raise ValueError("Precision targets overlap the FLARE reference")
        meta["dataset"] = f"source_{role}"
        rows.append(meta.rename(columns={"sample": "sample_id"}))
        subset = out / f"input/{role}.vcf.gz"
        subprocess.run(["bcftools", "view", "-s", ",".join(ids), "-Oz", "-o", str(subset),
                        str(vcf)], check=True)
        subprocess.run(["bcftools", "index", "-t", str(subset)], check=True)
        subsets.append(subset)
        tracked.extend([vcf, metadata])
    target = out / "input/targets.vcf.gz"
    subprocess.run(["bcftools", "merge", "-m", "none", "-Oz", "-o", str(target),
                    *map(str, subsets)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(target)], check=True)
    meta = pd.concat(rows, ignore_index=True)
    sample_ids = subprocess.run(["bcftools", "query", "-l", str(target)], check=True,
                               capture_output=True, text=True).stdout.splitlines()
    if len(sample_ids) != len(meta) or set(sample_ids) != set(meta.sample_id):
        raise ValueError("Target coverage changed during VCF merge")
    meta.to_csv(out / "samples.tsv", sep="\t", index=False)
    tracked.extend([
        target, out / "samples.tsv", ref, source / "source_reference.ref.panel",
        Path("/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map"),
        root / f"results_v3_sweep_100/models/{REP}.final.model",
        root / f"results_v3_sweep_100/models/{REP}.config.json",
        root / f"results_v3_sweep_100/stage3_full/candidate_{REP}/model.joblib",
    ])
    write_json(out / "input_audit.json", {
        "source_blob": blob, "changed_jar_entries": [CLASS], "reference_overlap": 0,
        "input_sha256": {str(p): digest(p) for p in tracked},
        "purpose": "Export-only precision diagnostic; no model training or promotion.",
    })
    print("Prepared 11 reference-disjoint targets; stock/precision exporter controls passed", flush=True)


def infer(root, out, arm):
    provenance = json.loads((out / "input_audit.json").read_text())
    jar = Path("/opt/flare/flare.jar") if arm == "stock" else out / f"{arm}.jar"
    for path in (jar, out / "input/targets.vcf.gz"):
        if digest(path) != provenance["input_sha256"][str(path)]:
            raise ValueError(f"Matched input changed: {path}")
    directory = out / arm
    directory.mkdir(exist_ok=False)
    config = json.loads((root / f"results_v3_sweep_100/models/{REP}.config.json").read_text())
    command = [
        "java", "-Xmx10g", "-jar", str(jar), f"gt={out / 'input/targets.vcf.gz'}",
        "ref=/data/flare2_test/source_model/source_reference.autosomes.vcf.gz",
        "ref-panel=/data/flare2_test/source_model/source_reference.ref.panel",
        "map=/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map",
        f"model={root / f'results_v3_sweep_100/models/{REP}.final.model'}",
        "update-p=false", "em=false", "probs=false", "array=false",
        f"min-maf={config['min_maf']}", f"min-mac={config['min_mac']}",
        "nthreads=2", f"seed={config['flare_seed']}", f"out={directory / 'run'}",
    ]
    write_json(directory / "command.json", command)
    subprocess.run(command, check=True)


def verify_inputs(provenance, recorded_script=None):
    replaced = False
    for path, expected in provenance["input_sha256"].items():
        actual_path = Path(path)
        if recorded_script is not None and actual_path == Path(__file__):
            actual_path = recorded_script
            replaced = True
        if digest(actual_path) != expected:
            raise ValueError(f"Input changed: {path}")
    if recorded_script is not None and not replaced:
        raise ValueError("Recorded-script override does not match a provenance entry")


def summarize(root, out, recorded_script=None):
    provenance = json.loads((out / "input_audit.json").read_text())
    verify_inputs(provenance, recorded_script)
    local = {arm: local_hash(out / arm / "run.anc.vcf.gz") for arm in ARMS}
    if len(set(local.values())) != 1:
        raise ValueError("Local ancestry differs: export-only control not established")
    frames = {arm: global_frame(out / arm / "run.global.anc.gz") for arm in ARMS}
    meta = pd.read_csv(out / "samples.tsv", sep="\t").set_index("sample_id")
    for frame in frames.values():
        if set(frame.index) != set(meta.index):
            raise ValueError("Inference sample coverage differs")
    frames = {arm: frame.loc[meta.index] for arm, frame in frames.items()}
    pd.testing.assert_frame_equal(frames["stock"], frames["recompiled"])
    np.testing.assert_allclose(frames["precision"].sum(axis=1), 1., atol=1e-12, rtol=0)
    np.testing.assert_allclose(np.rint(frames["precision"].to_numpy() * 1000) / 1000,
                               frames["stock"].to_numpy(), atol=1e-12, rtol=0)
    raw_difference = np.abs(frames["precision"].to_numpy() - frames["stock"].to_numpy())
    if raw_difference.max() > .000501:
        raise ValueError("Precision output exceeds expected three-decimal rounding interval")
    saved = joblib.load(root / f"results_v3_sweep_100/stage3_full/candidate_{REP}/model.joblib")
    saved["decoder"].set_params(n_jobs=2)
    cols = saved["feature_columns"]
    features, predictions = {}, {}
    results = meta.copy()
    truth = np.column_stack([meta.source_group.eq(s) for s in SOURCES]).astype(float)
    for arm, frame in frames.items():
        x = frame[cols].to_numpy(float)
        x /= x.sum(axis=1, keepdims=True)
        features[arm] = fractions(x, arm)
        pred = frozen_prediction(saved, x)
        predictions[arm] = pred
        results[f"{arm}_TV"] = np.abs(pred - truth).sum(axis=1) / 2
        for j, source in enumerate(SOURCES):
            results[f"{arm}_pred_{source}"] = pred[:, j]
        pd.DataFrame(x, index=meta.index, columns=cols).to_csv(
            out / f"{arm}.features.tsv", sep="\t", index_label="sample_id")
    np.testing.assert_allclose(predictions["stock"], predictions["recompiled"],
                               atol=1e-12, rtol=0)
    historical_differences = {}
    historical_hashes = {}
    for dataset, subset in meta.groupby("dataset"):
        path = root / f"results_v3_sweep_100/features/dev/{REP}__{dataset}.features.tsv"
        historical_hashes[str(path)] = digest(path)
        previous = pd.read_csv(path, sep="\t").set_index("sample_id").loc[subset.index, cols]
        ids = meta.index.get_indexer(subset.index)
        historical_differences[dataset] = float(
            np.abs(previous.to_numpy(float) - features["stock"][ids]).max())
    results["normalized_feature_TV_change"] = (
        np.abs(features["stock"] - features["precision"]).sum(axis=1) / 2)
    results["frozen_prediction_TV_change"] = (
        np.abs(predictions["stock"] - predictions["precision"]).sum(axis=1) / 2)
    pairs = []
    for donor in TRAINING:
        i, j = meta.index.get_loc("HGDP00720"), meta.index.get_loc(donor)
        pairs.append({
            "validation_donor": "HGDP00720", "training_donor": donor,
            **{f"{arm}_feature_TV": float(np.abs(x[i] - x[j]).sum() / 2)
               for arm, x in features.items()},
            **{f"{arm}_exact_collision": bool(np.array_equal(x[i], x[j]))
               for arm, x in features.items()},
        })
    results.to_csv(out / "sample_predictions.tsv", sep="\t", index_label="sample_id")
    write_json(out / "summary.json", {
        "n_samples": len(meta), "local_ancestry_hashes": local,
        "stock_matches_recompiled": True, "input_hashes_verified": True,
        "rounding_precision_vectors_reproduces_stock": True,
        "input_audit_sha256": digest(out / "input_audit.json"),
        "summary_implementation_sha256": digest(Path(__file__)),
        "archived_inference_script": str(recorded_script) if recorded_script else None,
        "max_raw_export_difference": float(raw_difference.max()),
        "known_collision_pairs": pairs,
        "mean_pure_TV": {arm: float(results[f"{arm}_TV"].mean()) for arm in ARMS},
        "mean_pure_TV_by_dataset": {
            dataset: {"n_samples": len(subset),
                      **{arm: float(subset[f"{arm}_TV"].mean()) for arm in ARMS}}
            for dataset, subset in results.groupby("dataset")
        },
        "max_normalized_feature_TV_change": float(results.normalized_feature_TV_change.max()),
        "max_frozen_prediction_TV_change": float(results.frozen_prediction_TV_change.max()),
        "stock_vs_historical_feature_max_abs_difference": historical_differences,
        "historical_feature_sha256": historical_hashes,
        "limitations": "Selected pure controls, including training donors; not a leaderboard. "
                        "Frozen decoder was trained on rounded features; no precision-aware retraining.",
    })
    print((out / "summary.json").read_text(), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", *ARMS, "summarize"))
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--recorded-script", type=Path,
                        help="Explicit archived implementation for pre-snapshot input audits")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.root, args.source, args.outdir)
    elif args.stage == "summarize":
        summarize(args.root, args.outdir, args.recorded_script)
    else:
        infer(args.root, args.outdir, args.stage)


if __name__ == "__main__":
    main()
