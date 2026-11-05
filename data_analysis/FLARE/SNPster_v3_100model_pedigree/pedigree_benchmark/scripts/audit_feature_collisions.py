#!/usr/bin/env python3
"""Measure exact recorded-feature collisions and check a genotype contrast."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct
import subprocess
import zipfile

import numpy as np
import pandas as pd

from screen_v1_decoders import digest, write_json


def collisions(frame, cols):
    if frame.sample_id.duplicated().any() or not len(frame):
        raise ValueError("Collision audit requires unique, nonempty sample IDs")
    rows, minority = [], 0
    for _, group in frame.groupby(cols, sort=False):
        counts = group.source_group.value_counts()
        minority += len(group) - counts.max()
        if len(counts) < 2:
            continue
        group_id = len(rows) + 1
        rows.append({
            "collision_group": group_id, "n_samples": len(group),
            "sources": ",".join(sorted(counts.index)),
            "source_counts": json.dumps(counts.to_dict(), sort_keys=True),
            "sample_ids": ",".join(group.sample_id),
            "includes_hgdp00720": bool(group.sample_id.eq("HGDP00720").any()),
        })
    return rows, {
        "n_samples": len(frame), "n_ambiguous_groups": len(rows),
        "n_samples_in_ambiguous_groups": sum(r["n_samples"] for r in rows),
        "exact_collision_minimum_mean_tv": float(minority / len(frame)),
        "interpretation": (
            "Finite-dataset lower bound for one deterministic decoder of these exact "
            "recorded features and pure-source labels, not a biological or validation-performance ceiling."
        ),
    }


def genotype_records(vcf, samples):
    result = subprocess.run([
        "bcftools", "query", "-r", "22", "-s", ",".join(samples),
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf),
    ], check=True, capture_output=True, text=True)
    records = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        key, gt = tuple(fields[:4]), fields[4:]
        if key in records or len(gt) != len(samples):
            raise ValueError("Duplicate/invalid genotype contrast records")
        if any(g not in ("0|0", "0|1", "1|0", "1|1") for g in gt):
            raise ValueError("Missing or unphased genotype in contrast")
        records[key] = gt
    if not records:
        raise ValueError("No genotype contrast records")
    return records


def disagreements(left, right):
    if not left or len(left) != len(right):
        raise ValueError("Genotype contrasts require aligned nonempty records")
    return sum(sorted(a.split("|")) != sorted(b.split("|")) for a, b in zip(left, right))


def run(root, source, audit_dir, out):
    summaries, all_rows, paths = [], [], []
    for index in range(10):
        directory = audit_dir / f"rep_{index:02d}"
        provenance = json.loads((directory / "input_audit.json").read_text())
        for kind in ("common_input_sha256", "representation_input_sha256"):
            for path, expected in provenance[kind].items():
                if digest(path) != expected:
                    raise ValueError(f"Representation audit input changed: {path}")
        pure_path = directory / "pure_predictions.tsv"
        pure = pd.read_csv(pure_path, sep="\t")
        rep = provenance["representation"]
        features = []
        for dataset in ("source_train", "source_validation"):
            path = root / f"results_v3_sweep_100/features/dev/{rep}__{dataset}.features.tsv"
            paths.append(path)
            features.append(pd.read_csv(path, sep="\t"))
        features = pd.concat(features, ignore_index=True)
        cols = [c for c in features if c.startswith("anc_")]
        if not cols:
            raise ValueError("Missing ancestry-state features")
        frame = features.merge(pure[["sample_id", "source_group", "dataset"]],
                               on="sample_id", how="outer", validate="one_to_one", indicator=True)
        if not frame._merge.eq("both").all():
            raise ValueError("Collision features and labels differ")
        if index == 0:
            base = frame.set_index("sample_id")
            target = base.loc["HGDP00720", cols].to_numpy(float)
            for donor in ("GRC11052377", "GRC11053498"):
                if not np.array_equal(target, base.loc[donor, cols].to_numpy(float)):
                    raise ValueError("Expected genotype-control feature collision was not reproduced")
        for scope in ("combined", "source_train", "source_validation"):
            subset = frame if scope == "combined" else frame[frame.dataset.eq(scope)]
            rows, metrics = collisions(subset, cols)
            summaries.append({"representation": rep, "scope": scope, **metrics})
            all_rows.extend({"representation": rep, "scope": scope, **row} for row in rows)
        paths.extend([directory / "input_audit.json", pure_path])
    pd.DataFrame(summaries).to_csv(out / "collision_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_rows).to_csv(out / "cross_source_collisions.tsv", sep="\t", index=False)
    val_vcf = source / "source_validation.autosomes.vcf.gz"
    train_vcf = source / "source_train.autosomes.vcf.gz"
    targets = ["GRC11052377", "GRC11053498"]
    left = genotype_records(val_vcf, ["HGDP00720"])
    right = genotype_records(train_vcf, targets)
    common = sorted(set(left) & set(right))
    if not common:
        raise ValueError("Genotype contrast marker sets do not overlap")
    checks = []
    for index, donor in enumerate(targets):
        count = disagreements([left[k][0] for k in common], [right[k][index] for k in common])
        if count == 0:
            raise ValueError("Feature-collision control has identical chr22 genotypes; investigate sample identity")
        checks.append({
            "validation_donor": "HGDP00720", "training_donor": donor,
            "chromosome": 22, "matched_sites": len(common),
            "unphased_genotype_disagreements": count,
            "unphased_concordance": 1 - count / len(common),
            "recorded_v1_features_identical": True,
        })
    write_json(out / "genotype_contrasts.json", checks)
    jar = Path("/opt/flare/flare.jar")
    with zipfile.ZipFile(jar) as archive:
        bytecode = archive.read("admix/EstimatedGlobalAncProportions.class")
        estimated = archive.read("admix/EstimatedAncestry.class")
    precision = {
        "flare_jar_sha256": digest(jar),
        "runtime_global_class_has_three_decimal_format": b"\x01\x00\x05#.###" in bytecode,
        "runtime_global_class_has_scale_1000_constant": b"\x04" + struct.pack(">f", 1000.) in bytecode,
        "runtime_global_class_has_rint_reference": b"rint" in bytecode,
        "runtime_local_class_has_global_ancestry_accumulator": b"globalAncProbs" in estimated,
        "verified_public_source": {
            "commit": "d6a6c8fd4479ebdb68529824acd11d2b6bbac4a1",
            "global_file_sha": "6ff16b37b23ccbd25f42b82cc503c75f84130183",
            "ancestry_file_sha": "67bd9cfbd77bcf0815ceeb3db6d7e349c8bf07b3",
        },
        "interpretation": (
            "Public FLARE 0.6.0 source averages posterior probabilities even with probs=false, "
            "then rounds global values to three decimals. Installed classes contain matching "
            "formatter/scale/method constants. Identical exported features do not prove identical "
            "unrounded posteriors; changing probs alone does not increase global-file precision."
        ),
    }
    if not all(precision[key] for key in (
        "runtime_global_class_has_three_decimal_format",
        "runtime_global_class_has_scale_1000_constant", "runtime_global_class_has_rint_reference",
        "runtime_local_class_has_global_ancestry_accumulator",
    )):
        raise ValueError("Installed FLARE precision evidence differs from verified source")
    write_json(out / "global_precision_evidence.json", precision)
    paths.extend([val_vcf, train_vcf, Path(__file__), jar])
    write_json(out / "input_audit.json", {
        "diagnostic_only": True, "input_sha256": {str(p): digest(p) for p in sorted(set(paths))},
    })
    print("Cross-source feature collisions and distinct-genotype controls verified", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    run(args.root, args.source, args.audit_dir, args.outdir)


if __name__ == "__main__":
    main()
