#!/usr/bin/env python3
"""Generate deliberately artificial, tiny mechanics-only inputs (not accuracy evidence)."""
import argparse
import csv
from pathlib import Path

import numpy as np

from core import write, write_vcf


def fixture(output, tools, *, models=None):
    output.mkdir()
    seed = 173
    rng = np.random.default_rng(seed)
    sources = ["B", "A"]
    keys = [("22", 501 + index * 1000, "A", "G") for index in range(256)]
    span = [0, 256000]
    samples, segments, calls = [], [], {}
    for role, size in (("reference", 24), ("train", 8), ("calibration", 4), ("validation", 4), ("test", 4)):
        for index in range(size):
            identity = f"{role}_{index}"
            source = sources[index % 2]
            mixed = role != "reference" and index % 4 >= 2
            samples.append({
                "sample_id": identity, "role": role, "family_id": f"family_{identity}", "donors": f"founder_{identity}",
                "source": source if not mixed else "", "fine": source, "track": "mixed" if mixed else "pure",
            })
            local = np.full((2, len(keys)), sources.index(source))
            if mixed:
                local[1, :] = 1 - local[0, 0]
            calls[identity] = rng.binomial(1, np.where(local == 0, .08, .92)).astype(np.int8)
            if role != "reference":
                for hap in (0, 1):
                    segments.append({"sample_id": identity, "chromosome": 22, "haplotype": hap,
                                     "start": span[0], "end": span[1], "source": sources[local[hap, 0]]})
    for name, records in (("samples", samples), ("segments", segments)):
        with (output / f"{name}.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    groups = {
        "reference": [row["sample_id"] for row in samples if row["role"] == "reference"],
        "training": [row["sample_id"] for row in samples if row["role"] == "train"],
        "query": [row["sample_id"] for row in samples if row["role"] in ("calibration", "validation", "test")],
    }
    for name, ids in groups.items():
        write_vcf(output / f"{name}.vcf.gz", keys, ids, np.stack([calls[identity] for identity in ids]))
    ids = groups["query"]
    write_vcf(output / "harmonized.vcf.gz", keys, ids, np.stack([calls[identity] for identity in ids]), phased=False)
    (output / "chr22.map").write_text("22\tleft\t0\t0\n22\tright\t2\t256000\n")
    if models is None:
        models = [
            {"id": "flare", "family": "flare", "params": {"em": False, "gen": 10, "states": 20}},
            {"id": "gnomix01", "family": "gnomix", "params": {"window_size_cM": .5, "crf_c2": .1}},
            {"id": "gnomix1", "family": "gnomix", "params": {"window_size_cM": .5, "crf_c2": 1}},
            {"id": "frequency", "family": "nnls", "params": {}},
            {"id": "pca_linear", "family": "pca_ridge", "params": {"n_components": 2, "alpha": 1}},
            {"id": "pca_forest", "family": "pca_rf", "params": {"n_components": 2, "n_estimators": 5, "min_samples_leaf": 1}},
            {"id": "pca_boost", "family": "pca_xgboost", "params": {"n_components": 2, "n_estimators": 5, "max_depth": 2}},
            {"id": "pca_clusters", "family": "pca_gmm", "params": {"n_components": 2, "k": 2}},
        ]
    config = {
        "version": 1, "seed": seed, "mode": "smoke", "sources": sources, "chromosomes": [22],
        "spans": {"22": span}, "samples": str(output / "samples.tsv"), "segments": str(output / "segments.tsv"),
        "reference": {"22": str(output / "reference.vcf.gz")}, "training": {"22": str(output / "training.vcf.gz")},
        "maps": {"22": str(output / "chr22.map")},
        "routes": [{"id": "array", "stage": "harmonized", "vcfs": {"22": str(output / "harmonized.vcf.gz")}},
                   {"id": "imputed", "stage": "imputed", "vcfs": {"22": str(output / "query.vcf.gz")}}],
        "models": models,
        "interpretations": [
            {"id": "native", "kind": "native", "params": {}},
            {"id": "simplex", "kind": "simplex", "params": {"l2": .001}},
            {"id": "linear", "kind": "ridge", "params": {"alpha": 1}},
            {"id": "forest", "kind": "rf", "params": {"n_estimators": 5, "min_samples_leaf": 1}},
            {"id": "boost", "kind": "xgboost", "params": {"n_estimators": 5, "max_depth": 2}},
        ],
        "tools": tools, "reference_filters": {"min_maf": 0, "min_mac": 0, "minimum_markers": 100},
        "allow_research_only": True, "independent_test": False,
    }
    write(output / "experiment.json", config)
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--flare-jar", required=True)
    parser.add_argument("--beagle-jar", required=True)
    parser.add_argument("--gnomix-root", required=True)
    args = parser.parse_args()
    fixture(args.output.resolve(), {"flare_jar": args.flare_jar, "beagle_jar": args.beagle_jar,
                                   "gnomix_root": args.gnomix_root})
