#!/usr/bin/env python3
"""Academic-only adapter to pinned upstream Gnomix; no ancestry simulation.

Build Dockerfile.gnomix-screen, then run this file with --outdir /results
--candidate gnomix_{fast,default,best}_{020,050} [--smoke].
The default stage fits and evaluates; --stage fit/evaluate separates them.
Smoke uses chromosome 22 only and is a compatibility test, never a ranking.

Gnomix uses equal-SNP-count windows, NOT equal-genetic-distance windows.
We use its exact target-cM-to-SNP conversion (max map cM), report actual spans, and
fail on unsupported density. XGB uses 75/31 windows (~15 cM); upstream CRF is
a whole-chromosome linear chain and ignores smooth_size. This limitation is
recorded, not silently replaced with an independent implementation.

predictions.tsv is marker-weighted (the primary frozen-baseline estimand).
predictions_physical.tsv uses marker midpoint cells extending to the supplied
chromosome span endpoints. Each includes both paired haplotypes equally.
No evaluation labels/segments are consumed, including by preprocessing.
Upstream's academic-research-only LICENSE.md applies; see third_party/gnomix.
"""

import argparse
import ast
import csv
from functools import lru_cache, partial
import gzip
import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import pickle
import re
import subprocess
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

SHA = "cd15f65eadfca8e2bbf209eb69a82978a200c5fd"
SEED = 20261013
SOURCES = ("EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA")
PROJECT = Path(__file__).resolve().parents[2]
UPSTREAM = PROJECT / "third_party/gnomix"
LOG = logging.getLogger("gnomix-screen")


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    path = Path(path)
    staging = path.with_name(path.name + ".pending")
    staging.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    staging.replace(path)


def config_for(candidate):
    match = re.fullmatch(r"gnomix_(fast|default|best)_(020|050)", candidate)
    if not match:
        raise ValueError(f"Unsupported candidate: {candidate}")
    mode, size = match.groups()
    return dict(candidate=candidate, architecture=mode, window_size_cM=.2 if size == "020" else .5,
                smooth_size=75 if size == "020" else 31, context_ratio=.5,
                calibrate=False, retrain_base=False, seed=SEED, sources=list(SOURCES))


def load_upstream():
    actual = subprocess.check_output(["git", "-C", str(UPSTREAM), "rev-parse", "HEAD"],
                                     text=True).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(UPSTREAM), "status", "--porcelain", "--untracked-files=all"], text=True)
    if actual != SHA or dirty:
        raise ValueError(f"Gnomix must be the clean pinned upstream {SHA}; got {actual}, {dirty!r}")
    if "academic research use only" not in (UPSTREAM / "LICENSE.md").read_text():
        raise ValueError("Missing upstream academic licensing notice")
    sys.path.insert(0, str(UPSTREAM))
    from src.model import Gnomix
    from src.preprocess import WindowedNpyData
    return Gnomix, WindowedNpyData


def records(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_records(path, labelled=False):
    with open(path, newline="") as handle:
        rows = [row if labelled else {"sample_id": row["sample_id"]}
                for row in csv.DictReader(handle, delimiter="\t")]
    ids = [row["sample_id"] for row in rows]
    if not ids or len(set(ids)) != len(ids) or any(not x for x in ids):
        raise ValueError(f"Empty/duplicate sample identifiers: {path}")
    if labelled and {row["source_group"] for row in rows} - set(SOURCES):
        raise ValueError(f"Unknown source in {path}")
    return rows


def validate_splits(manifest):
    if manifest["sources"] != list(SOURCES):
        raise ValueError("Unexpected source order")
    if manifest.get("seed", SEED) != SEED:
        raise ValueError("Unexpected screen seed")
    groups = {}
    for name in ("reference", "pure_train"):
        groups[name] = {r["sample_id"] for r in sample_records(
            manifest[name]["labels_tsv"], labelled=True)}
        found = {r["source_group"] for r in sample_records(manifest[name]["labels_tsv"], True)}
        if found != set(SOURCES):
            raise ValueError(f"{name} must include all seven sources")
    groups["pedigree_train"] = set()
    for batch in manifest["pedigree_train"]["batches"]:
        ids = {r["sample_id"] for r in sample_records(batch["samples_tsv"])}
        if groups["pedigree_train"] & ids:
            raise ValueError("Duplicate pedigree-training sample across batches")
        groups["pedigree_train"].update(ids)
    for name in ("pure", "pedigree", "real"):
        groups[name] = set()
        for batch in manifest["evaluation"][name]["batches"]:
            ids = {r["sample_id"] for r in sample_records(batch["samples_tsv"])}
            if groups[name] & ids:
                raise ValueError(f"Duplicate evaluation {name} sample")
            groups[name].update(ids)
    for left, a in groups.items():
        for right, b in groups.items():
            if left < right and a & b:
                raise ValueError(f"Sample overlap {left}/{right}: {sorted(a & b)[:3]}")
    return {key: len(value) for key, value in groups.items()}


@lru_cache(maxsize=128)
def vcf_header(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                contigs.append(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
            if line.startswith("#CHROM"):
                return tuple(line.rstrip().split("\t")[9:]), tuple(contigs)
    raise ValueError(f"VCF header missing: {path}")


def phased_gt(text):
    fields = text.split("|")
    if len(fields) != 2 or any(x not in ("0", "1", ".") for x in fields):
        raise ValueError(f"Require phased biallelic diploid GT, got {text!r}")
    return [2 if x == "." else int(x) for x in fields]


def variants(path, chrom, ids):
    """Stream selected genotypes in requested sample order; never infer phase."""
    header, contigs = vcf_header(str(path))
    if len(set(ids)) != len(ids) or set(ids) - set(header):
        raise ValueError(f"Missing/duplicate requested samples in {path}")
    region = str(chrom) if str(chrom) in contigs else "chr" + str(chrom)
    if region not in contigs:
        raise ValueError(f"Chromosome {chrom} absent from contig header: {path}")
    command = ["bcftools", "query", "-H", "-r", region, "-s", ",".join(ids),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    reorder = None
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                selected = [re.sub(r"^\[\d+\]", "", field).removesuffix(":GT")
                            for field in fields[4:]]
                if set(selected) != set(ids) or len(selected) != len(ids):
                    raise ValueError("bcftools query sample header mismatch")
                reorder = np.asarray([selected.index(x) for x in ids])
                continue
            if reorder is None:
                raise ValueError("bcftools query did not emit its sample-axis header")
            if len(fields) != 4 + len(ids):
                raise ValueError("Unexpected bcftools GT column count")
            key = (fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3])
            yield key, [fields[4 + i] for i in reorder]
        error = process.stderr.read()
        if process.wait():
            raise ValueError(f"bcftools query failed for {path}: {error.strip()}")
    finally:
        if process.poll() is None:
            process.terminate()
        process.wait()
        process.stdout.close()
        process.stderr.close()


def reference_markers(path, ids, chrom, filters):
    keys = []
    last = -1
    min_maf, min_mac = filters["min_maf"], filters["min_mac"]
    if not 0 <= min_maf <= .5 or min_mac < 0:
        raise ValueError("Invalid reference-only allele filters")
    for key, genotypes in variants(path, chrom, ids):
        if key[1] <= last:
            raise ValueError("VCF positions must be unique and increasing")
        last = key[1]
        if "," in key[3] or key[3] in (".", "*"):
            continue
        alleles = np.asarray([phased_gt(x) for x in genotypes], dtype=np.int8)
        called = alleles != 2
        total = int(called.sum())
        alt = int((alleles == 1).sum())
        mac = min(alt, total - alt)
        if total and mac >= min_mac and mac / total >= min_maf:
            keys.append(key)
    if not keys:
        raise ValueError(f"No retained markers for chromosome {chrom}")
    return keys


def check_shared_markers(path, chrom, keys):
    shared = []
    with open(path) as handle:
        for line in handle:
            chm, position, ref, alt = line.rstrip().split("\t")
            if chm.removeprefix("chr") == str(chrom):
                shared.append((str(chrom), int(position), ref, alt))
    if keys != shared:
        raise ValueError(f"Reference-filtered markers differ from frozen shared markers on chr{chrom}")


def load_haplotypes(path, ids, chrom, keys, output=None):
    shape = (2 * len(ids), len(keys))
    result = (np.lib.format.open_memmap(output, mode="w+", dtype=np.int8, shape=shape)
              if output else np.empty(shape, dtype=np.int8))
    wanted = {key[1]: (i, key) for i, key in enumerate(keys)}
    seen = np.zeros(len(keys), dtype=bool)
    for key, genotypes in variants(path, chrom, ids):
        item = wanted.get(key[1])
        if item is None:
            continue
        index, expected = item
        if key != expected or seen[index]:
            raise ValueError(f"Exact chr/POS/REF/ALT mismatch or duplicate: {key} vs {expected}")
        result[:, index] = np.asarray([phased_gt(x) for x in genotypes], dtype=np.int8).reshape(-1)
        seen[index] = True
    if not np.all(seen):
        raise ValueError(f"{path}: missing {int((~seen).sum())} required markers on chr{chrom}")
    if isinstance(result, np.memmap):
        result.flush()
    return result


def map_points(path, chrom):
    values = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].removeprefix("chr") != str(chrom):
                continue
            if len(fields) != 4:
                raise ValueError("Genetic map must have PLINK chrom/rsid/cM/pos_bp columns")
            values.append((int(fields[3]), float(fields[2])))
    data = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.all(np.isfinite(data)):
        raise ValueError(f"Missing/non-finite map for chr{chrom}")
    data = data[np.argsort(data[:, 0], kind="stable")]
    if np.any(np.diff(data[:, 0]) <= 0) or np.any(np.diff(data[:, 1]) < 0):
        raise ValueError(f"Map must have unique increasing positions and nondecreasing cM: chr{chrom}")
    return data


def window_layout(keys, points, span, config):
    pos = np.asarray([key[1] for key in keys], dtype=np.int64)
    if span[0] > pos[0] or span[1] <= pos[-1] or span[0] >= span[1]:
        raise ValueError("Retained markers must be inside the supplied half-open chromosome span")
    if pos[0] < points[0, 0] or pos[-1] > points[-1, 0]:
        raise ValueError("Genetic map does not bracket all retained markers (no extrapolation)")
    cm = np.interp(pos, points[:, 0], points[:, 1])
    genetic_length = float(points[:, 1].max())
    if genetic_length <= 0:
        raise ValueError("Nonpositive retained-marker genetic span")
    # Upstream LAIDataset uses max(map cM), not last-minus-first retained cM.
    size = int(round(config["window_size_cM"] * len(pos) / genetic_length))
    if size < 2:
        raise ValueError("SNP density too low for a >=2-SNP window and context_ratio=.5")
    count = len(pos) // size
    if count < 2 or (config["architecture"] != "fast" and count < 2 * config["smooth_size"]):
        raise ValueError(f"Unsupported sparse/short chromosome: W={count}, S={config['smooth_size']}")
    starts = np.arange(count, dtype=np.int64) * size
    ends = np.r_[starts[1:], len(pos)]
    midpoint = (pos[:-1].astype(float) + pos[1:]) / 2
    physical_cells = np.diff(np.r_[float(span[0]), midpoint, float(span[1])])
    marker_weights = ends - starts
    physical_weights = np.asarray([physical_cells[a:b].sum() for a, b in zip(starts, ends)])
    if np.any(marker_weights <= 0) or np.any(physical_weights <= 0):
        raise ValueError("Empty or nonpositive-span window")
    genetic_spans = [float(cm[b - 1] - cm[a]) for a, b in zip(starts, ends)]
    return dict(M=size, W=count, starts=starts, ends=ends, marker_weights=marker_weights,
                physical_weights=physical_weights, genetic_span_cM=genetic_length,
                window_genetic_spans_cM=genetic_spans, positions=pos)


def segment_labels(path, ids, chrom, positions):
    """One assignment per retained marker, preserving sample/haplotype pairing."""
    wanted = set(ids)
    groups = {}
    for row in records(path):
        if row["sample_id"] not in wanted or row["chrom"].removeprefix("chr") != str(chrom):
            continue
        hap = int(row["haplotype"])
        start, end = int(row["start_bp"]), int(row["end_bp"])
        if hap not in (1, 2) or start < 0 or end <= start or row["source_group"] not in SOURCES:
            raise ValueError(f"Invalid ancestry segment: {row}")
        groups.setdefault((row["sample_id"], hap), []).append(
            (start, end, SOURCES.index(row["source_group"])))
    result = np.empty((2 * len(ids), len(positions)), dtype=np.int8)
    for i, sample in enumerate(ids):
        for hap in (1, 2):
            values = sorted(groups.get((sample, hap), []))
            covered = np.zeros(len(positions), dtype=np.int8)
            previous_end = -1
            for start, end, source in values:
                if start < previous_end:
                    raise ValueError(f"Overlapping segments: {sample}/hap{hap}/chr{chrom}")
                previous_end = end
                left, right = np.searchsorted(positions, [start, end], side="left")
                result[2 * i + hap - 1, left:right] = source
                covered[left:right] += 1
            if not np.all(covered == 1):
                raise ValueError(f"Incomplete haplotype truth: {sample}/hap{hap}/chr{chrom}")
    return result


def window_labels(labels, layout):
    # Upstream scipy.stats.mode semantics: lowest numeric class wins ties.
    return np.asarray([[np.bincount(row[a:b], minlength=len(SOURCES)).argmax()
                        for a, b in zip(layout["starts"], layout["ends"])]
                       for row in labels], dtype=np.int16)


def ordered_probabilities(model, base):
    result = np.asarray(model.smooth.predict_proba(base), dtype=float)
    order = [int(x) for x in model.smooth.model.classes_]
    if sorted(order) != list(range(len(SOURCES))):
        raise ValueError(f"Unexpected smoother class axis {order}")
    result = result[:, :, np.argsort(order)]
    if (result.shape != (len(base), model.W, len(SOURCES)) or not np.all(np.isfinite(result))
            or np.any(result < 0) or not np.allclose(result.sum(axis=2), 1, atol=2e-6)):
        raise ValueError("Invalid posterior probability tensor")
    return result / result.sum(axis=2, keepdims=True)


def aggregate_probabilities(probabilities, layout):
    if len(probabilities) % 2:
        raise ValueError("Unpaired haplotypes")
    paired = probabilities.reshape(-1, 2, layout["W"], len(SOURCES)).mean(axis=1)
    primary = np.einsum("nwa,w->na", paired, layout["marker_weights"])
    physical = np.einsum("nwa,w->na", paired, layout["physical_weights"])
    return primary, physical


def initialize_model(config, keys, layout, ncores, gnomix_class):
    model = gnomix_class(C=len(keys), M=layout["M"], A=len(SOURCES), S=config["smooth_size"],
                         mode=config["architecture"], context_ratio=.5, calibrate=False,
                         seed=SEED, n_jobs=ncores, missing_encoding=2,
                         snp_pos=layout["positions"], snp_ref=np.asarray([k[2] for k in keys]),
                         snp_alt=np.asarray([k[3] for k in keys]),
                         population_order=np.asarray(SOURCES), mode_filter=0, verbose=False)
    # Sequential disk-window training bounds memory and avoids nested process pools.
    model.base.base_multithread = False
    for estimator in model.base.models:
        estimator.set_params(random_state=SEED)
    if config["architecture"] == "best":
        from src.Base.string_kernel import (CovRSK_DP_triangular_numbers,
                                            CovRSK_DP_triangular_numbers_multithread)
        kernel = (partial(CovRSK_DP_triangular_numbers, seed=SEED) if layout["M"] < 500 else
                  partial(CovRSK_DP_triangular_numbers_multithread, seed=SEED, n_jobs=ncores))
        model.base.kernel = kernel
        for estimator in model.base.models:
            estimator.set_params(kernel=kernel)
    return model


def chunks(rows, batch_size):
    for start in range(0, len(rows), batch_size):
        yield rows[start:start + batch_size]


def fit_chromosome(manifest, config, chrom, directory, ncores, batch_size, upstream):
    gnomix_class, windowed_class = upstream
    reference = manifest["reference"]
    ref_rows = sample_records(reference["labels_tsv"], labelled=True)
    ref_ids = [r["sample_id"] for r in ref_rows]
    filters = manifest.get("reference_filters", dict(min_maf=.01, min_mac=20))
    keys = reference_markers(reference["vcf"], ref_ids, chrom, filters)
    if manifest.get("markers_tsv"):
        check_shared_markers(manifest["markers_tsv"], chrom, keys)
    layout = window_layout(keys, map_points(manifest["genetic_map"], chrom),
                           manifest["chromosome_spans"][str(chrom)], config)
    LOG.info("chr%s: reference=%d, markers=%d, M=%d SNPs, W=%d, smoother=%s S=%d",
             chrom, len(ref_ids), len(keys), layout["M"], layout["W"],
             config["architecture"], config["smooth_size"])
    work = directory / "working"
    work.mkdir(exist_ok=True)
    x_path, y_path = work / "reference_X.npy", work / "reference_labels.npy"
    x = load_haplotypes(reference["vcf"], ref_ids, chrom, keys, x_path)
    del x
    y = np.lib.format.open_memmap(y_path, mode="w+", dtype=np.int8,
                                shape=(2 * len(ref_ids), len(keys)))
    y[:] = np.repeat([SOURCES.index(r["source_group"]) for r in ref_rows], 2)[:, None]
    y.flush()
    del y
    model = initialize_model(config, keys, layout, ncores, gnomix_class)
    data = windowed_class([str(x_path)], [str(y_path)], len(keys), layout["M"])
    begin = time.monotonic()
    model.base.train_windowed(data, verbose=False)
    LOG.info("chr%s: base trained in %.1fs", chrom, time.monotonic() - begin)
    base_probabilities, labels = [], []
    pure = manifest["pure_train"]
    training = [(pure, sample_records(pure["labels_tsv"], True), True)]
    training.extend((batch, sample_records(batch["samples_tsv"]), False)
                    for batch in manifest["pedigree_train"]["batches"])
    for batch, rows, is_pure in training:
        for part in chunks(rows, batch_size):
            ids = [r["sample_id"] for r in part]
            x = load_haplotypes(batch["vcf"], ids, chrom, keys)
            if is_pure:
                y_window = np.repeat([SOURCES.index(r["source_group"]) for r in part], 2)
                y_window = np.repeat(y_window[:, None], layout["W"], axis=1)
            else:
                y_window = window_labels(
                    segment_labels(batch["segments_tsv"], ids, chrom, layout["positions"]), layout)
            base_probabilities.append(model.base.predict_proba(x).astype(np.float32))
            labels.append(y_window)
            del x
    b, y = np.concatenate(base_probabilities), np.concatenate(labels)
    del base_probabilities, labels
    if set(np.unique(y)) != set(range(len(SOURCES))):
        raise ValueError("Smoother training does not contain all seven sources")
    model.smooth.train(b, y)
    del b, y
    # No Gnomix.train() re-training/calibration/evaluation stage: the supported
    # base/smooth stage APIs keep the reference-only base distribution frozen.
    bundle = dict(model=model, keys=keys, layout=layout, config=config, chromosome=chrom,
                  upstream_sha=SHA)
    model_path = directory / "models" / f"chr{chrom}.pkl"
    model_path.parent.mkdir(exist_ok=True)
    with open(model_path.with_suffix(".pending"), "wb") as handle:
        pickle.dump(bundle, handle, protocol=pickle.HIGHEST_PROTOCOL)
    model_path.with_suffix(".pending").replace(model_path)
    metadata = {k: v for k, v in layout.items()
                if k not in ("starts", "ends", "positions", "marker_weights", "physical_weights")}
    metadata.update(chromosome=chrom, n_markers=len(keys), n_reference=len(ref_ids),
                    marker_sha256=hashlib.sha256(
                        "".join("\t".join(map(str, k)) + "\n" for k in keys).encode()).hexdigest(),
                    model_sha256=digest(model_path), fit_seconds=time.monotonic() - begin,
                    effective_smoother=("whole-chromosome CRF (smooth_size ignored by upstream)"
                                        if config["architecture"] == "fast" else
                                        f"{config['smooth_size']} adjacent equal-SNP windows"))
    write_json(directory / "models" / f"chr{chrom}.json", metadata)
    x_path.unlink()
    y_path.unlink()
    work.rmdir()
    LOG.info("chr%s: saved model; total fit %.1fs", chrom, metadata["fit_seconds"])
    return bundle


def evaluate_chromosome(manifest, bundle, reloaded, batch_size):
    if (reloaded["config"] != bundle["config"] or reloaded["keys"] != bundle["keys"]
            or reloaded["upstream_sha"] != SHA):
        raise ValueError("Reloaded model metadata differs")
    primary, physical = {}, {}
    chrom, layout, keys = bundle["chromosome"], bundle["layout"], bundle["keys"]
    max_delta = 0.
    for dataset in ("pure", "pedigree", "real"):
        for batch in manifest["evaluation"][dataset]["batches"]:
            # Deliberately read only sample IDs, never evaluation ancestry truth.
            rows = sample_records(batch["samples_tsv"])
            for part in chunks(rows, batch_size):
                ids = [row["sample_id"] for row in part]
                x = load_haplotypes(batch["vcf"], ids, chrom, keys)
                before = ordered_probabilities(bundle["model"],
                                               bundle["model"].base.predict_proba(x))
                after = ordered_probabilities(reloaded["model"],
                                              reloaded["model"].base.predict_proba(x))
                delta = float(np.max(np.abs(before - after)))
                max_delta = max(max_delta, delta)
                if not np.allclose(before, after, rtol=0, atol=1e-10):
                    raise ValueError(f"Model reload changed chr{chrom} predictions by {delta}")
                a, b = aggregate_probabilities(after, layout)
                for i, sample in enumerate(ids):
                    key = (dataset, sample)
                    if key in primary:
                        raise ValueError(f"Duplicate sample in evaluation batches: {key}")
                    primary[key], physical[key] = a[i], b[i]
                del x, before, after
    LOG.info("chr%s: evaluated %d samples; reload max error %.3g", chrom, len(primary), max_delta)
    return primary, physical, max_delta


def input_hashes(manifest):
    paths = set()
    def collect(value):
        if isinstance(value, dict):
            for child in value.values():
                collect(child)
        elif isinstance(value, list):
            for child in value:
                collect(child)
        elif isinstance(value, str) and value.startswith("/") and Path(value).is_file():
            paths.add(value)
    collect(manifest)
    for path in list(paths):
        if path.endswith(".vcf.gz"):
            paths.update(str(index) for suffix in (".csi", ".tbi")
                         if (index := Path(path + suffix)).is_file())
    return {path: digest(path) for path in sorted(paths)}


def write_predictions(path, totals, denominator):
    with open(path.with_suffix(".pending"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset", "sample_id"] + ["pred_" + s for s in SOURCES])
        for (dataset, sample), values in sorted(totals.items()):
            values = values / denominator
            if not np.all(np.isfinite(values)) or not np.isclose(values.sum(), 1, atol=2e-6):
                raise ValueError("Invalid genomic posterior aggregation")
            writer.writerow([dataset, sample] + [f"{x:.15g}" for x in values])
    path.with_suffix(".pending").replace(path)


def implementation_fingerprint(path):
    excluded = {"run", "main", "implementation_fingerprint", "saved_chromosome"}
    tree = ast.parse(Path(path).read_text())
    tree.body = [node for node in tree.body
                 if not (isinstance(node, ast.FunctionDef) and node.name in excluded)
                 and not (isinstance(node, ast.Import)
                          and [alias.name for alias in node.names] == ["ast"])]
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest()


def saved_chromosome(manifest, config, chrom, directory):
    path = directory / "models" / f"chr{chrom}.pkl"
    metadata = json.loads(path.with_suffix(".json").read_text())
    if metadata["chromosome"] != chrom or digest(path) != metadata["model_sha256"]:
        raise ValueError(f"Invalid saved chromosome {chrom} checksum/identity")
    with path.open("rb") as handle:
        bundle = pickle.load(handle)
    if (bundle["config"] != config or bundle["chromosome"] != chrom
            or bundle["upstream_sha"] != SHA):
        raise ValueError(f"Saved chromosome {chrom} has incompatible model configuration")
    keys = bundle["keys"]
    if manifest.get("markers_tsv"):
        check_shared_markers(manifest["markers_tsv"], chrom, keys)
    marker_sha = hashlib.sha256(
        "".join("\t".join(map(str, key)) + "\n" for key in keys).encode()).hexdigest()
    if metadata["n_markers"] != len(keys) or metadata["marker_sha256"] != marker_sha:
        raise ValueError(f"Saved chromosome {chrom} marker metadata changed")
    expected = window_layout(keys, map_points(manifest["genetic_map"], chrom),
                             manifest["chromosome_spans"][str(chrom)], config)
    if (set(bundle["layout"]) != set(expected)
            or any(not np.array_equal(bundle["layout"][key], value)
                   for key, value in expected.items())):
        raise ValueError(f"Saved chromosome {chrom} window layout changed")
    return bundle, metadata


def run(outdir, candidate, smoke=False, stage="run", ncores=6, batch_size=64,
        resume=False, verify_only=False):
    if not 1 <= ncores <= 6 or not 1 <= batch_size <= 256:
        raise ValueError("Require 1<=ncores<=6 and 1<=batch_size<=256")
    if (resume or verify_only) and stage != "run":
        raise ValueError("Resume/verification requires the combined run stage")
    if verify_only and not resume:
        raise ValueError("Verification requires --resume")
    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                      MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    np.random.seed(SEED)
    outdir = Path(outdir)
    config = config_for(candidate)
    manifest_path = outdir / "cohort.json"
    manifest = json.loads(manifest_path.read_text())
    counts = validate_splits(manifest)
    upstream = load_upstream()
    configurations_path = outdir / "configurations.json"
    if configurations_path.exists():
        declared = [c for c in json.loads(configurations_path.read_text())
                    if c["candidate"] == candidate]
        if len(declared) != 1 or any(declared[0].get(k) != v for k, v in config.items()):
            raise ValueError("Candidate does not match frozen configurations.json")
    directory = outdir / ("smoke" if smoke else "candidates") / candidate
    chromosomes = [22] if smoke else list(range(1, 23))
    LOG.info("Hashing cohort inputs before %s; samples=%s", stage, counts)
    hashes = input_hashes(manifest)
    audit = dict(candidate=candidate, config=config, upstream_sha=SHA,
                 upstream_tree_sha=subprocess.check_output(
                     ["git", "-C", str(UPSTREAM), "rev-parse", "HEAD^{tree}"], text=True).strip(),
                 upstream_license_sha256=digest(UPSTREAM / "LICENSE.md"),
                 adapter_sha256=digest(__file__), manifest_sha256=digest(manifest_path),
                 execution_fingerprint=implementation_fingerprint(__file__),
                 requirements_sha256=digest(PROJECT / "requirements_gnomix_screen.txt"),
                 dockerfile_sha256=digest(PROJECT / "Dockerfile.gnomix-screen"),
                 configurations_sha256=(digest(configurations_path)
                                        if configurations_path.exists() else None),
                 input_sha256=hashes, sample_counts=counts, chromosomes=chromosomes,
                 smoke_only=smoke, status="running", stage=stage, ncores=ncores,
                 batch_size=batch_size, no_synthetic_mosaics=True,
                 calibrate=False, retrain_base=False, evaluation_truth_used_in_fit=False,
                 reference_filters=manifest.get("reference_filters",
                                                dict(min_maf=.01, min_mac=20)),
                 primary_aggregation="retained-marker-weighted paired-haplotype posterior mean",
                 secondary_aggregation="physical midpoint cells to supplied chromosome endpoints",
                 windowing="upstream equal-SNP count round(target_cM*C/max_map_cM); last remainder",
                 crf_span_limitation="upstream CRF is a whole-chromosome chain; smooth_size unused",
                 dependencies={name: importlib.metadata.version(name) for name in
                               ("numpy", "scipy", "scikit-learn", "pandas", "xgboost-cpu",
                                "sklearn-crfsuite", "scikit-allel", "uncertainty-calibration")},
                 chromosome_audits=[])
    audit_path = directory / "audit.json"
    checkpoints = []
    previous = None
    if stage == "evaluate" or (resume and audit_path.exists()):
        previous = json.loads(audit_path.read_text())
        for key in ("candidate", "config", "upstream_sha", "upstream_tree_sha",
                    "upstream_license_sha256", "manifest_sha256", "requirements_sha256",
                    "dockerfile_sha256", "configurations_sha256", "input_sha256",
                    "sample_counts", "chromosomes", "smoke_only", "dependencies",
                    "ncores", "batch_size"):
            if previous[key] != audit[key]:
                raise ValueError(f"Fit/evaluate provenance mismatch: {key}")
        if previous["adapter_sha256"] != audit["adapter_sha256"]:
            original = outdir / "revisions/model_family_gnomix.before_resume.py"
            fingerprint = previous.get("execution_fingerprint")
            if fingerprint is None:
                if not original.is_file() or digest(original) != previous["adapter_sha256"]:
                    raise ValueError("Original adapter revision is missing or changed")
                fingerprint = implementation_fingerprint(original)
            if fingerprint != audit["execution_fingerprint"]:
                raise ValueError("Saved models' fitting/evaluation implementation changed")
        for chrom in chromosomes:
            path = directory / "models" / f"chr{chrom}.pkl"
            if path.with_suffix(".json").exists():
                bundle, metadata = saved_chromosome(manifest, config, chrom, directory)
                checkpoints.append(chrom)
                del bundle
            elif path.exists():
                LOG.warning("chr%s has an uncommitted model without metadata; will refit", chrom)
        if checkpoints != chromosomes[:len(checkpoints)]:
            raise ValueError("Saved chromosomes are not a contiguous checkpoint prefix")
        if not {item["chromosome"] for item in previous["chromosome_audits"]} <= set(checkpoints):
            raise ValueError("An audited chromosome checkpoint is missing")
        if stage == "evaluate" and checkpoints != chromosomes:
            raise ValueError("Standalone evaluation requires all fitted chromosomes")
        complete = previous["status"] == "complete" and previous["stage"] != "fit"
        if complete:
            output = json.loads((directory / "output_audit.json").read_text())
            if output != previous or checkpoints != chromosomes or not output["model_roundtrip_verified"]:
                raise ValueError("Completed candidate has inconsistent output provenance")
            expected_ids = {(name, row["sample_id"])
                            for name, group in manifest["evaluation"].items()
                            for batch in group["batches"]
                            for row in sample_records(batch["samples_tsv"])}
            for name, key in (("predictions.tsv", "prediction_sha256"),
                              ("predictions_physical.tsv", "physical_prediction_sha256")):
                path = directory / name
                if digest(path) != output[key]:
                    raise ValueError(f"Completed candidate output changed: {name}")
                rows = records(path)
                ids = {(row["dataset"], row["sample_id"]) for row in rows}
                values = np.asarray([[float(row["pred_" + source]) for source in SOURCES]
                                     for row in rows])
                if (ids != expected_ids or len(rows) != len(expected_ids)
                        or not np.all(np.isfinite(values)) or np.any(values < 0)
                        or not np.allclose(values.sum(axis=1), 1, atol=2e-6)):
                    raise ValueError(f"Completed candidate predictions invalid: {name}")
            if resume:
                LOG.warning("%s: verified completed candidate; no fitting required", candidate)
                return directory
        if verify_only:
            LOG.warning("%s: verified %d saved chromosome checkpoints; no fitting requested",
                        candidate, len(checkpoints))
            return directory
        audit["resume_history"] = previous.get("resume_history", []) + [
            dict(audit_sha256=digest(audit_path), adapter_sha256=previous["adapter_sha256"],
                 status=previous["status"], saved_chromosomes=checkpoints)]
    elif verify_only:
        raise ValueError("No candidate audit exists to verify")
    elif resume and directory.exists() and any((directory / "models").glob("*")):
        raise ValueError("Cannot resume model files without their fit provenance")
    directory.mkdir(parents=True, exist_ok=True)
    for filename in ("predictions.tsv", "predictions_physical.tsv", "output_audit.json"):
        (directory / filename).unlink(missing_ok=True)
    logfile = logging.FileHandler(directory / "adapter.log")
    logfile.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    LOG.addHandler(logfile)
    write_json(audit_path, audit)
    primary_total, physical_total = {}, {}
    marker_total, span_total = 0, 0.
    try:
        with threadpool_limits(limits=1):
            for chrom in chromosomes:
                if (outdir / "pause_requested.json").exists():
                    audit["status"] = "paused"
                    write_json(audit_path, audit)
                    LOG.warning("%s: paused at saved chromosome boundary", candidate)
                    raise SystemExit(75)
                if stage == "evaluate" or chrom in checkpoints:
                    bundle, metadata = saved_chromosome(manifest, config, chrom, directory)
                    LOG.info("chr%s: reusing verified saved model; no fitting", chrom)
                else:
                    bundle = fit_chromosome(manifest, config, chrom, directory,
                                            ncores, batch_size, upstream)
                model_path = directory / "models" / f"chr{chrom}.pkl"
                metadata = json.loads(model_path.with_suffix(".json").read_text())
                if digest(model_path) != metadata["model_sha256"]:
                    raise ValueError("Model file changed since fitting")
                with open(model_path, "rb") as handle:
                    reloaded = pickle.load(handle)
                if stage != "fit":
                    primary, physical, delta = evaluate_chromosome(
                        manifest, bundle, reloaded, batch_size)
                    if primary_total and set(primary_total) != set(primary):
                        raise ValueError("Evaluation sample coverage changed between chromosomes")
                    for key in primary:
                        primary_total[key] = primary_total.get(key, np.zeros(len(SOURCES))) + primary[key]
                        physical_total[key] = physical_total.get(key, np.zeros(len(SOURCES))) + physical[key]
                    marker_total += len(bundle["keys"])
                    span_total += float(bundle["layout"]["physical_weights"].sum())
                    metadata["reload_max_absolute_error"] = delta
                audit["chromosome_audits"].append(metadata)
                del bundle, reloaded
                write_json(audit_path, audit)
                if (outdir / "pause_requested.json").exists():
                    audit["status"] = "paused"
                    write_json(audit_path, audit)
                    LOG.warning("%s: paused after saving/evaluating chr%s", candidate, chrom)
                    raise SystemExit(75)
        # Rehash after the run; never accept outputs from changing cohort files.
        if (input_hashes(manifest) != hashes or digest(manifest_path) != audit["manifest_sha256"]
                or digest(__file__) != audit["adapter_sha256"]
                or (configurations_path.exists() and
                    digest(configurations_path) != audit["configurations_sha256"])):
            raise ValueError("Cohort inputs changed during fitting/evaluation")
        if stage != "fit":
            write_predictions(directory / "predictions.tsv", primary_total, marker_total)
            write_predictions(directory / "predictions_physical.tsv", physical_total, span_total)
            audit.update(model_roundtrip_verified=True, n_predictions=len(primary_total),
                         n_markers=marker_total, physical_span_bp=span_total,
                         prediction_sha256=digest(directory / "predictions.tsv"),
                         physical_prediction_sha256=digest(directory / "predictions_physical.tsv"))
        audit["status"] = "complete"
        write_json(audit_path, audit)
        if stage != "fit":
            write_json(directory / "output_audit.json", audit)
        return directory
    except Exception as error:
        audit.update(status="failed", error=str(error))
        write_json(audit_path, audit)
        raise
    finally:
        LOG.removeHandler(logfile)
        logfile.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--stage", choices=("run", "fit", "evaluate"), default="run")
    parser.add_argument("--smoke", action="store_true", help="Chromosome 22 compatibility only; not ranking")
    parser.add_argument("--resume", action="store_true",
                        help="Verify and reuse completed candidates/chromosomes; fit only missing ones")
    parser.add_argument("--verify-only", action="store_true",
                        help="With --resume, verify saved work without fitting or changing outputs")
    parser.add_argument("--ncores", type=int, default=6)
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    LOG.info("Academic non-commercial research only; pinned Gnomix %s", SHA)
    run(args.outdir, args.candidate, args.smoke, args.stage, args.ncores, args.batch_size,
        args.resume, args.verify_only)


if __name__ == "__main__":
    main()
