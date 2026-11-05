#!/usr/bin/env python3
"""Source-anchored FLARE/Gnomix development screen on audited, reused cohorts."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from realized_ancestry import load_bp
from screen_precision_retraining import common
from screen_v1_decoders import (
    SOURCES as ORIGINAL_SOURCES, SEED, digest, fractions, frozen_prediction, write_json,
)
from select_unified_ancestry import weighted_quantile
from fit_evaluate_candidate import evaluate_constraints

SOURCES = ("EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA")
FLARE_MODES = tuple(f"flare_{structure}_{policy}"
                    for structure in ("direct", "fine", "subclusters")
                    for policy in ("fixed", "refined"))
GNOMIX_MODES = tuple(f"gnomix_{architecture}_{window}"
                     for architecture in ("fast", "default", "best")
                     for window in ("020", "050"))
MODES = FLARE_MODES + GNOMIX_MODES


def broad(source):
    if source not in ORIGINAL_SOURCES:
        raise ValueError(f"Unknown original source: {source}")
    return "EAS_SEA" if source in ("EAS", "SEA") else source


def collapse(values):
    values = fractions(values, "Eight-source proportions")
    return fractions(np.column_stack([
        values[:, ORIGINAL_SOURCES.index("EAS")] + values[:, ORIGINAL_SOURCES.index("SEA")]
        if source == "EAS_SEA" else values[:, ORIGINAL_SOURCES.index(source)]
        for source in SOURCES
    ]), "Seven-source proportions")


def samples(vcf):
    result = subprocess.run(["bcftools", "query", "-l", str(vcf)],
                            check=True, capture_output=True, text=True)
    ids = result.stdout.splitlines()
    if not ids or len(ids) != len(set(ids)):
        raise ValueError(f"Empty/duplicate VCF sample IDs: {vcf}")
    return ids


def run(command, directory, name):
    write_json(directory / f"{name}.command.json", list(map(str, command)))
    started = time.monotonic()
    with (directory / f"{name}.stdout.log").open("w") as stream:
        subprocess.run(list(map(str, command)), check=True, stdout=stream,
                       stderr=subprocess.STDOUT)
    print(f"{name} completed in {time.monotonic() - started:.1f}s", flush=True)


def table(frame, path):
    if frame.sample_id.isna().any() or frame.sample_id.duplicated().any():
        raise ValueError(f"Invalid sample metadata: {path}")
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def indexed_input(vcf, inputs, name):
    vcf = Path(vcf)
    if Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file():
        return vcf
    copied = inputs / f"{name}.vcf.gz"
    shutil.copyfile(vcf, copied)
    if digest(copied) != digest(vcf):
        raise ValueError("Reused genotype copy changed during transfer")
    run(["bcftools", "index", "-t", copied], inputs, f"index_{name}")
    return copied


def marker_table(vcf):
    import io
    result = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
                             str(vcf)], check=True, capture_output=True, text=True)
    return pd.read_csv(io.StringIO(result.stdout), sep="\t",
                       names=["chrom", "pos", "ref", "alt"])


def canonical_genotypes(vcf, reference, markers, inputs, name):
    target = inputs / f"{name}.canonical.vcf.gz"
    original = marker_table(vcf)
    run(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz", "-o", target,
         vcf, reference], inputs, f"canonical_{name}")
    run(["bcftools", "index", "-t", target], inputs, f"index_canonical_{name}")
    actual = marker_table(target)
    pd.testing.assert_frame_equal(actual, markers.reset_index(drop=True))
    if samples(target) != samples(vcf):
        raise ValueError("Canonicalization changed genotype sample identities/order")
    return str(target), {
        "original_sha256": digest(vcf), "canonical_sha256": digest(target),
        "original_records": len(original), "canonical_records": len(actual),
        "original_duplicate_position_records": int(original.duplicated(["chrom", "pos"]).sum()),
        "exact_ref_alt_identity_verified": True,
    }


def export_segments(directory, batch, selected, output):
    plan = directory / f"pedigree/plans/plans/{batch}/pedigrees.json"
    layout = json.loads(plan.read_text())
    spans, totals, tracked = {}, {}, [plan]
    fields = ["sample_id", "chrom", "haplotype", "start_bp", "end_bp", "source_group"]
    with output.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for chrom in range(1, 23):
            bp_path = directory / f"pedigree/chromosomes/pdb_{batch}_chr{chrom}.bp"
            tracked.append(bp_path)
            bp = load_bp(bp_path, chrom)
            ancestry = {}
            for founder, info in layout["founders"].items():
                for hap in ("h0", "h1"):
                    ids = {hid for _, _, hid in bp[founder][hap]}
                    if len(ids) != 1 or next(iter(ids)) in ancestry:
                        raise ValueError("Founder haplotype assignment is not unique")
                    ancestry[next(iter(ids))] = broad(info["source"])
            chrom_span = None
            for sample in selected.sample_id:
                if set(bp.get(sample, {})) != {"h0", "h1"}:
                    raise ValueError(f"Missing two-haplotype truth: {sample}, chr{chrom}")
                total = totals.setdefault(sample, np.zeros(len(SOURCES)))
                for number, hap in enumerate(("h0", "h1"), 1):
                    intervals = bp[sample][hap]
                    span = [intervals[0][0], intervals[-1][1] + 1]
                    if chrom_span is not None and span != chrom_span:
                        raise ValueError("Haplotype physical spans differ")
                    chrom_span = span
                    for index, (start, end, hid) in enumerate(intervals):
                        if index and start != intervals[index - 1][1] + 1:
                            raise ValueError("Ancestry segments have a gap or overlap")
                        source = ancestry[hid]
                        total[SOURCES.index(source)] += end - start + 1
                        writer.writerow(dict(zip(fields, [
                            sample, chrom, number, start, end + 1, source,
                        ])))
            spans[str(chrom)] = chrom_span
    expected = collapse(selected[[f"realized_{s}" for s in ORIGINAL_SOURCES]])
    actual = np.vstack([totals[sample] / totals[sample].sum() for sample in selected.sample_id])
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)
    return spans, tracked


def prepare(root, source, out):
    data, training, training_audit, tracked = common(root, source)
    saved, cols, pedigree, pure, validation, real, donors, *_, provenance = data
    inputs = out / "inputs"
    inputs.mkdir(exist_ok=False)
    ref_vcf = source / "source_reference.autosomes.vcf.gz"
    ref_ids = samples(ref_vcf)
    authority = pd.read_csv(source / "ancestry_role_audit.tsv", sep="\t")
    labels = authority.loc[authority["dataset"].eq("reference")
                           & authority["sample"].isin(ref_ids)].copy()
    if (labels["sample"].duplicated().any() or set(labels["sample"]) != set(ref_ids)
            or not labels.role.eq("source").all()
            or not labels.source_group.isin(ORIGINAL_SOURCES).all()):
        raise ValueError("Every actual reference donor needs an unambiguous pure-source label")
    original_panel = pd.read_csv(source / "source_reference.ref.panel", sep=r"\s+",
                                 names=["sample", "fine"])
    check = labels[["sample", "fine"]].merge(original_panel, on="sample", how="outer",
                                            validate="one_to_one", suffixes=("_audit", "_panel"))
    if len(check) != len(ref_ids) or not check.fine_audit.eq(check.fine_panel).all():
        raise ValueError("Reference-panel identities/fine labels differ from authoritative metadata")
    eval_donors = set().union(*donors.values())
    train_ids = set(pure.sample_id)
    if (set(ref_ids) & (train_ids | set(validation.sample_id) | eval_donors | set(real.sample_id))
            or train_ids & (set(validation.sample_id) | eval_donors)):
        raise ValueError("Reference, training and evaluation donors are not disjoint")
    labels = labels.rename(columns={"sample": "sample_id"})
    labels["original_source_group"] = labels.source_group
    labels["source_group"] = labels.source_group.map(broad)
    labels = labels.set_index("sample_id").loc[ref_ids].reset_index()
    reference_labels = table(labels[["sample_id", "source_group", "fine",
                                    "original_source_group"]], inputs / "reference_labels.tsv")
    pure_labels = pure[["sample_id", "source_group", "fine"]].copy()
    pure_labels["original_source_group"] = pure_labels.source_group
    pure_labels["source_group"] = pure_labels.source_group.map(broad)
    pure_labels_path = table(pure_labels, inputs / "pure_train_labels.tsv")
    pure_vcf = source / "source_train.autosomes.vcf.gz"
    if set(samples(pure_vcf)) != train_ids:
        raise ValueError("Pure-training genotype and label coverage differ")
    panel = inputs / "reference_broad.panel"
    labels[["sample_id", "source_group"]].to_csv(panel, sep="\t", index=False, header=False)
    retained = inputs / "reference.prefiltered.vcf.gz"
    run(["bcftools", "view", "-m2", "-M2", "-v", "snps", "-q", ".005:minor",
         "-c", "50:minor", "-Oz", "-o", retained, ref_vcf], inputs, "filter_reference")
    run(["bcftools", "index", "-t", retained], inputs, "index_reference")
    marker_path = inputs / "markers.tsv"
    with marker_path.open("w") as stream:
        subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
                        str(retained)], check=True, stdout=stream)
    markers = pd.read_csv(marker_path, sep="\t", names=["chrom", "pos", "ref", "alt"])
    if (set(markers.chrom.astype(str)) != {str(c) for c in range(1, 23)}
            or markers.duplicated(["chrom", "pos"]).any()):
        raise ValueError("Retained reference requires unique positions on all 22 autosomes")
    cohort = {
        "sources": list(SOURCES), "seed": SEED,
        "genetic_map": "/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map",
        "reference_filters": {"min_maf": .005, "min_mac": 50},
        "markers_tsv": str(marker_path),
        "reference": {"vcf": str(retained), "labels_tsv": reference_labels,
                      "fine_panel": str(source / "source_reference.ref.panel"),
                      "broad_panel": str(panel)},
        "pure_train": {"vcf": str(pure_vcf), "labels_tsv": pure_labels_path},
        "pedigree_train": {"batches": []},
        "evaluation": {"pure": {"batches": []}, "pedigree": {"batches": []},
                       "real": {"batches": []}},
        "chromosome_spans": None,
    }
    all_spans = None
    for kind, directory, frame in (
        ("pedigree_train", root / "results_v3_pedigree_training_screen_01/cohort", training),
        ("pedigree", root / "results_v3_pedigree_pilot_10", pedigree),
    ):
        vcfs = sorted((directory / "pedigree/panel_genomes").glob("*.autosomes.vcf.gz"))
        if not vcfs:
            raise ValueError("Reusable pedigree genomes are missing")
        covered = set()
        for original_vcf in vcfs:
            batch = original_vcf.name.removesuffix(".autosomes.vcf.gz")
            tracked.append(original_vcf)
            vcf = indexed_input(original_vcf, inputs, f"{kind}_{batch}")
            ids = samples(vcf)
            selected = frame[frame.sample_id.isin(ids)].copy()
            if set(selected.sample_id) != set(ids) or covered & set(ids):
                raise ValueError("Pedigree batch must contain exactly the expected descendants")
            covered.update(ids)
            selected["scenario"] = selected["class"]
            sample_path = inputs / f"{kind}_{batch}.samples.tsv"
            segment_path = inputs / f"{kind}_{batch}.segments.tsv"
            table(selected.drop(columns=cols), sample_path)
            spans, paths = export_segments(directory, batch, selected, segment_path)
            tracked.extend(paths)
            if all_spans is not None and spans != all_spans:
                raise ValueError("Training and development truth have different chromosome spans")
            all_spans = spans
            row = {"vcf": str(vcf), "samples_tsv": str(sample_path),
                   "segments_tsv": str(segment_path), "batch": batch}
            target = (cohort["pedigree_train"] if kind == "pedigree_train"
                      else cohort["evaluation"]["pedigree"])
            target["batches"].append(row)
        if covered != set(frame.sample_id):
            raise ValueError(f"Incomplete {kind} genotype coverage")
    cohort["chromosome_spans"] = all_spans
    for kind, frame, vcf in (
        ("pure", validation, source / "source_validation.autosomes.vcf.gz"),
        ("real", real, source / "real_user_validation/lonni/lonni.autosomes.vcf.gz"),
    ):
        tracked.append(vcf)
        vcf = indexed_input(vcf, inputs, kind)
        if set(samples(vcf)) != set(frame.sample_id):
            raise ValueError(f"{kind} evaluation genotype coverage differs")
        metadata = frame.drop(columns=cols).copy()
        if kind == "pure":
            metadata["scenario"] = "reference_source"
            metadata["family_id"] = "ref:" + metadata.sample_id
            for s in ORIGINAL_SOURCES:
                metadata[f"realized_{s}"] = metadata[s]
        sample_path = inputs / f"{kind}.samples.tsv"
        table(metadata, sample_path)
        cohort["evaluation"][kind]["batches"].append({
            "vcf": str(vcf), "samples_tsv": str(sample_path), "batch": kind,
        })
    mask = np.array([all_spans[str(c)][0] <= p < all_spans[str(c)][1]
                     for c, p in zip(markers.chrom, markers.pos)])
    in_span = markers[mask].reset_index(drop=True)
    real_markers = marker_table(cohort["evaluation"]["real"]["batches"][0]["vcf"])
    availability = in_span.merge(real_markers.drop_duplicates(), how="left", indicator=True,
                                 validate="one_to_one")
    common_markers = availability[availability._merge.eq("both")].drop(columns="_merge").reset_index(drop=True)
    missing_real = availability[availability._merge.eq("left_only")].drop(columns="_merge")
    if set(common_markers.chrom) != set(range(1, 23)):
        raise ValueError("Shared real-compatible marker panel must cover all 22 autosomes")
    dropped = markers[~mask]
    intervals = inputs / "common_markers.regions.tsv"
    common_markers[["chrom", "pos", "pos"]].to_csv(intervals, sep="\t", header=False, index=False)
    common_ref = inputs / "reference.retained.vcf.gz"
    run(["bcftools", "view", "-R", intervals, "-Oz", "-o", common_ref, retained],
        inputs, "clip_reference")
    run(["bcftools", "index", "-t", common_ref], inputs, "index_common_reference")
    pd.testing.assert_frame_equal(marker_table(common_ref), common_markers)
    common_markers.to_csv(marker_path, sep="\t", header=False, index=False)
    cohort["reference"]["vcf"] = str(common_ref)
    cohort["pure_train"]["vcf"], pure_audit = canonical_genotypes(
        pure_vcf, common_ref, common_markers, inputs, "pure_train")
    canonicalization = {"pure_train": pure_audit}
    for kind, group in [("pedigree_train", cohort["pedigree_train"]),
                        *cohort["evaluation"].items()]:
        for batch in group["batches"]:
            name = f"{kind}_{batch['batch']}"
            batch["vcf"], canonicalization[name] = canonical_genotypes(
                batch["vcf"], common_ref, common_markers, inputs, name)
    write_json(inputs / "canonicalization_audit.json", {
        "reference_prefiltered_markers": len(markers),
        "reference_common_markers": len(common_markers),
        "dropped_outside_training_truth_span": dropped.to_dict("records"),
        "dropped_absent_real_regression_markers_n": len(missing_real),
        "selection_uses_evaluation_labels_or_allele_frequencies": False,
        "policy": "reference-only frequency filters; training map spans; user-approved real-compatible SNP identities; exact REF/ALT intersection",
        "genotypes": canonicalization,
    })
    merged_training = inputs / "refinement_training.vcf.gz"
    training_vcfs = [Path(cohort["pure_train"]["vcf"]),
                     *[Path(b["vcf"]) for b in cohort["pedigree_train"]["batches"]]]
    run(["bcftools", "merge", "-m", "none", "-Oz", "-o", merged_training,
         *training_vcfs], inputs, "merge_refinement")
    run(["bcftools", "index", "-t", merged_training], inputs, "index_refinement")
    if set(samples(merged_training)) != train_ids | set(training.sample_id):
        raise ValueError("Refinement genotypes have unexpected coverage")
    train_truth = pd.concat([
        pd.DataFrame(collapse(pure[list(ORIGINAL_SOURCES)]), columns=SOURCES,
                     index=pure.sample_id),
        pd.DataFrame(collapse(training[[f"realized_{s}" for s in ORIGINAL_SOURCES]]),
                     columns=SOURCES, index=training.sample_id),
    ]).rename_axis("sample_id").reset_index()
    cohort["refinement"] = {"vcf": str(merged_training),
                            "truth_tsv": table(train_truth, inputs / "refinement_truth.tsv")}
    jar = root / "results_v3_global_precision_01/precision.jar"
    old = json.loads((root / "results_v3_global_precision_01/input_audit.json").read_text())
    if digest(jar) != old["input_sha256"]["/results/precision.jar"]:
        raise ValueError("Verified isolated precision jar changed")
    cohort["flare_jar"] = str(jar)
    write_json(out / "cohort.json", cohort)
    configs = []
    for candidate in MODES:
        config = {"candidate": candidate, "sources": list(SOURCES), "seed": SEED}
        if candidate in FLARE_MODES:
            _, structure, policy = candidate.split("_")
            config.update(family="FLARE", structure=structure, refinement=policy,
                          components_per_source=2 if structure == "subclusters" else 1,
                          generation=10, uniform_broad_prior=True,
                          calibration=False, decoder="sum source-anchored components")
        else:
            _, architecture, window = candidate.split("_")
            config.update(family="Gnomix", architecture=architecture,
                          window_size_cM=.2 if window == "020" else .5,
                          smooth_size=75 if window == "020" else 31,
                          context_ratio=.5, retrain_base=False, calibrate=False)
        configs.append(config)
    write_json(out / "configurations.json", configs)
    baseline = out / "baseline_cached"
    baseline.mkdir()
    baseline_rows = []
    for kind, frame in (("pure", validation), ("pedigree", pedigree), ("real", real)):
        prediction = collapse(frozen_prediction(saved, frame[cols].to_numpy(float)))
        result = pd.DataFrame(prediction, columns=[f"pred_{s}" for s in SOURCES])
        result.insert(0, "sample_id", frame.sample_id.to_numpy())
        result.insert(0, "dataset", kind)
        baseline_rows.append(result)
    pd.concat(baseline_rows).to_csv(baseline / "predictions.tsv", sep="\t", index=False)
    tracked.extend([ref_vcf, source / "source_reference.ref.panel",
                    source / "ancestry_role_audit.tsv", pure_vcf, jar,
                    Path(cohort["genetic_map"]), Path(__file__),
                    root / "config/unified_selection.json", root / "config/real_benchmarks.tsv"])
    tracked.extend(Path(b["vcf"]) for group in [cohort["pedigree_train"],
                   *cohort["evaluation"].values()] for b in group["batches"])
    tracked.extend(inputs.glob("*"))
    tracked.extend([out / "cohort.json", out / "configurations.json"])
    provenance["input_sha256"].update({str(p): digest(p) for p in tracked})
    provenance.update({
        "sources": list(SOURCES), "n_candidates": 12, "n_reference_base": len(ref_ids),
        "n_retained_markers": len(common_markers), "training_n_descendants": len(training),
        "training_n_families": training_audit["n_families"],
        "uses_old_mosaic_training_examples": False,
        "evaluation": "previously inspected individual/founder-disjoint development cohorts",
        "population_disjoint": False, "gnomix_use": "user-confirmed academic noncommercial research",
        "primary_aggregation": "marker-weighted posterior mean, matching frozen FLARE baseline",
        "secondary_aggregation": "physical-midpoint-weighted posterior mean (no cached baseline)",
        "pgs_calibration_changed": False,
    })
    write_json(out / "input_audit.json", provenance)
    score(out, "baseline_cached")
    print(f"Prepared 12 candidates: {len(ref_ids)} reference, {len(pure)} pure training, "
          f"{len(training)} pedigree training, {len(validation)} pure validation, "
          f"{len(pedigree)} development descendants; {len(common_markers)} common markers", flush=True)


def verify(out):
    audit = json.loads((out / "input_audit.json").read_text())
    for path, expected in audit["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Model-family input changed: {path}")
    return json.loads((out / "cohort.json").read_text())


def flare_command(cohort, gt, panel, prefix, model=None, refined=False, priors=None,
                  panel_probs=False, probabilities=False):
    command = [
        "java", "-Xmx32g", "-jar", cohort["flare_jar"],
        f"gt={gt}", f"ref={cohort['reference']['vcf']}", f"ref-panel={panel}",
        f"map={cohort['genetic_map']}", f"out={prefix}", "array=false",
        "min-maf=0.005", "min-mac=50", "nthreads=6", f"seed={SEED}",
        f"em={str(refined).lower()}", f"update-p={str(refined).lower()}",
        f"probs={str(probabilities).lower()}",
    ]
    if model is not None:
        command.append(f"model={model}")
    if priors is not None:
        command.append(f"gt-ancestries={priors}")
    if panel_probs:
        command.extend(["panel-probs=true", "panel-cm=0.5"])
    return command


def read_model(path):
    lines = [line.split() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) < 7:
        raise ValueError("Incomplete FLARE model")
    names, panels = lines[:2]
    n = len(names)
    if len(lines) != 2 * n + 5:
        raise ValueError("FLARE model line count differs from ancestry count")
    model = {
        "ancestries": names, "panels": panels, "generation": float(lines[2][0]),
        "mu": np.array(lines[3], float),
        "p": np.array(lines[4:4 + n], float),
        "theta": np.array(lines[4 + n:4 + 2 * n], float),
        "rho": np.array(lines[-1], float),
    }
    fractions(model["p"], "Copying matrix")
    if ((model["p"] < 0).any() or (model["p"] > 1).any()
            or model["theta"].shape != model["p"].shape
            or not np.isfinite(model["theta"]).all()
            or np.any((model["theta"] < 0) | (model["theta"] > 1))
            or len(model["rho"]) != n or not np.isfinite(model["rho"]).all()
            or (model["rho"] <= 0).any()):
        raise ValueError("Invalid FLARE theta/rho parameters")
    return model


def write_model(path, model):
    n = len(model["ancestries"])
    model = dict(model, generation=10., mu=np.full(n, 1. / n))
    fractions(model["p"], "Copying matrix")
    if (model["p"] < 0).any() or (model["p"] > 1).any():
        raise ValueError("FLARE model probabilities require strict [0,1] bounds")
    lines = ["# Seven broad sources; equal broad priors; training-only copying parameters.",
             "\t".join(model["ancestries"]), "\t".join(model["panels"]),
             str(model["generation"]), "\t".join(map(str, model["mu"]))]
    for matrix in ("p", "theta"):
        lines.extend("\t".join(map(str, row)) for row in model[matrix])
    lines.append("\t".join(map(str, model["rho"])))
    Path(path).write_text("\n\n".join(lines) + "\n")
    restored = read_model(path)
    np.testing.assert_allclose(restored["p"], model["p"], rtol=0, atol=0)


def copying_profiles(path, pure_ids, pure_sources, reservoir_size=128):
    rng = np.random.default_rng(SEED)
    with Path(path).open() as stream:
        nrefhaps = int(next(stream).strip().split("=")[1])
        panels = next(stream).split()[4:]
        nhaps, npanels = 2 * len(pure_ids), len(panels)
        counts = np.zeros(nhaps, dtype=int)
        reservoirs = np.zeros((nhaps, reservoir_size, npanels + 1), float)
        sums = np.zeros((nhaps, npanels + 1), float)
        old_location, seen = None, set()
        for raw in stream:
            if raw.startswith("#") or not raw.strip():
                continue
            bits = raw.split(maxsplit=4)
            location = (bits[0], bits[1].split(":")[0])
            hap = int(bits[2])
            if location != old_location:
                if old_location is not None and seen != set(range(nhaps)):
                    raise ValueError("Incomplete copying-profile haplotype coverage")
                old_location, seen = location, set()
            if hap not in range(nhaps) or hap in seen:
                raise ValueError("Duplicate/out-of-range copying-profile haplotype")
            seen.add(hap)
            probabilities = np.fromstring(bits[4], sep=" ")
            rho = float(bits[3])
            if (len(probabilities) != npanels or not np.isfinite(probabilities).all()
                    or (probabilities < 0).any() or probabilities.sum() <= 0
                    or not np.isfinite(rho) or rho <= 0):
                raise ValueError("Invalid copying-profile probabilities/rho")
            values = np.r_[probabilities / probabilities.sum(), rho]
            sums[hap] += values
            counts[hap] += 1
            slot = counts[hap] - 1 if counts[hap] <= reservoir_size else rng.integers(counts[hap])
            if slot < reservoir_size:
                reservoirs[hap, slot] = values
        if seen != set(range(nhaps)) or counts.min() <= 0 or counts.min() != counts.max():
            raise ValueError("Copying profiles require equal complete coverage of all training haplotypes")
    source_indices = np.repeat([SOURCES.index(pure_sources[s]) for s in pure_ids], 2)
    if set(source_indices) != set(range(len(SOURCES))):
        raise ValueError("Copying-profile training donors must cover all seven sources")
    means = np.vstack([sums[source_indices == i].sum(axis=0)
                       / counts[source_indices == i].sum() for i in range(len(SOURCES))])
    retained = min(reservoir_size, counts.min())
    return nrefhaps, panels, means, reservoirs[:, :retained], source_indices


def empirical_centers(observations, labels, components):
    if (observations.ndim != 2 or len(observations) != len(labels)
            or not np.isfinite(observations).all() or (observations < 0).any()
            or set(labels) != set(range(components))):
        raise ValueError("Invalid observations or empty assigned clusters")
    return np.vstack([observations[labels == i].mean(axis=0) for i in range(components)])


def initialize(out, reuse_profiles=False):
    cohort = verify(out)
    directory = out / "initialization"
    directory.mkdir(exist_ok=reuse_profiles)
    ref = cohort["reference"]
    if not reuse_profiles:
        run(flare_command(cohort, cohort["pure_train"]["vcf"], ref["broad_panel"],
                          directory / "direct"), directory, "direct")
        model = read_model(directory / "direct.model")
        if set(model["ancestries"]) != set(SOURCES):
            raise ValueError("Native broad FLARE initialization lacks one of seven sources")
        write_model(directory / "direct.initial.model", model)
        run(flare_command(cohort, cohort["pure_train"]["vcf"], ref["fine_panel"],
                          directory / "profiles", panel_probs=True), directory, "profiles")
    pure = pd.read_csv(cohort["pure_train"]["labels_tsv"], sep="\t")
    ids = samples(cohort["pure_train"]["vcf"])
    nrefhaps, panels, means, reservoir, indices = copying_profiles(
        directory / "profiles.panels", ids, pure.set_index("sample_id").source_group.to_dict())
    lam = 1 / (np.log(nrefhaps) + .5)
    theta = lam / (2 * lam + 2 * nrefhaps)
    for structure, components in (("fine", 1), ("subclusters", 2)):
        names, parameters, rhos = [], [], []
        for i, source in enumerate(SOURCES):
            if components == 1:
                centers, rates = means[i:i + 1, :-1], means[i:i + 1, -1]
            else:
                observations = reservoir[indices == i].reshape(-1, len(panels) + 1)
                cluster = KMeans(n_clusters=2, n_init=10, random_state=SEED,
                                 max_iter=300).fit(observations[:, :-1])
                labels = cluster.labels_
                if set(labels) != {0, 1}:
                    raise ValueError(f"Within-source clustering collapsed: {source}")
                assigned = empirical_centers(observations, labels, 2)
                centers, rates = assigned[:, :-1], assigned[:, -1]
            for j in range(components):
                names.append(source if components == 1 else f"{source}__{j + 1}")
                parameters.append(centers[j] / centers[j].sum())
                rhos.append(rates[j])
        model = {"ancestries": names, "panels": panels, "p": np.vstack(parameters),
                 "theta": np.full((len(names), len(panels)), theta), "rho": np.array(rhos)}
        write_model(directory / f"{structure}.initial.model", model)
    write_json(directory / "initialization_audit.json", {
        "label_source": "pure-training donors only; evaluation excluded",
        "subclusters": "two within-source kmeans components; no cross-source clustering",
        "reservoir_per_haplotype": reservoir.shape[1],
        "n_reference_haplotypes": nrefhaps,
        "input_sha256": {str(p): digest(p) for p in directory.glob("*.initial.model")},
        "copying_profile_sha256": digest(directory / "profiles.panels"),
    })


def recover(previous, out):
    audit = json.loads((previous / "revisions/original_input_audit.json").read_text())
    original = previous / "revisions/original_screen_model_families.py"
    if digest(original) != audit["input_sha256"][str(Path(__file__))]:
        raise ValueError("Archived original implementation does not match the old input audit")
    for name in ("inputs", "initialization", "baseline", "baseline_cached"):
        shutil.copytree(previous / name, out / name)
    for name in ("cohort.json", "configurations.json"):
        shutil.copyfile(previous / name, out / name)
    for path, expected in audit["input_sha256"].items():
        actual = Path(path)
        if actual == Path(__file__):
            actual = original
        if digest(actual) != expected:
            raise ValueError(f"Recovery input changed: {path}")
    shutil.copyfile(previous / "revisions/original_input_audit.json",
                    out / "recovered_input_audit.json")
    audit["input_sha256"][str(Path(__file__))] = digest(__file__)
    audit["recovery"] = {
        "source": str(previous), "original_implementation_sha256": digest(original),
        "original_input_audit_sha256": digest(previous / "revisions/original_input_audit.json"),
        "change": "strict model bounds; empirical assigned-cluster means avoid tiny negative kmeans centers",
        "prepared_genotypes_truth_and_baseline_unchanged": True,
    }
    write_json(out / "input_audit.json", audit)
    old = json.loads((out / "initialization/initialization_audit.json").read_text())
    if digest(out / "initialization/profiles.panels") != old["copying_profile_sha256"]:
        raise ValueError("Recovered raw copying profiles differ")
    for source in ("direct", "fine"):
        path = out / f"initialization/{source}.initial.model"
        if digest(path) != old["input_sha256"][str(path)]:
            raise ValueError("Recovered valid initial model changed")
    initialize(out, reuse_profiles=True)
    for source in ("direct", "fine"):
        path = out / f"initialization/{source}.initial.model"
        if digest(path) != old["input_sha256"][str(path)]:
            raise ValueError("Centroid repair altered a previously valid model")
    reused = []
    for candidate in ("flare_direct_fixed", "flare_direct_refined",
                      "flare_fine_fixed", "flare_fine_refined", "gnomix_fast_020"):
        source = previous / f"smoke/{candidate}"
        filename = "output_audit.json" if candidate.startswith("gnomix") else "persistence_audit.json"
        check = json.loads((source / filename).read_text())
        if check.get("model_roundtrip_verified") is not True or check.get("smoke_only") is not True:
            raise ValueError("Recovery requires a complete, roundtrip-verified smoke")
        if digest(source / "predictions.tsv") != check["prediction_sha256"]:
            raise ValueError("Recovered smoke prediction checksum differs")
        shutil.copytree(source, out / f"smoke/{candidate}")
        reused.append(candidate)
    write_json(out / "recovery_audit.json", {
        "all_previous_inputs_verified": True, "original_source": str(previous),
        "reused_smokes": reused, "reused_frozen_baseline": True,
        "reused_raw_copying_profiles": True,
        "new_subcluster_model_sha256": digest(out / "initialization/subclusters.initial.model"),
    })


def component_source(name):
    source = name.split("__")[0]
    if source not in SOURCES:
        raise ValueError(f"Unanchored FLARE component: {name}")
    return source


def ancestry_priors(truth, model, output):
    frame = pd.read_csv(truth, sep="\t").set_index("sample_id")
    fractions(frame[list(SOURCES)], "Refinement truth")
    names = model["ancestries"]
    counts = {s: sum(component_source(n) == s for n in names) for s in SOURCES}
    if min(counts.values()) <= 0 or len(set(counts.values())) != 1:
        raise ValueError("Source components must have equal complete broad coverage")
    priors = pd.DataFrame({
        n: frame[component_source(n)] / counts[component_source(n)] for n in names
    }, index=frame.index)
    fractions(priors, "Source-constrained training priors")
    priors.rename_axis("SAMPLE").to_csv(output, sep="\t")


def midpoint_weights(positions, span):
    positions = np.asarray(positions, float)
    if (not len(positions) or np.any(np.diff(positions) <= 0)
            or span[0] > positions[0] or span[1] <= positions[-1]):
        raise ValueError("Markers must be unique, ordered and inside the physical truth span")
    edges = np.r_[span[0], (positions[:-1] + positions[1:]) / 2, span[1]]
    return np.diff(edges)


def physical_predictions(path, names, spans, expected_markers):
    sample_ids, sums, marker_index = None, None, 0
    chromosome = None
    previous = None
    left_edge = None
    total_span = 0.
    component_indices = [SOURCES.index(component_source(n)) for n in names]
    with gzip.open(path, "rt") as stream:
        for raw in stream:
            if raw.startswith("##"):
                continue
            fields = raw.rstrip().split("\t")
            if raw.startswith("#CHROM"):
                sample_ids = fields[9:]
                sums = np.zeros((len(sample_ids), len(SOURCES)))
                continue
            chrom, pos = fields[0].removeprefix("chr"), int(fields[1])
            if marker_index >= len(expected_markers):
                raise ValueError("Local ancestry output contains unexpected extra markers")
            expected = expected_markers.iloc[marker_index]
            if (chrom, pos, fields[3], fields[4]) != (
                    str(expected.chrom), int(expected.pos), expected.ref, expected.alt):
                raise ValueError("FLARE retained marker identities differ from the shared reference")
            marker_index += 1
            format_fields = fields[8].split(":")
            first, second = format_fields.index("ANP1"), format_fields.index("ANP2")
            vectors = np.fromstring(",".join(
                value.split(":")[index] for value in fields[9:] for index in (first, second)
            ), sep=",").reshape(len(sample_ids), 2, len(names))
            if (not np.isfinite(vectors).all() or (vectors < 0).any()
                    or np.any(vectors.sum(axis=2) <= 0)):
                raise ValueError("Invalid local ancestry posteriors")
            vectors /= vectors.sum(axis=2, keepdims=True)
            collapsed = np.zeros_like(sums)
            for index, source_index in enumerate(component_indices):
                collapsed[:, source_index] += vectors[:, :, index].mean(axis=1)
            if chromosome != chrom:
                if previous is not None:
                    end = spans[chromosome][1]
                    sums += previous[1] * (end - left_edge)
                    total_span += end - spans[chromosome][0]
                chromosome, previous, left_edge = chrom, (pos, collapsed), spans[chrom][0]
                if pos < left_edge:
                    raise ValueError("Posterior marker is outside ancestry truth span")
            else:
                if pos <= previous[0]:
                    raise ValueError("Unordered local ancestry markers")
                boundary = (previous[0] + pos) / 2
                sums += previous[1] * (boundary - left_edge)
                previous, left_edge = (pos, collapsed), boundary
    if previous is None or marker_index != len(expected_markers):
        raise ValueError("Incomplete FLARE local ancestry marker coverage")
    end = spans[chromosome][1]
    if previous[0] >= end:
        raise ValueError("Posterior marker is outside ancestry truth span")
    sums += previous[1] * (end - left_edge)
    total_span += end - spans[chromosome][0]
    return sample_ids, fractions(sums / total_span, "Physical posterior mean")


def subset_chromosome(cohort, directory, chrom):
    result = json.loads(json.dumps(cohort))
    paths = {cohort["reference"]["vcf"], cohort["pure_train"]["vcf"],
             cohort["refinement"]["vcf"]}
    paths.update(b["vcf"] for group in [cohort["pedigree_train"],
                 *cohort["evaluation"].values()] for b in group["batches"])
    lookup = {}
    for index, original in enumerate(sorted(paths)):
        target = directory / f"chr{chrom}_{index:02d}.vcf.gz"
        run(["bcftools", "view", "-r", str(chrom), "-Oz", "-o", target, original],
            directory, f"subset_{index:02d}")
        run(["bcftools", "index", "-t", target], directory, f"index_{index:02d}")
        lookup[original] = str(target)
    for group in ("reference", "pure_train", "refinement"):
        result[group]["vcf"] = lookup[result[group]["vcf"]]
    for group in [result["pedigree_train"], *result["evaluation"].values()]:
        for batch in group["batches"]:
            batch["vcf"] = lookup[batch["vcf"]]
    return result


def flare_candidate(out, candidate, smoke=False):
    cohort = verify(out)
    if candidate not in FLARE_MODES:
        raise ValueError("Unknown FLARE candidate")
    directory = out / ("smoke" if smoke else "candidates") / candidate
    directory.mkdir(parents=True, exist_ok=False)
    if smoke:
        cohort = subset_chromosome(cohort, directory, 22)
    _, structure, policy = candidate.split("_")
    original = out / f"initialization/{structure}.initial.model"
    model = read_model(original)
    panel = cohort["reference"]["broad_panel" if structure == "direct" else "fine_panel"]
    selected = directory / "model.final.model"
    if policy == "refined":
        priors = directory / "training.priors.tsv"
        ancestry_priors(cohort["refinement"]["truth_tsv"], model, priors)
        run(flare_command(cohort, cohort["refinement"]["vcf"], panel, directory / "training",
                          model=original, refined=True, priors=priors), directory, "training")
        model = read_model(directory / "training.model")
        if model["ancestries"] != read_model(original)["ancestries"]:
            raise ValueError("Refinement changed anchored ancestry identities")
    write_model(selected, model)
    names = model["ancestries"]
    markers = pd.read_csv(cohort["markers_tsv"], sep="\t",
                          names=["chrom", "pos", "ref", "alt"])
    if smoke:
        markers = markers[markers.chrom.eq(22)].reset_index(drop=True)
    primary, physical = [], []
    for dataset, group in cohort["evaluation"].items():
        for batch in group["batches"]:
            task = f"{dataset}_{batch['batch']}"
            run(flare_command(cohort, batch["vcf"], panel, directory / task,
                              model=selected, probabilities=True), directory, task)
            globals_path = directory / f"{task}.global.anc.gz"
            frame = pd.read_csv(globals_path, sep=r"\s+").set_index("SAMPLE")
            if list(frame.columns) != names:
                raise ValueError("FLARE global columns differ from anchored model")
            values = np.zeros((len(frame), len(SOURCES)))
            for n in names:
                values[:, SOURCES.index(component_source(n))] += frame[n].to_numpy(float)
            fractions(values, "Marker posterior mean")
            expected_ids = samples(batch["vcf"])
            if set(frame.index) != set(expected_ids):
                raise ValueError("FLARE prediction sample coverage differs from evaluation")
            row = pd.DataFrame(values, columns=[f"pred_{s}" for s in SOURCES], index=frame.index)
            row = row.rename_axis("sample_id").reset_index()
            row.insert(0, "dataset", dataset)
            primary.append(row)
            ids, values = physical_predictions(directory / f"{task}.anc.vcf.gz", names,
                                                cohort["chromosome_spans"], markers)
            row = pd.DataFrame(values, columns=[f"pred_{s}" for s in SOURCES])
            row.insert(0, "sample_id", ids)
            row.insert(0, "dataset", dataset)
            physical.append(row)
    pd.concat(primary).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    pd.concat(physical).to_csv(directory / "predictions_physical.tsv", sep="\t", index=False)
    write_json(directory / "persistence_audit.json", {
        "model_roundtrip_verified": True, "smoke_only": smoke,
        "ancestry_identity": "source-constrained pure and pedigree training priors",
        "refinement_parameters_retained": "p, theta, rho; evaluation T=10, broad priors uniform",
        "prediction_sha256": digest(directory / "predictions.tsv"),
        "physical_prediction_sha256": digest(directory / "predictions_physical.tsv"),
        "model_sha256": digest(selected),
    })
    if not smoke:
        score(out, candidate)


def matched_baseline(root, source, out):
    cohort = verify(out)
    saved, cols, *_ = common(root, source)[0]
    directory = out / "baseline"
    directory.mkdir(exist_ok=False)
    model = root / "results_v3_sweep_100/models/v1_frozen_k18.final.model"
    frames = []
    for dataset, group in cohort["evaluation"].items():
        for batch in group["batches"]:
            task = f"{dataset}_{batch['batch']}"
            run(flare_command(cohort, batch["vcf"], cohort["reference"]["fine_panel"],
                              directory / task, model=model), directory, task)
            features = pd.read_csv(directory / f"{task}.global.anc.gz",
                                   sep=r"\s+").set_index("SAMPLE")
            if set(features.columns) != set(cols):
                raise ValueError("Matched frozen-v1 feature identities changed")
            x = np.rint(features[cols].to_numpy(float) * 1000) / 1000
            if (x.sum(axis=1) <= 0).any():
                raise ValueError("Frozen feature quantization collapsed a sample")
            x /= x.sum(axis=1, keepdims=True)
            prediction = collapse(frozen_prediction(saved, x))
            frame = pd.DataFrame(prediction, columns=[f"pred_{s}" for s in SOURCES])
            frame.insert(0, "sample_id", features.index)
            frame.insert(0, "dataset", dataset)
            frames.append(frame)
    pd.concat(frames).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    write_json(directory / "persistence_audit.json", {
        "frozen_models_retrained": False,
        "input_control": "same canonical marker set as twelve new candidates",
        "decoder_input": "historical three-decimal globals, normalized",
        "prediction_sha256": digest(directory / "predictions.tsv"),
    })
    score(out, "baseline")


def unified_seven(frame, cfg):
    p = fractions(frame[[f"pred_{s}" for s in SOURCES]], "Combined prediction")
    y = fractions(frame[[f"truth_{s}" for s in SOURCES]], "Combined truth")
    d = frame.copy()
    d["tv"] = np.abs(p - y).sum(axis=1) / 2
    d["off"] = (p * (y < cfg["absent_ancestry_max"])).sum(axis=1)
    d["combination"] = ["+".join(s for s, v in zip(SOURCES, row)
                                  if v >= cfg["absent_ancestry_max"]) for row in y]
    scenarios = []
    for name, group in d.groupby("scenario"):
        weights = pd.Series(0., index=group.index)
        ncomb = group.combination.nunique()
        for _, combination in group.groupby("combination"):
            for _, family in combination.groupby("family_id"):
                weights.loc[family.index] = 1 / (ncomb * combination.family_id.nunique() * len(family))
        scenarios.append({
            "scenario": name, "n": len(group),
            "mean_tv": float(np.dot(group.tv, weights.loc[group.index])),
            "p95_tv": weighted_quantile(group.tv.to_numpy(), weights.loc[group.index].to_numpy(), .95),
            "off_target": float(np.dot(group.off, weights.loc[group.index])),
        })
    minor = (y >= cfg["minor_min"]) & (y <= cfg["minor_max"])
    return {
        "n_reference": int(d.scenario.eq("reference_source").sum()),
        "n_pedigree": int(d.scenario.ne("reference_source").sum()),
        "balanced_mean_tv": float(np.mean([r["mean_tv"] for r in scenarios])),
        "balanced_p95_tv": float(np.mean([r["p95_tv"] for r in scenarios])),
        "balanced_off_target_mass": float(np.mean([r["off_target"] for r in scenarios])),
        "minor_2_5_n": int(minor.sum()),
        "minor_2_5_recall": float((p[minor] >= cfg["minor_detection_threshold"]).mean())
        if minor.any() else None, "scenarios": scenarios,
    }


def score(out, candidate):
    cohort = json.loads((out / "cohort.json").read_text())
    directory = out / (candidate if candidate.startswith("baseline") else f"candidates/{candidate}")
    predictions = pd.read_csv(directory / "predictions.tsv", sep="\t")
    if predictions.duplicated(["dataset", "sample_id"]).any():
        raise ValueError("Duplicate candidate predictions")
    reports, combined = {}, []
    for dataset, group in cohort["evaluation"].items():
        metadata = pd.concat([pd.read_csv(b["samples_tsv"], sep="\t")
                              for b in group["batches"]], ignore_index=True)
        pred = predictions[predictions.dataset.eq(dataset)].drop(columns="dataset")
        merged = metadata.merge(pred, on="sample_id", how="outer", validate="one_to_one",
                                indicator=True)
        if not merged._merge.eq("both").all():
            raise ValueError("Candidate has incomplete/unexpected evaluation sample coverage")
        p = fractions(merged[[f"pred_{s}" for s in SOURCES]], "Evaluation prediction")
        if dataset == "real":
            reports["real"] = merged[["sample_id", *[f"pred_{s}" for s in SOURCES]]].to_dict("records")
            continue
        y = collapse(merged[[f"realized_{s}" for s in ORIGINAL_SOURCES]])
        merged["TV_distance"] = np.abs(p - y).sum(axis=1) / 2
        merged["off_target_mass"] = (p * (y < .01)).sum(axis=1)
        merged["dominant_correct"] = p.argmax(axis=1) == y.argmax(axis=1)
        reports[dataset] = {
            "n": len(merged), "mean_tv": float(merged.TV_distance.mean()),
            "p95_tv": float(merged.TV_distance.quantile(.95)),
            "off_target_mass": float(merged.off_target_mass.mean()),
            "dominant_accuracy": float(merged.dominant_correct.mean()),
        }
        if dataset == "pure":
            reports["pure"]["by_original_source"] = {
                str(s): {"n": len(g), "mean_tv": float(g.TV_distance.mean())}
                for s, g in merged.groupby("source_group")
            }
            reports["pure"]["by_population"] = {
                str(s): {"n": len(g), "mean_tv": float(g.TV_distance.mean())}
                for s, g in merged.groupby("fine")
            }
        else:
            reports["pedigree"]["family_macro_tv"] = float(
                merged.groupby("family_id").TV_distance.mean().mean())
        for i, source in enumerate(SOURCES):
            merged[f"truth_{source}"] = y[:, i]
        combined.append(merged[["sample_id", "family_id", "scenario",
                                *[f"pred_{s}" for s in SOURCES], *[f"truth_{s}" for s in SOURCES]]])
        merged.drop(columns="_merge").to_csv(directory / f"{dataset}_predictions_scored.tsv",
                                              sep="\t", index=False)
    if len(predictions) != sum(len(pd.read_csv(b["samples_tsv"], sep="\t"))
                               for g in cohort["evaluation"].values() for b in g["batches"]):
        raise ValueError("Unexpected prediction dataset")
    cfg = json.loads(Path("/workspace/config/unified_selection.json").read_text())
    unified = unified_seven(pd.concat(combined, ignore_index=True), cfg)
    if (unified["n_reference"] < cfg["min_reference_samples"]
            or unified["n_pedigree"] < cfg["min_pedigree_samples"]
            or not set(cfg["required_scenarios"]).issubset(r["scenario"] for r in unified["scenarios"])):
        raise ValueError("Evaluation does not cover unchanged required scenarios")
    gates = {k: unified[k] <= limit for k, limit in cfg["max_gates"].items()}
    gates.update({k: unified[k] >= limit for k, limit in cfg["min_gates"].items()})
    if unified["minor_2_5_n"] >= cfg["minor_gate_min_support"]:
        gates["minor_2_5_recall"] = unified["minor_2_5_recall"] >= cfg["minor_recall_min"]
    real = reports["real"]
    if len(real) != 1:
        raise ValueError("Expected the one frozen real-user regression case")
    real_config = pd.read_csv("/workspace/config/real_benchmarks.tsv", sep="\t")
    if len(real_config) != 1:
        raise ValueError("This fixed screen requires exactly the audited real regression benchmark")
    constraints = json.loads(real_config.iloc[0].constraints_json)
    converted = {}
    for expr, rule in constraints.items():
        parts = expr.split("+")
        if ("EAS" in parts) != ("SEA" in parts):
            raise ValueError("A separate EAS or SEA gate cannot be evaluated with merged outputs")
        if "EAS" in parts:
            parts = [p for p in parts if p not in ("EAS", "SEA")] + ["EAS_SEA"]
        if not set(parts).issubset(SOURCES):
            raise ValueError("Unknown real regression ancestry group")
        converted["+".join(parts)] = rule
    real_pass, real_details = evaluate_constraints(
        {s: real[0][f"pred_{s}"] for s in SOURCES}, converted)
    if cfg["require_real_hard_gate"]:
        gates["real_world_regression"] = real_pass
    reports["real_constraints"] = real_details
    baseline = reports if candidate.startswith("baseline") else json.loads(
        (out / "baseline/metrics.json").read_text())
    gates["pure_no_baseline_regression"] = reports["pure"]["mean_tv"] <= baseline["pure"]["mean_tv"] + 1e-10
    gates["pedigree_no_baseline_regression"] = reports["pedigree"]["mean_tv"] <= baseline["pedigree"]["mean_tv"] + 1e-10
    reports.update(candidate=candidate, sources=list(SOURCES), unified=unified, gates=gates,
                   all_gates_pass=all(gates.values()), development_screen_only=True,
                   minor_gate_evaluable=unified["minor_2_5_n"] >= cfg["minor_gate_min_support"])
    write_json(directory / "metrics.json", reports)
    print(f"{candidate}: pure TV={reports['pure']['mean_tv']:.6f}, "
          f"pedigree TV={reports['pedigree']['mean_tv']:.6f}, "
          f"gates={reports['all_gates_pass']}", flush=True)


def summarize(out):
    rows = []
    for candidate in ("baseline_cached", "baseline", *MODES):
        directory = out / (candidate if candidate.startswith("baseline") else f"candidates/{candidate}")
        if not (directory / "metrics.json").is_file():
            raise ValueError(f"Candidate not completed: {candidate}")
        metrics = json.loads((directory / "metrics.json").read_text())
        rows.append({
            "candidate": candidate, "pure_mean_tv": metrics["pure"]["mean_tv"],
            "pedigree_mean_tv": metrics["pedigree"]["mean_tv"],
            "balanced_mean_tv": metrics["unified"]["balanced_mean_tv"],
            "balanced_p95_tv": metrics["unified"]["balanced_p95_tv"],
            "balanced_off_target_mass": metrics["unified"]["balanced_off_target_mass"],
            "all_gates_pass": metrics["all_gates_pass"],
            "failed_gates": ",".join(k for k, passed in metrics["gates"].items() if not passed),
        })
    pd.DataFrame(rows).to_csv(out / "model_family_screen.tsv", sep="\t", index=False)
    write_json(out / "completion.json", {
        "n_candidates": 12, "n_baselines": 2, "production_promotion": False,
        "all_candidates_evaluated": True, "task": "coarser seven-source development screen",
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True,
                        choices=["prepare", "recover", "initialize", "matched-baseline", "flare", "score", "summarize"])
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--root", type=Path, default=Path("/workspace"))
    parser.add_argument("--source", type=Path, default=Path("/data/flare2_test/source_model"))
    parser.add_argument("--candidate", choices=MODES)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--previous", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.root, args.source, args.outdir)
    elif args.stage == "recover":
        if args.previous is None:
            raise ValueError("Recovery requires an explicitly archived previous run")
        recover(args.previous, args.outdir)
    elif args.stage == "initialize":
        initialize(args.outdir)
    elif args.stage == "matched-baseline":
        matched_baseline(args.root, args.source, args.outdir)
    elif args.stage == "flare":
        flare_candidate(args.outdir, args.candidate, args.smoke)
    elif args.stage == "score":
        verify(args.outdir)
        score(args.outdir, args.candidate)
    else:
        verify(args.outdir)
        summarize(args.outdir)


if __name__ == "__main__":
    main()
