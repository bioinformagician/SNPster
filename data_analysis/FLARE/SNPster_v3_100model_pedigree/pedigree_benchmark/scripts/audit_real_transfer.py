#!/usr/bin/env python3
"""Audit frozen Gnomix transfer without fitting or changing benchmark inputs."""

import argparse
import json
import logging
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from screen_model_families import verify


def query(path, template, regions=None):
    command = ["bcftools", "query", "-f", template]
    if regions is not None:
        command += ["-R", str(regions)]
    command.append(str(path))
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return [line.split("\t") for line in result.stdout.splitlines()]


def randomize_phase(x, seed):
    if x.ndim != 2 or len(x) % 2 or np.any(x > 1) or np.any(x < 0):
        raise ValueError("Phase control requires complete paired biallelic haplotypes")
    paired = x.reshape(-1, 2, x.shape[1])
    swap = np.random.default_rng(seed).integers(0, 2, (len(paired), x.shape[1])).astype(bool)
    result = np.where(swap[:, None, :], paired[:, ::-1], paired).reshape(x.shape)
    if not np.array_equal(result.reshape(paired.shape).sum(axis=1), paired.sum(axis=1)):
        raise ValueError("Phase control changed diploid genotypes")
    return result


def genotype_qc(x):
    paired = x.reshape(-1, 2, x.shape[1])
    complete = (paired < 2).all(axis=1)
    heterozygous = complete & (paired[:, 0] != paired[:, 1])
    return {
        "sites": np.full(len(paired), x.shape[1], dtype=np.int64),
        "missing_sites": (~complete).sum(axis=1),
        "heterozygous_sites": heterozygous.sum(axis=1),
        "heterozygous_01": (heterozygous & (paired[:, 0] == 0)).sum(axis=1),
        "heterozygous_10": (heterozygous & (paired[:, 0] == 1)).sum(axis=1),
        "alt_alleles": ((paired == 1) & complete[:, None]).sum(axis=(1, 2)),
    }


def genotype_posterior_qc(gp):
    if (gp.ndim != 2 or gp.shape[1] != 3 or not np.all(np.isfinite(gp))
            or np.any(gp < 0) or np.any(gp > 1)):
        raise ValueError("Invalid genotype posterior probabilities")
    sums = gp.sum(axis=1)
    # Beagle's two-decimal GP export can round three probabilities by 0.005 each.
    if np.any(np.abs(sums - 1) > .015 + 1e-8):
        raise ValueError("Genotype posterior sums exceed export-rounding tolerance")
    return sums


def thin_sites(keys, eligible, spacing=100_000):
    if len(keys) != len(eligible):
        raise ValueError("Thinning mask differs from marker axis")
    result = np.zeros(len(keys), dtype=bool)
    previous_bin = None
    for index, (key, selected) in enumerate(zip(keys, eligible)):
        current = (key[0], key[1] // spacing)
        if selected and current != previous_bin:
            result[index] = True
            previous_bin = current
    return result


def crf_in_blocks(model, base, size):
    if size < 1:
        raise ValueError("CRF block diagnostic requires a positive block size")
    order = np.argsort([int(value) for value in model.smooth.model.classes_])
    result = np.concatenate(
        [model.smooth.model.predict_proba(base[:, start:start + size])[:, :, order]
         for start in range(0, base.shape[1], size)], axis=1)
    if (result.shape != base.shape or not np.all(np.isfinite(result))
            or np.any(result < 0) or not np.allclose(result.sum(axis=2), 1, atol=2e-6)):
        raise ValueError("Invalid block diagnostic posterior")
    return result / result.sum(axis=2, keepdims=True)


def frequency_fits(frequencies, dosages):
    if (frequencies.ndim != 2 or dosages.ndim != 2
            or dosages.shape[1] != len(frequencies)
            or frequencies.shape[1] != len(adapter.SOURCES)
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies < 0) or np.any(frequencies > 1)
            or np.any(dosages < 0) or np.any(dosages > 2)):
        raise ValueError("Invalid inputs to phase-independent frequency diagnostic")
    informative = np.ptp(frequencies, axis=1) >= .03
    f = frequencies[informative]
    if len(f) < 100:
        raise ValueError("Insufficient informative markers for frequency diagnostic")
    gram = f.T @ f / len(f)
    target = (dosages[:, informative] / 2) @ f / len(f)
    result = []
    for vector in target:
        fit = minimize(
            lambda q: q @ gram @ q - 2 * vector @ q,
            np.full(len(adapter.SOURCES), 1 / len(adapter.SOURCES)),
            jac=lambda q: 2 * (gram @ q - vector),
            constraints={"type": "eq", "fun": lambda q: q.sum() - 1,
                         "jac": lambda q: np.ones_like(q)},
            bounds=[(0, 1)] * len(adapter.SOURCES), method="SLSQP",
            options={"ftol": 1e-12, "maxiter": 1000})
        if not fit.success:
            raise ValueError(f"Frequency diagnostic did not converge: {fit.message}")
        if np.any(fit.x < -1e-10) or not np.isclose(fit.x.sum(), 1, atol=1e-8):
            raise ValueError("Frequency diagnostic returned invalid fractions")
        result.append(fit.x)
    return np.asarray(result), int(informative.sum())


def reference_frequencies(manifest, chrom, groups, output):
    tags = "".join(f"\t%INFO/AF_{source}" for source in adapter.SOURCES)
    command = ["bcftools", "+fill-tags", manifest["reference"]["vcf"],
               "-r", str(chrom), "-Ou", "--", "-t", "AF", "-S", str(groups)]
    with (output / f"reference_frequencies_chr{chrom}.stderr.log").open("w") as errors:
        first = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        try:
            result = subprocess.run(
                ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT" + tags + "\n"],
                stdin=first.stdout, check=True, capture_output=True, text=True)
            first.stdout.close()
            if first.wait():
                raise ValueError(f"Reference frequency calculation failed on chr{chrom}")
        finally:
            if first.poll() is None:
                first.terminate()
                first.wait()
            first.stdout.close()
    rows = [line.split("\t") for line in result.stdout.splitlines()]
    keys = [(r[0], int(r[1]), r[2], r[3]) for r in rows]
    values = np.asarray([[float(value) for value in row[4:]] for row in rows])
    if values.shape != (len(keys), len(adapter.SOURCES)) or not np.all(np.isfinite(values)):
        raise ValueError("Incomplete reference source frequencies")
    return keys, values


def source_rows(manifest):
    reference = adapter.sample_records(manifest["reference"]["labels_tsv"], True)
    training = adapter.sample_records(manifest["pure_train"]["labels_tsv"], True)
    validation = [row for batch in manifest["evaluation"]["pure"]["batches"]
                  for row in adapter.records(batch["samples_tsv"])]
    return pd.DataFrame([dict(role=role, **row) for role, rows in
                         (("reference", reference), ("pure_train", training),
                          ("pure_validation", validation)) for row in rows])


def run(screen, output, upstream, candidate="gnomix_fast_050"):
    if candidate not in ("gnomix_fast_020", "gnomix_fast_050"):
        raise ValueError("Transfer audit currently requires an LR+CRF candidate")
    if (output / "audit.json").exists():
        raise ValueError("Select a new audit output; existing results are immutable")
    output.mkdir(parents=True, exist_ok=True)
    manifest = verify(screen)
    adapter.load_upstream()
    directory = screen / "candidates" / candidate
    original = json.loads((directory / "output_audit.json").read_text())
    if original["status"] != "complete" or not original["model_roundtrip_verified"]:
        raise ValueError("Audit requires a completed, persisted candidate")
    hashes = adapter.input_hashes(manifest)
    if hashes != original["input_sha256"]:
        raise ValueError("Candidate input hashes differ from current inputs")
    audit = {
        "status": "running", "candidate": candidate, "no_training": True,
        "production_changed": False, "acceptance_gates_changed": False,
        "script_sha256": adapter.digest(__file__),
        "adapter_sha256": adapter.digest(adapter.__file__),
        "input_sha256": hashes, "output_audit_sha256": adapter.digest(directory / "output_audit.json"),
        "upstream_sha256": {}, "chromosomes": [], "limitations": [
            "This audit has no independent real-sample ancestry or phase truth.",
            "Phase randomization is a stress control, not measured switch-error ground truth.",
            "Frequency fits are reference-dependent diagnostic projections, not a certified model.",
            "Observed means absence of INFO/IMP; upstream array concordance is not established.",
            "Single-sample DR2 is summarized by genotype; low DR2 alone does not establish bad calls.",
        ],
    }
    adapter.write_json(output / "audit.json", audit)
    populations = source_rows(manifest)
    (populations.groupby(["role", "source_group", "fine"]).size().rename("n").reset_index()
     .to_csv(output / "population_coverage.tsv", sep="\t", index=False))
    reference = populations[populations.role.eq("reference")]
    groups = output / "reference_groups.tsv"
    reference[["sample_id", "source_group"]].to_csv(groups, sep="\t", header=False, index=False)
    batches, all_rows = [], []
    for dataset, group in manifest["evaluation"].items():
        for batch in group["batches"]:
            rows = adapter.records(batch["samples_tsv"])
            batches.append((dataset, batch, rows))
            all_rows.extend(dict(dataset=dataset, **row) for row in rows)
    identities = [(r["dataset"], r["sample_id"]) for r in all_rows]
    if len(set(identities)) != len(identities):
        raise ValueError("Repeated evaluation identities")
    real_index = [i for i, row in enumerate(all_rows) if row["dataset"] == "real"]
    if len(real_index) != 1:
        raise ValueError("Expected exactly one real evaluation sample")
    matched = [i for i, row in enumerate(all_rows)
               if row["dataset"] == "pedigree"
               and abs(float(row["realized_EUR"]) - .5) < 1e-8
               and abs(float(row["realized_EAS"]) + float(row["realized_SEA"]) - .5) < 1e-8]
    phase_indices = np.asarray(real_index + matched, dtype=int)
    metadata = pd.DataFrame(all_rows)
    metadata[["dataset", "sample_id"]].to_csv(output / "sample_axes.tsv", sep="\t", index=False)
    totals, marker_total, chromosome_rows = {}, 0, []
    qc_total = {name: np.zeros(len(all_rows), dtype=np.int64) for name in
                ("sites", "missing_sites", "heterozygous_sites", "heterozygous_01",
                 "heterozygous_10", "alt_alleles")}
    frequency_parts, dosage_parts, observed_parts, thinned_parts, observed_thinned_parts = (
        [], [], [], [], [])
    raw_summary, model_hashes = [], {}
    r = real_index[0]
    try:
        with threadpool_limits(limits=1):
            for chrom in range(1, 23):
                bundle, model_metadata = adapter.saved_chromosome(
                    manifest, adapter.config_for(candidate), chrom, directory)
                model_hashes[str(chrom)] = model_metadata["model_sha256"]
                model, keys, layout = bundle["model"], bundle["keys"], bundle["layout"]
                if any(not np.array_equal(m.classes_, np.arange(len(adapter.SOURCES)))
                       for m in model.base.models):
                    raise ValueError("Unexpected base-classifier ancestry axis")
                positions = output / f"chr{chrom}.positions.tsv"
                positions.write_text("".join(f"{key[0]}\t{key[1]}\n" for key in keys))
                arrays = [adapter.load_haplotypes(batch["vcf"], [row["sample_id"] for row in rows],
                                                  chrom, keys) for _, batch, rows in batches]
                x = np.concatenate(arrays)
                del arrays
                qc = genotype_qc(x)
                for name, values in qc.items():
                    qc_total[name] += values
                if np.any(qc["missing_sites"]):
                    raise ValueError("Missing genotypes require a separately designed diagnostic")
                real_vcf = manifest["evaluation"]["real"]["batches"][0]["vcf"]
                info = query(real_vcf, "%CHROM\t%POS\t%REF\t%ALT\t%INFO/IMP\t%INFO/DR2"
                             "[\t%GT\t%GP]\n", positions)
                info_keys = [(row[0], int(row[1]), row[2], row[3]) for row in info]
                if info_keys != keys:
                    raise ValueError("Real genotype information lost marker alignment")
                observed = np.asarray([row[4] != "1" for row in info])
                dr2 = np.asarray([float(row[5]) if row[5] != "." else np.nan for row in info])
                gp = np.asarray([[float(value) for value in row[7].split(",")] for row in info])
                if gp.shape != (len(keys), 3):
                    raise ValueError("Unexpected genotype posterior count")
                gp_sums = genotype_posterior_qc(gp)
                canonical_gt = [row[6] for row in info]
                source = Path("/data/flare2_test/source_model/real_user_validation/lonni/"
                              "lonni.autosomes.vcf.gz")
                source_rows_gt = query(source, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", positions)
                source_map = {tuple(row[:4]): row[4] for row in source_rows_gt}
                mismatches = sum(source_map.get(tuple(row[:4])) != row[6] for row in info)
                upstream_path = upstream / f"qc_output/IMPID325.chr{chrom}.split_QC.vcf.gz"
                if not upstream_path.is_file():
                    raise ValueError(f"Upstream chromosome file is absent: {upstream_path}")
                audit["upstream_sha256"][str(upstream_path)] = adapter.digest(upstream_path)
                upstream_rows = query(upstream_path, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", positions)
                upstream_map = {tuple(row[:4]): row[4] for row in upstream_rows}
                upstream_mismatches = sum(upstream_map.get(tuple(row[:4])) != row[6] for row in info)
                if mismatches or upstream_mismatches:
                    raise ValueError("Raw-to-canonical genotype or phased-allele concordance failed")
                dosage = x.reshape(-1, 2, len(keys)).sum(axis=1)
                gt_info = np.asarray([adapter.phased_gt(gt) for gt in canonical_gt]).T
                if not np.array_equal(gt_info, x[2 * r:2 * r + 2]):
                    raise ValueError("Real sample axis or GT interpretation changed")
                heterozygous = dosage[r] == 1
                for label, mask in (("observed", observed), ("imputed", ~observed),
                                    ("homozygous", ~heterozygous), ("heterozygous", heterozygous)):
                    raw_summary.append(dict(
                        chrom=chrom, stratum=label, n_sites=int(mask.sum()),
                        dr2_available=int(np.isfinite(dr2[mask]).sum()),
                        dr2_sum=float(np.nansum(dr2[mask])),
                        dr2_zero=int((dr2[mask] == 0).sum()),
                        gp_max_sum=float(gp[mask].max(axis=1).sum()) if mask.any() else 0.,
                        gp_rounded_sum_sites=int((np.abs(gp_sums[mask] - 1) > 2e-4).sum()),
                        source_phased_gt_mismatches=mismatches,
                        upstream_phased_gt_mismatches=upstream_mismatches))
                base = model.base.predict_proba(x)
                smooth = adapter.ordered_probabilities(model, base)
                if not np.allclose(crf_in_blocks(model, base, layout["W"]), smooth,
                                   rtol=0, atol=1e-10):
                    raise ValueError("Whole-chromosome block diagnostic differs from saved inference")
                for size in (20, 100):
                    values = crf_in_blocks(model, base, size)
                    paired, _ = adapter.aggregate_probabilities(values, layout)
                    stage = f"crf_blocks_{size}"
                    if stage not in totals:
                        totals[stage] = np.zeros_like(paired)
                    totals[stage] += paired
                for stage, values, indices in (("base", base, np.arange(len(all_rows))),
                                               ("smooth", smooth, np.arange(len(all_rows)))):
                    paired, _ = adapter.aggregate_probabilities(values, layout)
                    key = stage
                    if key not in totals:
                        totals[key] = np.zeros_like(paired)
                    totals[key] += paired
                    for index in real_index + matched:
                        chromosome_rows.append(dict(
                            dataset=all_rows[index]["dataset"], sample_id=all_rows[index]["sample_id"],
                            chrom=chrom, stage=stage, haplotype="diploid", n_markers=len(keys),
                            **dict(zip(adapter.SOURCES, paired[index] / len(keys)))))
                    for hap in (0, 1):
                        fractions = np.einsum("wa,w->a", values[2 * r + hap],
                                              layout["marker_weights"])
                        chromosome_rows.append(dict(
                            dataset="real", sample_id=all_rows[r]["sample_id"], chrom=chrom,
                            stage=stage, haplotype=hap + 1, n_markers=len(keys),
                            **dict(zip(adapter.SOURCES, fractions / len(keys)))))
                selected = x.reshape(-1, 2, len(keys))[phase_indices].reshape(-1, len(keys))
                randomized = randomize_phase(selected, adapter.SEED + chrom)
                changed_base = model.base.predict_proba(randomized)
                changed_smooth = adapter.ordered_probabilities(model, changed_base)
                for stage, values in (("random_phase_base", changed_base),
                                      ("random_phase_smooth", changed_smooth)):
                    paired, _ = adapter.aggregate_probabilities(values, layout)
                    if stage not in totals:
                        totals[stage] = np.zeros_like(paired)
                    totals[stage] += paired
                swapped = x[2 * r:2 * r + 2][::-1]
                swapped_p = adapter.ordered_probabilities(model, model.base.predict_proba(swapped))
                if not np.allclose(swapped_p, smooth[2 * r:2 * r + 2][::-1], rtol=0, atol=1e-10):
                    raise ValueError("Swapping entire haplotypes changed the sample's posteriors")
                fkeys, frequencies = reference_frequencies(manifest, chrom, groups, output)
                if fkeys != keys:
                    raise ValueError("Reference frequency marker order changed")
                thinned = thin_sites(keys, np.ones(len(keys), dtype=bool))
                frequency_parts.append(frequencies)
                dosage_parts.append(dosage)
                observed_parts.append(observed)
                thinned_parts.append(thinned)
                observed_thinned_parts.append(thin_sites(keys, observed))
                marker_total += len(keys)
                audit["chromosomes"].append(dict(chromosome=chrom, n_markers=len(keys),
                                                  observed_sites=int(observed.sum())))
                adapter.write_json(output / "audit.json", audit)
                print(f"chr{chrom}: genotype provenance, posteriors and frequency controls verified",
                      flush=True)
                del bundle, x, base, smooth, changed_base, changed_smooth
        predictions = []
        for stage, values in totals.items():
            rows = all_rows if not stage.startswith("random_phase") else [
                all_rows[index] for index in phase_indices]
            for row, fractions in zip(rows, values / marker_total):
                predictions.append(dict(dataset=row["dataset"], sample_id=row["sample_id"],
                                        stage=stage, **dict(zip(adapter.SOURCES, fractions))))
        frame = pd.DataFrame(predictions)
        saved = pd.read_csv(directory / "predictions.tsv", sep="\t")
        replicated = frame[frame.stage.eq("smooth")].merge(
            saved, on=["dataset", "sample_id"], validate="one_to_one")
        error = float(np.max(np.abs(replicated[list(adapter.SOURCES)].to_numpy()
                                    - replicated[["pred_" + s for s in adapter.SOURCES]].to_numpy())))
        if len(replicated) != len(saved) or error > 1e-10:
            raise ValueError(f"Audit did not reproduce the saved candidate: {error}")
        frame.to_csv(output / "stage_predictions.tsv", sep="\t", index=False)
        pd.DataFrame(chromosome_rows).to_csv(output / "chromosome_predictions.tsv", sep="\t", index=False)
        qc = pd.DataFrame(qc_total)
        qc.insert(0, "sample_id", [row["sample_id"] for row in all_rows])
        qc.insert(0, "dataset", [row["dataset"] for row in all_rows])
        qc.to_csv(output / "genotype_qc.tsv", sep="\t", index=False)
        pd.DataFrame(raw_summary).to_csv(output / "real_input_qc.tsv", sep="\t", index=False)
        frequencies = np.concatenate(frequency_parts)
        dosages = np.concatenate(dosage_parts, axis=1)
        observed = np.concatenate(observed_parts)
        thinned = np.concatenate(thinned_parts)
        observed_thinned = np.concatenate(observed_thinned_parts)
        frequency_predictions = []
        for label, mask in (("all_markers", np.ones(marker_total, dtype=bool)),
                            ("observed_sites", observed), ("imputed_sites", ~observed),
                            ("thinned_100kb", thinned),
                            ("observed_thinned_100kb", observed_thinned)):
            fits, informative = frequency_fits(frequencies[mask], dosages[:, mask])
            for row, fractions in zip(all_rows, fits):
                frequency_predictions.append(dict(
                    dataset=row["dataset"], sample_id=row["sample_id"], panel=label,
                    n_markers=int(mask.sum()), n_informative=informative,
                    **dict(zip(adapter.SOURCES, fractions))))
        pd.DataFrame(frequency_predictions).to_csv(
            output / "frequency_diagnostics.tsv", sep="\t", index=False)
        metrics = []
        truths = []
        for row in all_rows:
            if row["dataset"] == "real":
                truths.append(np.zeros(len(adapter.SOURCES)))
            else:
                values = {s: float(row["realized_" + s]) for s in
                          ("EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA")}
                truths.append([values[s] if s != "EAS_SEA" else values["EAS"] + values["SEA"]
                               for s in adapter.SOURCES])
        truth = np.asarray(truths)
        for stage, values in totals.items():
            indices = phase_indices if stage.startswith("random_phase") else np.arange(len(all_rows))
            for dataset in ("pure", "pedigree", "real"):
                selected = np.asarray([all_rows[index]["dataset"] == dataset for index in indices])
                if not selected.any():
                    continue
                p = values[selected] / marker_total
                result = dict(stage=stage, dataset=dataset, n=int(selected.sum()),
                              mean_prediction=dict(zip(adapter.SOURCES, p.mean(axis=0))))
                if dataset != "real":
                    tv = np.abs(p - truth[indices[selected]]).sum(axis=1) / 2
                    result.update(mean_tv=float(tv.mean()), p95_tv=float(np.quantile(tv, .95)))
                metrics.append(result)
            if not stage.startswith("random_phase"):
                tv = np.abs(values[matched] / marker_total - truth[matched]).sum(axis=1) / 2
                metrics.append(dict(
                    stage=stage, dataset="matched_50_50_eur_asian", n=len(matched),
                    mean_tv=float(tv.mean()), p95_tv=float(np.quantile(tv, .95)),
                    mean_prediction=dict(zip(adapter.SOURCES,
                                             (values[matched] / marker_total).mean(axis=0)))))
        frequency_frame = pd.DataFrame(frequency_predictions).merge(
            metadata[["dataset", "sample_id"]].assign(_index=np.arange(len(all_rows))),
            on=["dataset", "sample_id"], validate="many_to_one")
        frequency_metrics = []
        for (panel, dataset), rows in frequency_frame.groupby(["panel", "dataset"]):
            if dataset == "real":
                continue
            p = rows[list(adapter.SOURCES)].to_numpy()
            tv = np.abs(p - truth[rows._index]).sum(axis=1) / 2
            frequency_metrics.append(dict(panel=panel, dataset=dataset, n=len(rows),
                                          mean_tv=float(tv.mean()),
                                          p95_tv=float(np.quantile(tv, .95))))
        audit.update(status="complete", n_markers=marker_total,
                     exact_saved_predictions_reproduced=True, maximum_saved_prediction_error=error,
                     matched_50_50_eur_asian_controls=[all_rows[i]["sample_id"] for i in matched],
                     model_sha256=model_hashes, stage_metrics=metrics,
                     frequency_diagnostic_metrics=frequency_metrics,
                     crf_block_control="20/100 consecutive SNP-count windows, unchanged fitted weights; not retrained candidates")
        if adapter.input_hashes(manifest) != hashes or adapter.digest(__file__) != audit["script_sha256"]:
            raise ValueError("Audit inputs or implementation changed during execution")
        for path, expected in audit["upstream_sha256"].items():
            if adapter.digest(path) != expected:
                raise ValueError("An upstream genotype file changed during the audit")
        for chrom, expected in model_hashes.items():
            if adapter.digest(directory / f"models/chr{chrom}.pkl") != expected:
                raise ValueError("A frozen model changed during the audit")
        audit["output_sha256"] = {p.name: adapter.digest(p) for p in output.glob("*.tsv")}
        adapter.write_json(output / "audit.json", audit)
        return audit
    except Exception as error:
        audit.update(status="failed", error=str(error))
        adapter.write_json(output / "audit.json", audit)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--candidate", default="gnomix_fast_050",
                        choices=("gnomix_fast_020", "gnomix_fast_050"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.screen, args.outdir, args.upstream, args.candidate)


if __name__ == "__main__":
    main()
