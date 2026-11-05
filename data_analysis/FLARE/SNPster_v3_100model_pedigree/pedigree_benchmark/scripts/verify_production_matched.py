#!/usr/bin/env python3
"""Independent coverage, Mendelian transmission, local aggregation and TV audit."""
import argparse
import gzip
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

from audit_production_ancestry import SOURCES, digest, write_json
from production_matched_preprocessing import PROJECT, ROOT, read, read_call_matrix, reference
from verify_sgdp_full_manifest import require_hashes


def vector(value):
    value = np.asarray(value, dtype=float)
    if (value.shape != (8,) or not np.isfinite(value).all() or (value < 0).any()
            or abs(value.sum() - 1) > 1e-10):
        raise ValueError("Invalid eight-source fractions in independent scorer")
    return value


def merged(value):
    value = vector(value)
    return value[[0, 1, 2, 4, 5, 6, 7]] + np.asarray([0, value[3], 0, 0, 0, 0, 0])


def error_pct(prediction, truth, collapse=False):
    a, b = vector(prediction), vector(truth)
    if collapse:
        a, b = merged(a), merged(b)
    return float(50 * np.abs(a - b).sum())


def gate(cohort_path, chromosomes, masks, reports, output):
    cohort = read(cohort_path)
    expected = set(itertools.product(chromosomes, masks))
    seen = set()
    for path in reports:
        record = read(path)
        key = (record["chromosome"], record["mask"])
        if (key in seen or key not in expected or record["status"] != "production_preprocessing_genotypes_and_coverage_verified"
                or record["people"] != len(cohort["people"]) or record["production_changed"]):
            raise ValueError("Incomplete, duplicated or invalid production preprocessing checks")
        require_hashes(record["input_sha256"])
        seen.add(key)
    if seen != expected:
        raise ValueError("All declared chromosomes and masks must pass before any model testing")
    write_json(output, dict(status="all_preprocessing_checks_passed", chromosomes=chromosomes, masks=masks,
                            reports_verified=len(reports), cohort_sha256=digest(cohort_path)))
    print(f"Production genotype gate passed for all {len(reports)} chromosome/mask batches.", flush=True)


def local_flare(path, panel):
    keys = [tuple(key) for key in panel["keys"]]
    from screen_model_families import read_model
    names = read_model(str(path).replace(".anc.vcf.gz", ".model"))["ancestries"]
    if set(names) != set(SOURCES) or len(names) != 8:
        raise ValueError("FLARE source names do not match eight-source contract")
    order = [names.index(source) for source in SOURCES]
    totals = {identity: np.zeros(8) for identity in panel["evaluation"]}
    previous, edge, count = None, float(panel["span"][0]), 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                identities = line.rstrip().split("\t")[9:]
                if identities != panel["evaluation"]:
                    raise ValueError("Independent local ancestry sample axis mismatch")
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if count >= len(keys) or key != keys[count]:
                raise ValueError("Independent local ancestry exact marker axis mismatch")
            columns = fields[8].split(":")
            values = []
            for call in fields[9:]:
                parts = call.split(":")
                copies = [np.asarray([float(x) for x in parts[columns.index(name)].split(",")])
                          for name in ("ANP1", "ANP2")]
                if any(len(copy) != 8 or not np.isfinite(copy).all() or (copy < 0).any()
                       or abs(copy.sum() - 1) > .005 for copy in copies):
                    raise ValueError("Independent FLARE probability validation failed")
                values.append(sum(copy[order] / copy.sum() for copy in copies) / 2)
            if previous is not None:
                boundary = (previous[0] + key[1]) / 2
                for identity, posterior in zip(identities, previous[1]):
                    totals[identity] += posterior * (boundary - edge)
                edge = boundary
            previous = (key[1], values)
            count += 1
    if count != len(keys) or previous is None:
        raise ValueError("Incomplete independent local ancestry evidence")
    for identity, posterior in zip(identities, previous[1]):
        totals[identity] += posterior * (panel["span"][1] - edge)
        totals[identity] /= panel["span"][1] - panel["span"][0]
    return totals


def local_gnomix(path, panel, window_cm):
    data = np.load(path)
    probabilities, supplied_weights = data["probabilities"], data["weights"]
    positions = np.asarray([key[1] for key in panel["keys"]], dtype=float)
    chromosome = panel["chromosome"]
    map_path = ROOT / "plink.GRCh38.map" / f"plink.chr{chromosome}.GRCh38.map"
    with map_path.open() as handle:
        cm = [float(fields[2]) for line in handle
              if (fields := line.split()) and fields[0].removeprefix("chr") == str(chromosome)]
    width = round(window_cm * len(positions) / max(cm))
    windows = len(positions) // width
    cells = np.diff(np.r_[panel["span"][0], (positions[:-1] + positions[1:]) / 2, panel["span"][1]])
    independent_weights = np.asarray([
        cells[index * width:(index + 1) * width if index < windows - 1 else len(cells)].sum()
        for index in range(windows)])
    np.testing.assert_array_equal(supplied_weights, independent_weights)
    if (probabilities.shape != (2 * len(panel["evaluation"]), windows, 8)
            or not np.isfinite(probabilities).all() or (probabilities < 0).any()
            or not np.allclose(probabilities.sum(axis=2), 1, rtol=0, atol=1e-10)):
        raise ValueError("Independent Gnomix posterior shape/probability check failed")
    result = {}
    for index, identity in enumerate(panel["evaluation"]):
        result[identity] = sum(
            (probabilities[2 * index, window] + probabilities[2 * index + 1, window]) *
            independent_weights[window] / 2 for window in range(windows))
        result[identity] /= panel["span"][1] - panel["span"][0]
    return result


def validate_qp(frequencies, dosage, prediction):
    keep = dosage >= 0
    if keep.sum() < 100:
        raise ValueError("Insufficient independently called phase-free markers")
    a, b = 2 * frequencies[keep], dosage[keep]
    x = vector(prediction)
    gradient = 2 * a.T @ (a @ x - b) / len(b)
    active = x > 1e-6
    dual = float(gradient[active].mean())
    if (np.max(np.abs(gradient[active] - dual)) > 1e-5
            or np.any(gradient[~active] < dual - 1e-5)):
        raise ValueError("Phase-independent prediction fails independent convex-optimality conditions")


def inheritance(cohort, run_dir, chromosomes, output):
    if cohort.get("truth_design") == "synthetic_reference_haplotype_development":
        return synthetic_inheritance(cohort, run_dir, chromosomes, output)
    cells = 0
    for chromosome in chromosomes:
        folder = run_dir / f"founders/chr{chromosome}/founders"
        metadata = read(folder / "founders.json")
        keys = [tuple(key) for key in metadata["keys"]]
        original = np.load(folder / "observed.npz")["calls"]
        phased_path = run_dir / f"founder_phase/chr{chromosome}/phased/imputed_chr{chromosome}.risk.vcf.gz"
        phased = read_call_matrix(phased_path, keys, metadata["founder_ids"], output, phased=True)
        np.testing.assert_array_equal(original.sum(axis=1), phased.sum(axis=1))
        truth_path = run_dir / f"truth/chr{chromosome}/truth_chr{chromosome}.vcf.gz"
        truth = read_call_matrix(truth_path, keys, [row["sample_id"] for row in cohort["people"]],
                                 output, phased=True)
        offsets = {row["sample_id"]: index for index, row in enumerate(metadata["people"])}
        for row, value in zip(cohort["people"], truth):
            donors = row["donors"]
            if row["kind"] == "pure":
                np.testing.assert_array_equal(value, phased[offsets[donors[0]]])
            else:
                if len(set(donors)) != 2:
                    raise ValueError("An F1 control reuses a donor as both parents")
                for copy_index, donor in enumerate(donors):
                    if not any(np.array_equal(value[copy_index], hap) for hap in phased[offsets[donor]]):
                        raise ValueError("F1 allele transmission does not match a complete parental homolog")
            cells += len(keys)
    return cells


def synthetic_inheritance(cohort, run_dir, chromosomes, output):
    from realized_ancestry import load_bp

    cells = 0
    sources = {row["sample_id"]: row["source"]
               for row in cohort["training_founders"] + cohort["evaluation_founders"]}
    for chromosome in chromosomes:
        folder = run_dir / f"founders/chr{chromosome}/founders"
        metadata = read(folder / "founders.json")
        keys = [tuple(key) for key in metadata["keys"]]
        original = np.load(folder / "template.npz")["calls"]
        recovered = read_call_matrix(reference(chromosome), keys, metadata["founder_ids"], output, phased=True)
        np.testing.assert_array_equal(original, recovered)
        truth = read_call_matrix(run_dir / f"truth/chr{chromosome}/truth_chr{chromosome}.vcf.gz",
                                 keys, [row["sample_id"] for row in cohort["people"]], output, phased=True)
        offsets = {identity: index for index, identity in enumerate(metadata["founder_ids"])}
        positions = np.asarray([key[1] for key in keys])
        span = cohort["chromosome_spans"][str(chromosome)]
        decoded = {}
        for index, person in enumerate(cohort["people"]):
            if person["kind"] == "pure":
                parts = [dict(haplotype=hap, donor=person["donors"][0], donor_haplotype=hap,
                              start=span[0], end=span[1], source=sources[person["donors"][0]])
                         for hap in (0, 1)]
            else:
                plan_key = person["plan_key"]
                if plan_key not in decoded:
                    role, batch = plan_key.split(":")
                    root = PROJECT / ("results_v3_pedigree_training_screen_01/cohort/pedigree"
                                      if role == "training" else "results_v3_pedigree_pilot_10/pedigree")
                    plan = read(root / f"plans/plans/{batch}/pedigrees.json")
                    bp = load_bp(root / f"chromosomes/pdb_{batch}_chr{chromosome}.bp", chromosome)
                    mapping = {}
                    for founder, info in plan["founders"].items():
                        for hap_index, hap in enumerate(("h0", "h1")):
                            ids = {part[2] for part in bp[founder][hap]}
                            if len(ids) != 1 or next(iter(ids)) in mapping:
                                raise ValueError("Independent founder homolog mapping is ambiguous")
                            mapping[next(iter(ids))] = (info["vcf_sample"], hap_index)
                    decoded[plan_key] = bp, mapping
                bp, mapping = decoded[plan_key]
                parts = []
                for hap_index, hap in enumerate(("h0", "h1")):
                    previous = span[0]
                    for start, inclusive_end, donor_hap in bp[person["sample_id"]][hap]:
                        first, end = max(start, span[0]), min(inclusive_end + 1, span[1])
                        if first >= end:
                            continue
                        if first != previous:
                            raise ValueError("Independent physical ancestry has a gap or overlap")
                        donor, copy = mapping[donor_hap]
                        parts.append(dict(haplotype=hap_index, donor=donor, donor_haplotype=copy,
                                          start=first, end=end, source=sources[donor]))
                        previous = end
                    if previous != span[1]:
                        raise ValueError("Independent physical ancestry is incomplete")
                if parts != cohort["local_ancestry"][str(chromosome)][person["sample_id"]]:
                    raise ValueError("Frozen local ancestry differs from original Ped-sim breakpoints")
            covered = np.zeros((2, len(keys)), dtype=int)
            amounts = np.zeros(8, dtype=np.int64)
            for part in parts:
                use = (positions >= part["start"]) & (positions < part["end"])
                np.testing.assert_array_equal(truth[index, part["haplotype"], use],
                    original[offsets[part["donor"]], part["donor_haplotype"], use])
                covered[part["haplotype"], use] += 1
                amounts[SOURCES.index(part["source"])] += part["end"] - part["start"]
            if not np.all(covered == 1) or amounts.sum() != 2 * (span[1] - span[0]):
                raise ValueError("Independent synthetic transmission/physical truth coverage failed")
            np.testing.assert_allclose(amounts / amounts.sum(),
                cohort["truth_by_chromosome"][str(chromosome)][person["sample_id"]], rtol=0, atol=1e-12)
            cells += 2 * len(keys)
    return cells


def aggregate_truth(cohort, person, chromosomes, lengths):
    if cohort.get("truth_design") == "synthetic_reference_haplotype_development":
        values = np.asarray([cohort["truth_by_chromosome"][str(int(chromosome))][person["sample_id"]]
                             for chromosome in chromosomes])
        return vector((values * np.asarray(lengths)[:, None]).sum(axis=0) / np.sum(lengths))
    sources = {row["sample_id"]: row["source"] for row in cohort["evaluation_founders"]}
    value = np.zeros(8)
    for donor in person["donors"]:
        value[SOURCES.index(sources[donor])] += 1 / len(person["donors"])
    return vector(value)


def score(cohort_path, run_dir, chromosomes, masks, folders, output):
    output.mkdir()
    cohort = read(cohort_path)
    require_hashes(cohort["input_sha256"])
    people = [row for row in cohort["people"] if row["role"] == "evaluation"]
    person_ids = [row["sample_id"] for row in people]
    candidates = cohort["protocol"]["candidates"]
    expected = set(itertools.product(candidates, masks, chromosomes, person_ids))
    seen, frames, bound = set(), [], {}
    for folder in folders:
        complete = read(folder / "complete.json")
        chromosome, candidate = complete["chromosome"], complete["candidate"]
        panel_path = run_dir / f"panels/chr{chromosome}/panel/panel.json"
        panel = read(panel_path)
        if (digest(folder / "predictions.tsv") != complete["predictions_sha256"]
                or digest(cohort_path) != complete["cohort_sha256"]
                or digest(panel_path) != complete["panel_sha256"]):
            raise ValueError("Published model/panel/cohort evidence differs from fitted inputs")
        frame = pd.read_csv(folder / "predictions.tsv", sep="\t", float_precision="round_trip")
        if len(frame) != len(masks) * len(people):
            raise ValueError("Incomplete per-model prediction rows")
        policy = candidates[candidate]
        for mask in masks:
            if policy["family"] == "flare":
                recalculated = local_flare(folder / f"{mask}.anc.vcf.gz", panel)
            elif policy["family"] == "gnomix":
                recalculated = local_gnomix(folder / f"{mask}.local_posteriors.npz", panel,
                                            policy["window_size_cM"])
            else:
                if policy["input"] == "production_phased_imputed":
                    data = np.load(run_dir / f"panels/chr{chromosome}/panel/{mask}.evaluation.npz")
                    dosage = data["haplotypes"].reshape(-1, 2, panel["marker_count"]).sum(axis=1)
                    frequencies = np.load(run_dir / f"panels/chr{chromosome}/panel/reference.npz")["frequencies"]
                else:
                    data = np.load(run_dir / f"panels/chr{chromosome}/panel/{mask}.observed.npz")
                    dosage, frequencies = data["dosage"], data["frequencies"]
            for index, identity in enumerate(person_ids):
                matches = frame[(frame["mask"] == mask) & (frame["sample_id"] == identity)]
                if len(matches) != 1:
                    raise ValueError("Duplicate or missing model/person/mask prediction")
                row = matches.iloc[0]
                key = (row["candidate"], row["mask"], int(row["chromosome"]), row["sample_id"])
                if key not in expected or key in seen or key[0] != candidate:
                    raise ValueError("Prediction Cartesian coverage is not exact")
                if int(row["span_bp"]) != panel["span"][1] - panel["span"][0]:
                    raise ValueError("Incorrect physical chromosome span in prediction")
                seen.add(key)
                values = vector(row[list(SOURCES)].to_numpy(dtype=float))
                if policy["family"] == "phasefree":
                    validate_qp(frequencies, dosage[index], values)
                else:
                    np.testing.assert_allclose(values, vector(recalculated[identity]), rtol=0, atol=1e-10)
        for path in folder.rglob("*"):
            if path.is_file():
                destination = run_dir / f"models/{candidate}/chr{chromosome}/model" / path.relative_to(folder)
                bound[str(destination)] = digest(path)
        frames.append(frame)
    if seen != expected:
        raise ValueError(f"Partial benchmark cannot be scored: {len(seen)}/{len(expected)} predictions")
    genotype_cells = inheritance(cohort, run_dir, chromosomes, output)
    frame = pd.concat(frames, ignore_index=True)
    frame.to_csv(output / "chromosome_predictions.tsv", sep="\t", index=False)
    metadata = {row["sample_id"]: row for row in people}
    rows = []
    for (candidate, mask, identity), subset in frame.groupby(["candidate", "mask", "sample_id"]):
        if sorted(subset["chromosome"].tolist()) != sorted(chromosomes):
            raise ValueError("Incomplete chromosome coverage during physical aggregation")
        lengths = subset["span_bp"].to_numpy()
        prediction = vector((subset[list(SOURCES)].to_numpy() * lengths[:, None]).sum(axis=0) / lengths.sum())
        truth = aggregate_truth(cohort, metadata[identity], subset["chromosome"], lengths)
        rows.append(dict(candidate=candidate, mask=mask, sample_id=identity,
                         kind=metadata[identity]["kind"], family_id=metadata[identity]["family_id"],
                         eight_error_pct=error_pct(prediction, truth),
                         seven_error_pct=error_pct(prediction, truth, collapse=True),
                         **{source: float(prediction[index]) for index, source in enumerate(SOURCES)}))
    per_mask = pd.DataFrame(rows)
    per_mask.to_csv(output / "person_mask_predictions.tsv", sep="\t", index=False)
    per_person = per_mask.groupby(["candidate", "sample_id", "kind", "family_id"], as_index=False)[
        ["eight_error_pct", "seven_error_pct"]].mean()
    per_person.to_csv(output / "person_errors.tsv", sep="\t", index=False)
    summary = per_person.groupby("candidate", as_index=False)[["eight_error_pct", "seven_error_pct"]].mean()
    summary.to_csv(output / "summary.tsv", sep="\t", index=False)
    per_person.groupby(["candidate", "kind"], as_index=False)[["eight_error_pct", "seven_error_pct"]].mean().to_csv(
        output / "diagnostic_pure_mixed.tsv", sep="\t", index=False)
    whole_genome = sorted(chromosomes) == list(range(1, 23)) and masks == cohort["protocol"]["masks"]
    write_json(output / "independent_verification.json", dict(
        status="complete_verified_development_comparison" if whole_genome else "smoke_passed_not_a_whole_genome_leaderboard",
        model_configurations=len(candidates), chromosomes=chromosomes, masks=masks,
        n_evaluation_people=len(people), n_evaluation_founders=len(cohort["evaluation_founders"]),
        independent_founder_family_components=None if cohort.get("truth_design") else 1,
        truth_design=cohort.get("truth_design", "observed_founder_controls"),
        requested_counts=cohort.get("requested_counts"), retained_counts=cohort.get("retained_counts"),
        predictions_verified=len(seen),
        mendelian_genotype_cells_verified=genotype_cells, independently_verified_model_sha256=bound,
        cohort_sha256=digest(cohort_path), production_changed=False,
        independent_final_test=False, production_promotion=False,
        headline="equal_person_mean_of_complete_mask_mean_tv_error_pct",
        reported_resolutions=["eight_category", "seven_category"]))
    print(summary.to_string(index=False), flush=True)
    print("Independent ancestry evidence, Mendelian transmission and both TV resolutions verified.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("gate", "score"))
    parser.add_argument("--cohort", required=True, type=Path)
    parser.add_argument("--chromosomes", type=int, nargs="+", required=True)
    parser.add_argument("--masks", nargs="+", required=True)
    parser.add_argument("--reports", type=Path, nargs="+")
    parser.add_argument("--models", type=Path, nargs="+")
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "gate":
        gate(args.cohort, args.chromosomes, args.masks, args.reports, args.output)
    else:
        score(args.cohort, args.run_dir, args.chromosomes, args.masks, args.models, args.output)
