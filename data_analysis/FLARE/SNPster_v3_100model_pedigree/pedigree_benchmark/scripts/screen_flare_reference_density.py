#!/usr/bin/env python3
"""Training-only paired FLARE reference representation and marker-density screen."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import benchmark_dense_observed as dense
import benchmark_observed_ancestry as common
import diagnose_flare_native as diagnostic
import model_family_gnomix as native
import pilot_matched_imputation as shared
from screen_model_families import physical_predictions, read_model, write_model
from verify_dense_observed import check_case, local_physical


PROJECT = Path(__file__).resolve().parents[2]
PROTOCOL = PROJECT / "config/flare_reference_density.json"
BASE = Path("/benchmark")
SOURCES = tuple(native.SOURCES)


def rank(identity, seed, tag):
    return hashlib.sha256(f"{seed}:{tag}:{identity}".encode()).digest()


def reference_sets(rows, cap, seed):
    if cap < 1 or len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Invalid reference cap or duplicate reference identity")
    fine, fine_sources = defaultdict(list), {}
    for row in rows:
        source, population = row["source_group"], row["fine"]
        if source not in SOURCES or not population:
            raise ValueError("Invalid source/fine reference label")
        if population in fine_sources and fine_sources[population] != source:
            raise ValueError("Fine population belongs to multiple broad sources")
        fine_sources[population] = source
        fine[population].append(row)
    capped = [row for population in sorted(fine) for row in sorted(
        fine[population], key=lambda r: rank(r["sample_id"], seed, "fine-cap"))[:cap]]
    counts = Counter(r["source_group"] for r in capped)
    if set(counts) != set(SOURCES):
        raise ValueError("Reference representation loses a source")
    control = [row for source in SOURCES for row in sorted(
        (r for r in rows if r["source_group"] == source),
        key=lambda r: rank(r["sample_id"], seed, "size-control"))[:counts[source]]]
    return {name: sorted(r["sample_id"] for r in values)
            for name, values in (("full", rows), ("finecap16", capped), ("size_control", control))}


def thin_indices(positions, count):
    positions = np.asarray(positions)
    if (count < 2 or count > len(positions) or np.any(np.diff(positions) <= 0)):
        raise ValueError("Thinning requires sorted distinct positions and a valid count")
    targets = np.linspace(positions[0], positions[-1], count)
    result, previous = [], -1
    for i, target in enumerate(targets):
        index = int(np.clip(np.searchsorted(positions, target), previous + 1, len(positions) - count + i))
        result.append(index)
        previous = index
    return result


def chromosome_truth(path, identities, chrom, span):
    pieces = defaultdict(lambda: defaultdict(list))
    with Path(path).open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["sample_id"] in identities and row["chrom"].removeprefix("chr") == str(chrom):
                source = row["source_group"]
                source = "EAS_SEA" if source in ("EAS", "SEA") else source
                if source not in SOURCES or row["haplotype"] not in ("1", "2"):
                    raise ValueError("Invalid pedigree segment ancestry/haplotype")
                pieces[row["sample_id"]][row["haplotype"]].append(
                    (int(row["start_bp"]), int(row["end_bp"]), source))
    result = {}
    for identity in identities:
        if set(pieces[identity]) != {"1", "2"}:
            raise ValueError(f"Missing chromosome-specific pedigree truth: {identity}")
        mass = np.zeros(len(SOURCES))
        for segments in pieces[identity].values():
            cursor = span[0]
            for start, end, source in sorted(segments):
                if start != cursor or not start < end <= span[1]:
                    raise ValueError(f"Pedigree truth has gaps/overlap/wrong span: {identity}")
                mass[SOURCES.index(source)] += end - start
                cursor = end
            if cursor != span[1]:
                raise ValueError("Pedigree truth does not cover the full chromosome span")
        result[identity] = (mass / (2 * (span[1] - span[0]))).tolist()
    return result


def markers(task):
    return [(r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            for r in native.records(task / "markers.tsv")]


def candidates(protocol):
    return [(f"{reference}_{'em' if em else 'fixed'}", reference, em)
            for reference, em in itertools.product(protocol["reference_representations"], protocol["em_modes"])]


def command(cohort, target, prefix, em, model=None):
    result = diagnostic.command(cohort, target, prefix, em=em, model=model)
    if any(str(arg).startswith("gt-ancestries=") for arg in result):
        raise ValueError("Target truth must not be passed to inference")
    return result


def prepare(out):
    if (out / "experiment.json").exists():
        raise ValueError("Fresh output required; prepared experiments are immutable")
    protocol = common.read(PROTOCOL)
    parent = common.read(BASE / "experiment.json")
    shared.verify_hashes(parent["implementation_sha256"])
    selection = common.read(BASE / "selection.json")
    pure = native.records(BASE / "inputs/pure_train_labels.tsv")
    reference = native.records(BASE / "inputs/reference_labels.tsv")
    diagnostic.source_audit(
        common.read(BASE / "tasks/lonni_chr22/cohort.json")["reference"]["vcf"],
        reference, pure, pure)
    jar_audit = diagnostic.jar_check(common.read(BASE / "tasks/lonni_chr22/cohort.json")["flare_jar"])
    representations = reference_sets(reference, protocol["reference_fine_cap"], protocol["selection_seed"])
    inputs = [PROTOCOL, BASE / "experiment.json", BASE / "selection.json",
              BASE / "inputs/pure_train_labels.tsv", BASE / "inputs/reference_labels.tsv",
              BASE / "inputs/reference.panel", Path("/data/flare2_test/source_model/ancestry_role_audit.tsv")]
    people = [dict(r, track="pure", family_id=r["sample_id"], scenario="pure") for r in pure]
    donors = common.read(PROJECT / "results_v3_phasefree_rich_01/donors.json")
    pure_ids, families = {r["sample_id"] for r in pure}, set()
    for batch in selection["pedigree_train"]:
        inputs += [Path(batch["samples_tsv"]), Path(batch["segments_tsv"])]
        for row in native.records(batch["samples_tsv"]):
            if not set(donors[row["family_id"]]) <= pure_ids:
                raise ValueError("Incomplete or non-training family founders")
            families.add(row["family_id"])
            people.append(dict(row, track="pedigree", batch=batch["batch"], fine="pedigree",
                               source_group="mixture", scenario=row["class"]))
    forbidden = set(selection["evaluation_donors"]) | set(selection["reference_ids"])
    if (len(pure) != 240 or len(people) != 316 or len(families) != 25
            or len({r["sample_id"] for r in people}) != len(people)
            or forbidden & {r["sample_id"] for r in people}):
        raise ValueError("Training-only identity/family contract changed")
    truth = {}
    for chrom in protocol["chromosomes"]:
        for mask in ("lonni", "frederik"):
            task = BASE / f"tasks/{mask}_chr{chrom}"
            checkpoint = common.read(task / "inputs_complete.json")
            shared.verify_hashes(checkpoint["output_sha256"])
            cohort = common.read(task / "cohort.json")
            inputs += [task / "inputs_complete.json", task / "cohort.json", task / "markers.tsv",
                       Path(cohort["reference"]["vcf"]), Path(cohort["pure_train"]["vcf"]),
                       Path(cohort["flare_jar"]), Path(cohort["genetic_map"]), common.BEAGLE]
            inputs += [Path(b["vcf"]) for b in cohort["pedigree_train"]["batches"]]
        span = cohort["chromosome_spans"][str(chrom)]
        values = {r["sample_id"]: [float(r["source_group"] == s) for s in SOURCES] for r in pure}
        for batch in selection["pedigree_train"]:
            ids = {r["sample_id"] for r in people if r.get("batch") == batch["batch"]}
            values.update(chromosome_truth(batch["segments_tsv"], ids, chrom, span))
        if set(values) != {r["sample_id"] for r in people}:
            raise ValueError("Incomplete training truth coverage")
        truth[str(chrom)] = values
    paths = [Path(__file__), Path(diagnostic.__file__), Path(dense.__file__), Path(common.__file__),
             Path(native.__file__), Path(shared.__file__),
             Path(__file__).with_name("screen_model_families.py"),
             Path(__file__).with_name("verify_dense_observed.py"),
             Path(__file__).with_name("validate_phasefree_arrays.py"),
             PROJECT / "run_flare_reference_density.sh"]
    output_inputs = out / "inputs"
    output_inputs.mkdir()
    by_id = {r["sample_id"]: r for r in reference}
    representation_counts = {}
    for name, ids in representations.items():
        pd.DataFrame([(identity, by_id[identity]["source_group"]) for identity in ids]).to_csv(
            output_inputs / f"{name}.panel", sep="\t", header=False, index=False)
        representation_counts[name] = dict(Counter(by_id[identity]["source_group"] for identity in ids))
    native.write_json(out / "selection.json", dict(people=people, references=representations, truth=truth))
    inputs += list(output_inputs.glob("*.panel")) + [out / "selection.json"]
    native.write_json(out / "experiment.json", dict(
        protocol=protocol, n_pure=240, n_descendants=76, n_families=25,
        reference_counts=representation_counts, jar_audit=jar_audit,
        input_sha256=shared.fingerprint_files(inputs), implementation_sha256=shared.fingerprint_files(paths),
        evaluation_people_included=False, independent_final_test=False))
    native.write_json(out / "state.json", dict(status="prepared_training_only"))


def smoke_people(selection, seed):
    pure = min((r for r in selection["people"] if r["source_group"] == "SAS"),
               key=lambda r: rank(r["sample_id"], seed, "smoke"))
    pedigrees = [r for r in selection["people"] if r["track"] == "pedigree"]
    descendant = min(pedigrees, key=lambda r: (
        -sum(x > 0 for x in selection["truth"]["22"][r["sample_id"]]),
        rank(r["sample_id"], seed, "smoke")))
    return [pure, descendant]


def evaluate(out, panel, chrom, smoke):
    config, selection = common.read(out / "experiment.json"), common.read(out / "selection.json")
    protocol = config["protocol"]
    people = smoke_people(selection, protocol["selection_seed"]) if smoke else selection["people"]
    task = BASE / f"tasks/{'frederik' if panel == 'frederik_thin' else panel}_chr{chrom}"
    cohort = common.read(task / "cohort.json")
    keys = markers(task)
    directory = out / ("smoke" if smoke else "tasks") / f"{panel}_chr{chrom}"
    directory.mkdir(parents=True)
    if panel == "frederik_thin":
        n = len(markers(BASE / f"tasks/lonni_chr{chrom}"))
        indices = thin_indices([key[1] for key in keys], n)
        keys = [keys[i] for i in indices]
    frame = pd.DataFrame(keys, columns=["chrom", "pos", "ref", "alt"])
    frame.to_csv(directory / "markers.tsv", sep="\t", index=False)
    inherited_cohort = {**cohort, "reference": dict(cohort["reference"])}
    refs = {}
    for representation, ids in selection["references"].items():
        work = directory / representation
        work.mkdir()
        vcf = (cohort["reference"]["vcf"] if representation == "full" and panel != "frederik_thin"
               else common.subset(cohort["reference"]["vcf"], ids, keys, work, "reference"))
        refs[representation] = {**cohort, "reference": dict(
            vcf=vcf, broad_panel=str(out / "inputs" / f"{representation}.panel"))}
    inherited_cohort["reference"] = refs["full"]["reference"]
    span = cohort["chromosome_spans"][str(chrom)]
    cases = []
    sources = {r["sample_id"]: cohort["pure_train"]["vcf"] for r in people if r["track"] == "pure"}
    batches = {b["batch"]: b["vcf"] for b in cohort["pedigree_train"]["batches"]}
    sources.update({r["sample_id"]: batches[r["batch"]] for r in people if r["track"] == "pedigree"})
    ids_by_path = defaultdict(list)
    for identity, path in sources.items():
        ids_by_path[path].append(identity)
    metadata = {r["sample_id"]: r for r in people}
    for path, ids in ids_by_path.items():
        all_keys = markers(task)
        matrix = native.load_haplotypes(path, ids, chrom, all_keys).reshape(len(ids), 2, -1)
        if panel == "frederik_thin":
            matrix = matrix[:, :, indices]
        for i, identity in enumerate(ids):
            work = directory / "cases" / identity
            work.mkdir(parents=True)
            case = dense.phase_case(work, identity, chrom, keys, matrix[i], inherited_cohort)
            n, switches, pairs = check_case((case, keys))
            if n != len(keys):
                raise ValueError("Independent phase-only GT verification failed")
            case.update(track=metadata[identity]["track"], family_id=metadata[identity]["family_id"],
                        phase_switches=switches, phase_pairs=pairs)
            cases.append(case)
            native.write_json(work / "verified.json", dict(
                **case, output_sha256=shared.fingerprint_files([
                    Path(case["inherited_phase"]), Path(case["array_phase_only"]),
                    work / "typed.unphased.vcf.gz"])))
            if len(cases) % 20 == 0:
                print(f"{panel} chr{chrom}: phased/independently verified {len(cases)}/{len(people)} training people", flush=True)
    native.write_json(directory / "cases.json", cases)
    templates = {}
    template_case = next(case for case in cases if case["track"] == "pure")
    for representation, current in refs.items():
        prefix = directory / representation / "initialization"
        shared.checked(command(current, template_case["inherited_phase"], prefix, False),
                       prefix.parent, "initialization")
        model = Path(str(prefix) + ".fixed.model")
        write_model(model, read_model(str(prefix) + ".model"))
        templates[representation] = model
    rows = []
    for candidate, representation, em in candidates(protocol):
        work = directory / candidate
        work.mkdir()
        for case in cases:
            identity = case["sample_id"]
            truth = np.asarray(selection["truth"][str(chrom)][identity], float)
            for stage in protocol["stages"]:
                prefix = work / f"{identity}_{stage}"
                seconds = shared.checked(command(
                    refs[representation], case[stage], prefix, em, templates[representation]), work, prefix.name)
                model = read_model(str(prefix) + ".model")
                ids, p = physical_predictions(str(prefix) + ".anc.vcf.gz", model["ancestries"],
                                               {str(chrom): span}, frame)
                if ids != [identity] or not np.isfinite(p).all() or not np.isclose(p.sum(), 1):
                    raise ValueError("Invalid native FLARE sample/fraction axis")
                independent = local_physical(str(prefix) + ".anc.vcf.gz", model["ancestries"], span, keys)
                np.testing.assert_allclose(p[0], independent, rtol=0, atol=1e-10)
                if identity == template_case["sample_id"] and stage == "inherited_phase" and not em:
                    initial = directory / representation / "initialization.anc.vcf.gz"
                    implicit = local_physical(initial, model["ancestries"], span, keys)
                    np.testing.assert_allclose(p[0], implicit, rtol=0, atol=1e-10)
                row = dict(panel=panel, chromosome=chrom, candidate=candidate, stage=stage,
                           sample_id=identity, track=case["track"], family_id=case["family_id"],
                           fine=metadata[identity]["fine"], scenario=metadata[identity]["scenario"],
                           source_group=metadata[identity]["source_group"], n_markers=len(keys),
                           span_bp=span[1] - span[0], seconds=seconds,
                           tv=float(np.abs(p[0] - truth).sum() / 2),
                           max_component_error=float(np.abs(p[0] - truth).max()),
                           off_target_mass=float(p[0][truth == 0].sum()),
                           **{f"pred_{s}": float(p[0, i]) for i, s in enumerate(SOURCES)},
                           **{f"truth_{s}": float(truth[i]) for i, s in enumerate(SOURCES)})
                native.write_json(str(prefix) + ".verified.json", dict(
                    prediction=row, independent_local_aggregation=True,
                    output_sha256=shared.fingerprint_files([
                        Path(str(prefix) + ".anc.vcf.gz"), Path(str(prefix) + ".model")])))
                rows.append(row)
            if len(rows) % 80 == 0:
                print(f"{panel} chr{chrom}: {len(rows)} posterior predictions independently verified", flush=True)
        print(f"{panel} chr{chrom}: {candidate} complete", flush=True)
    expected = {(candidate, stage, r["sample_id"]) for candidate, _, _ in candidates(protocol)
                for stage in protocol["stages"] for r in people}
    if len(rows) != len(expected) or {
            (r["candidate"], r["stage"], r["sample_id"]) for r in rows} != expected:
        raise ValueError("Incomplete/duplicate training method coverage")
    pd.DataFrame(rows).to_csv(directory / "predictions.tsv", sep="\t", index=False)
    native.write_json(directory / "complete.json", dict(
        status="complete_training_only_block", n_people=len(people), n_predictions=len(rows),
        n_typed_gt_preserved=sum(case["n_typed_gt_preserved"] for case in cases),
        independent_local_aggregations=len(rows), no_changed_gt=True,
        output_sha256=shared.fingerprint_files([
            directory / "predictions.tsv", directory / "cases.json", directory / "markers.tsv"])))
    return str(directory)


def summarize(out, smoke):
    config = common.read(out / "experiment.json")["protocol"]
    blocks = [(panel, chrom) for chrom in ([22] if smoke else config["chromosomes"]) for panel in config["panels"]]
    frames = []
    for panel, chrom in blocks:
        task = out / ("smoke" if smoke else "tasks") / f"{panel}_chr{chrom}"
        shared.verify_hashes(common.read(task / "complete.json")["output_sha256"])
        frames.append(pd.read_csv(task / "predictions.tsv", sep="\t", float_precision="round_trip"))
    frame = pd.concat(frames, ignore_index=True)
    summaries, fine = [], []
    axes = ["panel", "chromosome", "candidate", "stage", "track"]
    for key, group in frame.groupby(axes):
        summaries.append(dict(zip(axes, key), n_people=len(group), mean_tv=group.tv.mean(),
                              p95_tv=group.tv.quantile(.95), max_tv=group.tv.max(),
                              all_components_within_10pp=(group.max_component_error <= .1 + 1e-12).mean(),
                              all_components_within_5pp=(group.max_component_error <= .05 + 1e-12).mean(),
                              off_target_mass=group.off_target_mass.mean(),
                              mean_max_component_error=group.max_component_error.mean(),
                              family_balanced_tv=group.groupby("family_id").tv.mean().mean(),
                              pure_source_balanced_tv=(group.groupby("source_group").tv.mean().mean()
                                                       if key[-1] == "pure" else np.nan),
                              inference_seconds=group.seconds.sum()))
    pure = frame[frame.track.eq("pure")]
    for key, group in pure.groupby(axes[:-1] + ["source_group", "fine"]):
        fine.append(dict(zip(axes[:-1] + ["source_group", "fine"], key),
                         n_people=len(group), mean_tv=group.tv.mean(), max_tv=group.tv.max()))
    directory = out / "smoke" if smoke else out
    frame.to_csv(directory / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(directory / "comparison.tsv", sep="\t", index=False)
    pd.DataFrame(fine).to_csv(directory / "fine_population_errors.tsv", sep="\t", index=False)
    native.write_json(directory / "completion.json", dict(
        status="complete_smoke" if smoke else "complete_training_only_three_chromosome_screen",
        n_unique_people=frame.sample_id.nunique(), n_predictions=len(frame),
        chromosomes=[22] if smoke else config["chromosomes"], whole_genome_score=False,
        independent_final_test=False, production_promotion=False,
        output_sha256=shared.fingerprint_files([
            directory / name for name in ("predictions.tsv", "comparison.tsv", "fine_population_errors.tsv")])))


def run(out, smoke):
    experiment = common.read(out / "experiment.json")
    shared.verify_hashes(experiment["input_sha256"])
    shared.verify_hashes(experiment["implementation_sha256"])
    protocol = experiment["protocol"]
    if not smoke:
        completed = common.read(out / "smoke/completion.json")
        if completed["status"] != "complete_smoke":
            raise ValueError("Verified real smoke required before full launch")
        shared.verify_hashes(completed["output_sha256"])
    native.write_json(out / "state.json", dict(status="running_smoke" if smoke else "running_training_only"))
    for chrom in ([22] if smoke else protocol["chromosomes"]):
        with ProcessPoolExecutor(max_workers=2) as pool:
            jobs = [pool.submit(evaluate, out, panel, chrom, smoke) for panel in protocol["panels"]]
            for job in as_completed(jobs):
                print(f"Verified block {job.result()}", flush=True)
    summarize(out, smoke)
    native.write_json(out / "state.json", dict(status="smoke_verified" if smoke else "complete_training_only"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "smoke", "run"), required=True)
    parser.add_argument("--outdir", type=Path, default=Path("/screen"))
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        try:
            if args.stage == "prepare":
                prepare(args.outdir)
            else:
                run(args.outdir, args.stage == "smoke")
        except Exception as error:
            native.write_json(args.outdir / "state.json", dict(status="failed", stage=args.stage, error=str(error)))
            raise
