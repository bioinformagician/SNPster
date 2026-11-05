"""Bounded GME mechanism screen and founder-grounded eight-source discrimination."""
import argparse
from collections import Counter, defaultdict
import copy
import csv
import shutil
from pathlib import Path

import numpy as np

import core
import models


def merged(source):
    return "EAS_SEA" if source in ("EAS", "SEA") else source


def write_rows(path, records, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def enriched_samples(original, previous, authority):
    original_by_id = {row["sample_id"]: row for row in original}
    for row in previous:
        if original_by_id.get(row["sample_id"]) != row:
            raise ValueError("Previous development subset differs from original role ownership")
    chosen = {row["sample_id"] for row in previous}
    for row in original:
        if (row["role"] in ("calibration", "validation") and row["track"] == "pure"
                and (row["source"] == "GME" or (
                    row["source"] == "EAS_SEA" and
                    (row["role"] == "validation" or authority[row["sample_id"]] == "SEA")))):
            chosen.add(row["sample_id"])
    result = []
    for row in original:
        if row["role"] == "test" or row["sample_id"] not in chosen:
            continue
        item = dict(row)
        if row["track"] == "pure":
            label = authority[row["sample_id"]]
            if merged(label) != row["source"]:
                raise ValueError("Original source disagrees with authoritative unmerged label")
            item["source"] = label
        result.append(item)
    if any(row["role"] == "test" for row in result):
        raise ValueError("Guided development must exclude test people")
    for role in ("reference", "train"):
        expected = [row for row in original if row["role"] == role]
        actual = [{**row, "source": merged(row["source"])} for row in result if row["role"] == role]
        if actual != expected:
            raise ValueError("Guided screen changes the original reference or training cohort")
    return result


def expand_inheritance(records, samples, original, authority, chrom, span):
    roles = {row["sample_id"]: row["role"] for row in original}
    by_id = {row["sample_id"]: row for row in samples}
    selected, tiles = [], defaultdict(list)
    for row in records:
        if row["sample_id"] not in by_id:
            continue
        identity, donor = row["sample_id"], row["donor_id"]
        sample = by_id[identity]
        if (row["chromosome"] != str(chrom) or donor not in authority
                or donor not in sample["donors"].split(",") or roles[donor] != sample["role"]
                or row["donor_haplotype"] not in ("0", "1")
                or merged(authority[donor]) != row["source"]):
            raise ValueError("Invalid founder label, chromosome, haplotype or cross-role inheritance")
        if sample["track"] == "pure" and (donor != identity or authority[donor] != sample["source"]):
            raise ValueError("Pure-source label differs from physical founder origin")
        item = {key: row[key] for key in ("sample_id", "chromosome", "haplotype", "start", "end", "source")}
        item["source"] = authority[donor]
        selected.append(item)
        tiles[(identity, row["haplotype"])].append((int(row["start"]), int(row["end"])))
    for identity in by_id:
        for hap in ("0", "1"):
            pieces = sorted(tiles[(identity, hap)])
            cursor = span[0]
            for start, end in pieces:
                if start != cursor or not start < end <= span[1]:
                    raise ValueError("Expanded founder truth has a gap, overlap or invalid span")
                cursor = end
            if cursor != span[1]:
                raise ValueError("Expanded founder truth lacks complete sample/haplotype coverage")
    return selected


def native_parameters(path):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    names = lines[0]
    return names, np.asarray([list(map(float, row)) for row in lines[4:4 + len(names)]]), float(lines[2][0])


def diagnose(old, original, previous, output):
    summary = core.read(old / "genome/results/summary/summary/summary.json")
    evidence = core.read(old / "genome/results/summary/summary/prediction_evidence.json")
    receipt = core.read(old / "genome/results/independent_verification.json")
    counts = {}
    for role in ("reference", "train", "calibration", "validation"):
        selected = [row for row in original if row["role"] == role and row["track"] == "pure"]
        counts[role] = {"sources": dict(Counter(row["source"] for row in selected)),
                        "GME_populations": dict(Counter(row["fine"] for row in selected if row["source"] == "GME"))}
    gme = []
    for row in previous:
        if row["source"] != "GME" or row["role"] not in ("calibration", "validation"):
            continue
        scores = {}
        for model in ("flare_full_em_g010", "flare_full_fixed_g002", "flare_cap16_em_g010"):
            values = [item["tv_percentage_points"] for item in evidence
                      if item["sample_id"] == row["sample_id"] and item["model"] == model
                      and item["interpreter"] == "native"]
            if len(values) != 2:
                raise ValueError("Diagnosis lacks both observed masks for a GME sample")
            scores[model] = float(np.mean(values))
        copying = []
        for chrom in (20, 22):
            for route in ("harmonized_lonni", "harmonized_frederik"):
                path = old / f"genome/results/predictions/{route}/chr{chrom}/flare_full_em_g010/prediction/{row['sample_id']}.model"
                if core.digest(path) != receipt["files"][str(path.relative_to(old / "genome/results"))]:
                    raise ValueError("Native parameter diagnosis differs from independently audited evidence")
                names, matrix, generations = native_parameters(path)
                copying.append({"chromosome": chrom, "route": route,
                                "GME_copying_row": dict(zip(names, matrix[names.index("GME")].tolist())),
                                "target_fitted_generations": generations})
        gme.append({**row, "native_mean_tv_pp": scores, "native_parameter_diagnostics": copying})
    core.write(output, {
        "counts": counts, "GME_development_people": gme,
        "old_summary_sha256": core.digest(old / "genome/results/summary/summary/summary.json"),
        "diagnosis": [
            "GME has 65 reference people spanning four populations, versus 12 OCE and 21 NAM; not an absent source.",
            "GME underprediction affects Bedouin, Druze and Palestinian development people and both masks.",
            "Fixed-generation 2 reduces pure-GME error; fine-population cap16 worsens it.",
            "Native reference-estimated GME copying support is shared across targets and weak relative to EUR/SAS copying.",
            "Counts, related reference structure, shared variation, finite IBS states and target EM may all contribute; causality is not isolated.",
            "Operational 100%-GME truth is a population annotation, not measured individual biological admixture.",
        ],
        "sweep_hypotheses": {
            "full_em100": "eight-source version of the previous overall winner",
            "full_fixed2": "previous challenger; improved pure GME but worse G1 admixture",
            "full_em280": "more IBS states without changing reference membership",
            "full_em100_updatep": "allow target EM to adapt reference-copying parameters",
            "balanced65_em100": "reduce large-source count imbalance, retaining all 65 GME; no duplicates or held-out references",
        },
        "test_used_for_diagnosis_or_tuning": False, "previous_winner_tv_pp": min(
            row["validation_mean_person_tv_pp"] for row in summary["validation_leaderboard"]),
    })


def plan(specification, output, durable):
    spec = core.read(specification)
    if (spec["version"] != 1 or spec["chromosomes"] != [20, 22]
            or spec["sources"] != ["EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"]
            or len(spec["models"]) != 5 or any(row["family"] != "flare" for row in spec["models"])):
        raise ValueError("This guided protocol requires exactly five FLARE settings, eight sources and chr20/22")
    old = Path(spec["old_run"])
    verified = core.read(old / "genome/results/independent_verification.json")
    if not verified["verified"]:
        raise ValueError("Guided follow-up requires the independently audited completed predecessor")
    from evaluation import verify
    verify(old / "genome/results/summary/summary")
    for relative in ("summary/summary/summary.json", "summary/summary/prediction_evidence.json",
                     "frozen/frozen/config.json"):
        if core.digest(old / "genome/results" / relative) != verified["files"][relative]:
            raise ValueError("Predecessor summary or frozen configuration differs from independent audit")
    config = core.read(old / "experiment.json")
    original = core.rows(config["samples"])
    previous = core.rows(old / "planning/plan/samples.tsv")
    release = Path(config["samples"]).parents[3]
    resource_receipt = core.read(release / "verification/independent.json")
    if not resource_receipt["verified"]:
        raise ValueError("Original standardized founder-transmission resource is not independently verified")
    for path in (release / "inventory/inventory.json", Path(config["samples"]), Path(config["segments"])):
        if core.digest(path) != resource_receipt["files_sha256"][str(path.relative_to(release))]:
            raise ValueError("Original source roles, physical truth or provenance changed")
    provenance = core.read(release / "inventory/inventory.json")
    authority_path = Path(spec["authority"])
    if core.digest(authority_path) != provenance["input_sha256"][str(authority_path)]:
        raise ValueError("Authoritative source labels changed since the standardized resource was frozen")
    labels = [row for row in core.rows(authority_path) if row["role"] == "source"]
    authority = {row["sample"]: row["source_group"] for row in labels}
    if len(authority) != len(labels):
        raise ValueError("Authoritative source-label identities are ambiguous")
    samples = enriched_samples(original, previous, authority)
    core.validate_roles(samples, spec["sources"], development_only=True)
    queries = [row for row in samples if row["role"] in ("calibration", "validation")]
    calls = len(queries) * len(config["routes"]) * len(spec["chromosomes"]) * len(spec["models"])
    if calls > spec["maximum_singleton_calls"] or spec["maximum_singleton_calls"] > 4000:
        raise ValueError("Guided screen exceeds its 4,000-singleton hard scope gate")
    output.mkdir()
    write_rows(output / "samples.tsv", samples, list(samples[0]))
    diagnose(old, original, previous, output / "diagnosis.json")
    segments, phase_assets, input_hashes = [], {}, {str(authority_path): core.digest(authority_path)}
    original_segments = defaultdict(list)
    sample_ids = {row["sample_id"] for row in samples}
    for row in core.rows(config["segments"]):
        if row["chromosome"] in map(str, spec["chromosomes"]) and row["sample_id"] in sample_ids:
            original_segments[(row["sample_id"], row["chromosome"], row["haplotype"])].append(row)
    for chrom in spec["chromosomes"]:
        inheritance = release / f"panels/observed/genotypes/all/chr{chrom}/master/inheritance.tsv"
        if core.digest(inheritance) != resource_receipt["files_sha256"][str(inheritance.relative_to(release))]:
            raise ValueError("Founder-haplotype transmissions changed after independent genotype verification")
        reconstructed = expand_inheritance(core.rows(inheritance), samples, original, authority, chrom,
                                           config["spans"][str(chrom)])
        input_hashes[str(inheritance)] = core.digest(inheritance)
        for row in reconstructed:
            if {**row, "source": merged(row["source"])} not in original_segments[
                    (row["sample_id"], str(chrom), row["haplotype"])]:
                raise ValueError("Eight-source physical truth does not reproduce an original merged segment")
        expected_count = sum(len(rows) for key, rows in original_segments.items() if key[1] == str(chrom))
        if len(reconstructed) != expected_count:
            raise ValueError("Expanded physical-truth segment coverage differs from the original release")
        segments.extend(reconstructed)
        for route in config["routes"]:
            root = old / "genome/results"
            panel = root / f"panels/{route['id']}/chr{chrom}/panel"
            phased = root / f"phased/{route['id']}/chr{chrom}/phased"
            required = [panel / "panel.json", phased / "phase.json", phased / "query.npy"]
            hashes = {}
            for path in required:
                relative = str(path.relative_to(root))
                expected = verified["files"][relative]
                if core.digest(path) != expected:
                    raise ValueError("Original reusable preprocessing differs from independent audit")
                hashes[str(path)] = expected
            phase_assets[f"{route['id']}/{chrom}"] = {"panel": str(panel), "phased": str(phased), "hashes": hashes}
    write_rows(output / "segments.tsv", segments,
               ["sample_id", "chromosome", "haplotype", "start", "end", "source"])
    result = copy.deepcopy(config)
    result.update(version=2, mode="development_screen", independent_test=False, sources=spec["sources"],
                  chromosomes=spec["chromosomes"], samples=str(durable / "samples.tsv"),
                  segments=str(durable / "segments.tsv"), models=spec["models"],
                  interpretations=[{"id": "native", "kind": "native", "params": {}}])
    for key in ("reference", "training", "maps", "spans"):
        result[key] = {str(chrom): result[key][str(chrom)] for chrom in spec["chromosomes"]}
    for route in result["routes"]:
        route["vcfs"] = {str(chrom): route["vcfs"][str(chrom)] for chrom in spec["chromosomes"]}
    core.validate_config(result)
    core.write(output / "config.json", result)
    core.write(output / "phase-reuse.json", {"assets": phase_assets, "seed": result["seed"],
                                            "reference_filters": result["reference_filters"]})
    write_rows(output / "source-labels.tsv", [
        {"sample_id": row["sample_id"], "role": row["role"], "source": row["source"], "fine": row["fine"]}
        for row in samples if row["track"] == "pure"], ["sample_id", "role", "source", "fine"])
    reference_rows = [row for row in samples if row["role"] == "reference"]
    reference_meta = {"reference": reference_rows, "sources": spec["sources"]}
    core.write(output / "provenance.json", {
        "source_specification_sha256": core.digest(specification), "original_run": str(old),
        "input_sha256": input_hashes, "truth_policy": "original founder donor and haplotype at each physical interval",
        "authority": str(authority_path),
        "merge_reconstructs_original_truth_exactly": True,
        "reference_training_people_and_genotypes_changed": False,
        "new_reference_source_labels": True, "old_merged_model_weights_reused_as_eight_source": False,
        "source_taxonomy": "existing frozen authority; CDX/Dai EAS; KHV/Cambodian/Lahu SEA; no geographic relabelling",
        "roles": dict(Counter(row["role"] for row in samples)),
        "pure_sources_by_role": {role: dict(Counter(row["source"] for row in queries
            if row["role"] == role and row["track"] == "pure")) for role in ("calibration", "validation")},
        "family_components": {role: len({row["family_id"] for row in queries if row["role"] == role})
                              for role in ("calibration", "validation")},
        "old_matched_development_ids": [row["sample_id"] for row in previous
                                       if row["role"] in ("calibration", "validation")],
        "singleton_calls": calls, "maximum_singleton_calls": spec["maximum_singleton_calls"],
        "reference_source_counts_by_model": {model["id"]: dict(Counter(
            reference_rows[index]["source"] for index in models.reference_subset(reference_meta, model["params"])))
            for model in spec["models"]},
        "automatic_genome_expansion": False, "test_evaluated": False,
        "partial_chromosome_development_evidence_only": True,
        "SEA_calibration_warning": "Only three pure SEA calibration people; native output only, no learned decoder sweep",
    })
    print(f"Guided plan: {len(queries)} development queries; {len(spec['models'])} FLARE settings; "
          f"chr{spec['chromosomes']}; {calls} singleton calls; no test or automatic genome expansion")


def phase_chunk(panel, reuse_path, start, end, output):
    meta = core.read(panel / "panel.json")
    reuse = core.read(reuse_path)
    if not 0 <= start < end <= len(meta["queries"]):
        raise ValueError("Guided phase range is outside the declared query axis")
    asset = reuse["assets"][f"{meta['route']['id']}/{meta['chromosome']}"]
    for filename, expected in asset["hashes"].items():
        if core.digest(filename) != expected:
            raise ValueError(f"Original phase reuse evidence changed: {filename}")
    old = core.read(Path(asset["panel"]) / "panel.json")
    phase = core.read(Path(asset["phased"]) / "phase.json")
    old_reference = [{**row, "source": merged(row["source"])} for row in meta["reference"]]
    source_reference = [{**row, "source": merged(row["source"])} for row in old["reference"]]
    if (old_reference != source_reference or meta["keys"] != old["keys"]
            or meta["config"]["seed"] != reuse["seed"]
            or meta["config"]["reference_filters"] != reuse["reference_filters"]
            or meta["route"]["vcfs"][str(meta["chromosome"])] != old["route"]["vcfs"][str(meta["chromosome"])]
            or meta["config"]["maps"][str(meta["chromosome"])] != old["config"]["maps"][str(meta["chromosome"])]
            or meta["config"]["tools"]["beagle_jar"] != old["config"]["tools"]["beagle_jar"]
            or not phase["phase_only"] or phase["changed_diploid_calls"] != 0
            or phase["query_sha256"] != core.digest(Path(asset["phased"]) / "query.npy")):
        raise ValueError("Guided phase reuse changes scientific genotypes, reference membership, keys or seed")
    core.verify_tools(meta)
    old_index = {row["sample_id"]: index for index, row in enumerate(old["queries"])}
    before = np.load(panel / "data.npz")["query"][start:end]
    previous = np.load(Path(asset["phased"]) / "query.npy", mmap_mode="r")
    samples = meta["queries"][start:end]
    output.mkdir()
    calls, reused = [], []
    for index, sample in enumerate(samples):
        identity = sample["sample_id"]
        prefix = output / f"{identity}.phase"
        if identity in old_index:
            source_row = old["queries"][old_index[identity]]
            if {**sample, "source": merged(sample["source"])} != {**source_row, "source": merged(source_row["source"])}:
                raise ValueError("Reused phasing has different sample, truth label or role ownership")
            value = previous[old_index[identity]]
            reused.append(identity)
            source = Path(asset["phased"])
            for path in source.glob(f"{identity}.*"):
                if path.is_file() and not path.is_symlink():
                    shutil.copyfile(path, output / path.name)
        else:
            target = output / f"{identity}.vcf.gz"
            core.write_vcf(target, meta["keys"], [identity], before[index:index + 1], phased=False)
            config = meta["config"]
            core.run(["java", core.java_heap(), "-jar", config["tools"]["beagle_jar"], f"gt={target}",
                      f"ref={panel / 'reference.vcf.gz'}", f"map={config['maps'][str(meta['chromosome'])]}",
                      f"out={prefix}", "impute=false", "nthreads=1", f"seed={config['seed']}"], prefix)
            keys, native = core.vcf(str(prefix) + ".vcf.gz", meta["chromosome"], [identity], True)
            if keys != [tuple(key) for key in meta["keys"]]:
                raise ValueError("Guided phase-only execution changed exact allele keys")
            value = native[0]
        if not np.array_equal(value.sum(axis=0), before[index].sum(axis=0)):
            raise ValueError("Guided phase-only execution changed observed diploid genotypes")
        calls.append(value)
    np.save(output / "query.npy", np.asarray(calls))
    core.write(output / "chunk.json", {
        "start": start, "end": end, "total": len(meta["queries"]), "samples": samples,
        "panel_sha256": core.digest(panel / "panel.json"), "data_sha256": core.digest(panel / "data.npz"),
        "query_sha256": core.digest(output / "query.npy"), "reuse_manifest_sha256": core.digest(reuse_path),
        "reused_singletons": reused,
        "artifacts": {path.name: core.digest(path) for path in output.iterdir() if path.is_file() and path.name != "query.npy"},
    })


def collect_phase(panel, folders, output):
    meta = core.read(panel / "panel.json")
    records = sorted(((path, core.read(path / "chunk.json")) for path in folders), key=lambda x: x[1]["start"])
    if not records:
        raise ValueError("Require complete guided phase chunks")
    output.mkdir()
    before = np.load(panel / "data.npz")["query"]
    cursor, values, reused = 0, [], []
    shared = {key: records[0][1][key] for key in ("panel_sha256", "data_sha256", "reuse_manifest_sha256")}
    if shared["panel_sha256"] != core.digest(panel / "panel.json") or shared["data_sha256"] != core.digest(panel / "data.npz"):
        raise ValueError("Phase chunks belong to a different genotype panel")
    for folder, record in records:
        start, end = record["start"], record["end"]
        if (start != cursor or not start < end <= len(meta["queries"]) or record["total"] != len(meta["queries"])
                or record["samples"] != meta["queries"][start:end]
                or any(record[key] != value for key, value in shared.items())
                or core.digest(folder / "query.npy") != record["query_sha256"]):
            raise ValueError("Incomplete, overlapping, changed or mismatched guided phase chunks")
        value = np.load(folder / "query.npy")
        if value.shape != before[start:end].shape or not np.array_equal(value.sum(axis=1), before[start:end].sum(axis=1)):
            raise ValueError("Collected phase chunks changed diploid genotype calls")
        for name, expected in record["artifacts"].items():
            if Path(name).name != name or (output / name).exists() or core.digest(folder / name) != expected:
                raise ValueError("Invalid, duplicated or changed native phase evidence")
            shutil.copyfile(folder / name, output / name)
        values.append(value)
        reused.extend(record["reused_singletons"])
        cursor = end
    if cursor != len(meta["queries"]):
        raise ValueError("Guided phase chunks have incomplete person coverage")
    np.save(output / "query.npy", np.concatenate(values))
    core.write(output / "phase.json", {"stage": "harmonized", "phase_only": True,
        "added_variants": 0, "changed_diploid_calls": 0, "query_sha256": core.digest(output / "query.npy"),
        "reused_singletons": reused, "query_chunks": len(records), **shared})


def report(results, design, output):
    from audit import verify_manifest
    verification = verify_manifest(results)
    if not verification["verified"]:
        raise ValueError("Guided report requires independently audited predictions")
    config = core.read(results / "frozen/frozen/config.json")
    provenance = core.read(design / "provenance.json")
    evidence = core.read(results / "summary/summary/prediction_evidence.json")
    samples = {row["sample_id"]: row for row in core.rows(config["samples"])}
    for path, expected in provenance["input_sha256"].items():
        if core.digest(path) != expected:
            raise ValueError("Authoritative label or founder-haplotype truth changed before reporting")
    if any(row["role"] == "test" for row in evidence):
        raise ValueError("Guided report must not include test results")
    authority = {row["sample"]: row["source_group"] for row in core.rows(provenance["authority"])
                 if row["role"] == "source"}
    denominator = 2 * sum(config["spans"][str(chrom)][1] - config["spans"][str(chrom)][0]
                          for chrom in config["chromosomes"])
    independent_truth = defaultdict(lambda: np.zeros(len(config["sources"])))
    query_ids = [identity for identity, row in samples.items() if row["role"] in ("calibration", "validation")]
    for path in provenance["input_sha256"]:
        if not path.endswith("/inheritance.tsv"):
            continue
        for row in core.rows(path):
            if row["sample_id"] in query_ids:
                source = authority[row["donor_id"]]
                independent_truth[row["sample_id"]][config["sources"].index(source)] += (
                    int(row["end"]) - int(row["start"])) / denominator
    expected_truth = core.probability(np.stack([independent_truth[identity] for identity in query_ids]))
    for route in config["routes"]:
        first_model = config["models"][0]["id"]
        aggregate = results / f"aggregated/{route['id']}/{first_model}/aggregate"
        axis = core.read(aggregate / "aggregate.json")["samples"]
        if [row["sample_id"] for row in axis] != query_ids:
            raise ValueError("Independent founder truth has a different query axis")
        np.testing.assert_allclose(np.load(aggregate / "aggregate.npz")["truth"], expected_truth, rtol=0, atol=1e-10)
    old = Path(provenance["original_run"])
    old_receipt = core.read(old / "genome/results/independent_verification.json")
    matched_ids = {row["sample_id"] for row in samples.values() if row["role"] == "validation"
                   and row["sample_id"] in provenance["old_matched_development_ids"]}
    old_baselines = {}
    for model in ("flare_full_em_g010", "flare_full_fixed_g002"):
        errors = defaultdict(list)
        for route in config["routes"]:
            weighted_values, weighted_truth, length = [], [], 0
            for chrom in config["chromosomes"]:
                folder = old / f"genome/results/predictions/{route['id']}/chr{chrom}/{model}/prediction"
                for name in ("prediction.json", "predictions.npz"):
                    path = folder / name
                    if core.digest(path) != old_receipt["files"][str(path.relative_to(old / "genome/results"))]:
                        raise ValueError("Matched previous-model prediction evidence changed")
                meta = core.read(folder / "prediction.json")
                data = np.load(folder / "predictions.npz")
                span = meta["span"][1] - meta["span"][0]
                weighted_values.append(span * data["values"])
                weighted_truth.append(span * data["truth"])
                length += span
            values, truth = np.sum(weighted_values, axis=0) / length, np.sum(weighted_truth, axis=0) / length
            for index, person in enumerate(meta["samples"]):
                if person["sample_id"] in matched_ids:
                    errors[person["sample_id"]].append(float(50 * np.abs(values[index] - truth[index]).sum()))
        if set(errors) != matched_ids or any(len(values) != len(config["routes"]) for values in errors.values()):
            raise ValueError("Matched old-cohort baseline lacks complete person/mask coverage")
        old_baselines[model] = {"people": len(errors),
            "mean_tv_pp": float(np.mean([np.mean(values) for values in errors.values()]))}
    output.mkdir()
    rows = []
    for model in config["models"]:
        selected = [row for row in evidence if row["role"] == "validation" and row["model"] == model["id"]]
        people = defaultdict(list)
        for row in selected:
            people[row["sample_id"]].append(row)
        metrics = {}
        for source in config["sources"]:
            pure = [row for row in selected if samples[row["sample_id"]]["source"] == source
                    and row["track"] == "pure"]
            index = config["sources"].index(source)
            metrics[source] = {
                "pure_people": len({row["sample_id"] for row in pure}),
                "mean_tv_pp": float(np.mean([row["tv_percentage_points"] for row in pure])) if pure else None,
                "mean_predicted_own_source_pp": float(100 * np.mean([row["predictions"][index] for row in pure])) if pure else None,
                "mean_confusion_pp": {other: float(100 * np.mean([row["predictions"][j] for row in pure]))
                                      for j, other in enumerate(config["sources"])} if pure else {},
            }
        eas, sea = config["sources"].index("EAS"), config["sources"].index("SEA")
        split_people = [row for row in selected if row["track"] == "pure"
                        and samples[row["sample_id"]]["source"] in ("EAS", "SEA")]
        decisions = defaultdict(list)
        class_decisions, full_class_decisions = defaultdict(list), defaultdict(list)
        for row in split_people:
            true_source = samples[row["sample_id"]]["source"]
            correct = row["predictions"][eas] > row["predictions"][sea] if true_source == "EAS" else row["predictions"][sea] > row["predictions"][eas]
            decisions[row["sample_id"]].append(bool(correct))
            class_decisions[true_source].append(bool(correct))
            values = np.asarray(row["predictions"])
            index = config["sources"].index(true_source)
            full_class_decisions[true_source].append(bool(values[index] > np.delete(values, index).max()))
        split_ancestry = [row for row in selected if row["truth"][eas] + row["truth"][sea] >= .02]
        matched = [row for row in selected if row["sample_id"] in provenance["old_matched_development_ids"]]
        def merged_tv(row):
            errors = np.asarray(row["predictions"]) - np.asarray(row["truth"])
            return float(50 * (np.abs(np.delete(errors, [eas, sea])).sum() + abs(errors[eas] + errors[sea])))
        matched_people = defaultdict(list)
        for row in matched:
            matched_people[row["sample_id"]].append(merged_tv(row))
        rows.append({
            "model": model["id"], "sources": metrics,
            "mean_tv_pp": float(np.mean([np.mean([row["tv_percentage_points"] for row in items])
                                        for items in people.values()])),
            "EAS_SEA_pure_per_route_fraction_correct": float(np.mean([value for values in decisions.values() for value in values])),
            "EAS_SEA_pure_fraction_correct_by_source": {source: float(np.mean(values))
                                                        for source, values in class_decisions.items()},
            "EAS_SEA_pure_balanced_fraction_correct": float(np.mean([
                np.mean(values) for values in class_decisions.values()])),
            "pure_EAS_SEA_full_eight_class_fraction_correct_by_source": {
                source: float(np.mean(values)) for source, values in full_class_decisions.items()},
            "true_EAS_SEA_at_least_2pct": {
                "people": len({row["sample_id"] for row in split_ancestry}),
                "EAS_mae_pp": float(100 * np.mean([abs(row["predictions"][eas] - row["truth"][eas])
                                                 for row in split_ancestry])) if split_ancestry else None,
                "SEA_mae_pp": float(100 * np.mean([abs(row["predictions"][sea] - row["truth"][sea])
                                                 for row in split_ancestry])) if split_ancestry else None,
            },
            "matched_old_development_validation_people": len(matched_people),
            "matched_old_development_merged_seven_mean_tv_pp": float(np.mean([
                np.mean(values) for values in matched_people.values()])) if matched_people else None,
            "EAS_SEA_pure_people": len(decisions),
            "discrimination_note": "EAS-versus-SEA probability comparison; ties fail; not a calibrated classification threshold",
        })
    summary = core.read(results / "summary/summary/summary.json")
    winner = min(summary["validation_leaderboard"], key=lambda row: (row["validation_mean_person_tv_pp"], row["model"]))
    core.write(output / "source-discrimination.json", {
        "models": rows, "overall_validation_selected_model": winner["model"],
        "previous_winning_configuration_eight_source_version": "full_em100",
        "partial_chromosome_development_only": True, "chromosomes": config["chromosomes"],
        "test_evaluated": False, "production_promotion": False,
        "not_the_previous_merged_model": True, "source_taxonomy": provenance["source_taxonomy"],
        "not_genome_wide_or_independent_test_accuracy": True,
        "independent_founder_truth_reconstruction": True, "independent_truth_people_checked": len(query_ids),
        "previous_seven_source_baselines_same_chromosomes_and_matched_people": old_baselines,
        "selection_note": "Overall validation TV on a deliberately GME/EAS/SEA-enriched development cohort; "
                          "matched old-cohort merged-seven TV is reported separately",
    })
    tabular = []
    for row in rows:
        for source, values in row["sources"].items():
            tabular.append({"model": row["model"], "source": source, "people": values["pure_people"],
                            "mean_tv_pp": values["mean_tv_pp"],
                            "mean_predicted_own_source_pp": values["mean_predicted_own_source_pp"],
                            **{f"predicted_{name}_pp": value for name, value in values["mean_confusion_pp"].items()}})
    write_rows(output / "pure-source-confusion.tsv", tabular, list(tabular[0]))
    core.write(output / "complete.json", {"summary_sha256": verification["summary_sha256"],
        "files": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
        "test_evaluated": False})
    print("Guided source-discrimination report passed; development only; no test or production promotion")


def cleanup(root, output):
    from audit import verify_manifest
    root = root.resolve()
    if (root / "guided-followup").read_text().strip() != "true" or not (root / "source.sha256").is_file():
        raise ValueError("Cleanup requires the exact dedicated guided run")
    receipt = verify_manifest(root / "screen/results")
    report_root = root / "reporting/discrimination"
    completion = core.read(report_root / "complete.json")
    if not receipt["verified"] or completion["summary_sha256"] != receipt["summary_sha256"]:
        raise ValueError("Cannot clean a guided run without its independent audit and matching final report")
    for name, expected in completion["files"].items():
        if core.digest(report_root / name) != expected:
            raise ValueError("Guided report changed before successful-cache cleanup")
    workdirs = [root / stage / "work" for stage in ("planning", "screen", "reporting")]
    if any(path.is_symlink() or path.resolve() != path for path in workdirs):
        raise ValueError("Refuse redirected guided work-cache directories")
    removed = []
    for path in workdirs:
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(str(path))
    verify_manifest(root / "screen/results")
    for name, expected in completion["files"].items():
        if core.digest(report_root / name) != expected:
            raise ValueError("Guided report changed after successful-cache cleanup")
    core.write(output, {"removed_exact_successful_caches": removed, "durable_outputs_reverified": True,
                        "predecessor_and_shared_inputs_untouched": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("plan", "phase-chunk", "collect-phase", "report", "cleanup"))
    for name in ("config", "panel", "reuse", "design", "results", "output", "durable"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--folders", nargs="+", type=Path)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    args = parser.parse_args()
    if args.stage == "plan":
        plan(args.config, args.output, args.durable)
    elif args.stage == "phase-chunk":
        phase_chunk(args.panel, args.reuse, args.start, args.end, args.output)
    elif args.stage == "collect-phase":
        collect_phase(args.panel, args.folders, args.output)
    elif args.stage == "report":
        report(args.results, args.design, args.output)
    else:
        cleanup(args.results, args.output)
