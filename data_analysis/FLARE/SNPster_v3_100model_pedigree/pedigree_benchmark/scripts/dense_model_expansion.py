#!/usr/bin/env python3
"""Model-only expansion; leave the active ten-model producers immutable."""
import argparse
from collections import Counter, defaultdict
import copy
import itertools
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from audit_production_ancestry import SOURCES, digest, write_json
import model_family_gnomix as native
from production_matched_models import physical_flare, validate_vector
from production_matched_preprocessing import PROJECT, read, read_call_matrix
from sgdp_genetic_qc import run
from verify_sgdp_full_manifest import require_hashes


def signature(policy):
    family = policy["family"]
    if family == "flare":
        # With min-maf=0 and min-mac=0, pinned FLARE's array flag only changes logging.
        return (family, policy["reference"], policy["em"], policy.get("gen", 10))
    if family == "gnomix":
        return (family, policy["reference"], policy["window_size_cM"], policy.get("crf_c2", 1))
    if family == "phasefree":
        return (family, policy["reference"], policy["input"])
    raise ValueError(f"Unsupported model family: {family}")


def candidate_grid(baseline, plan):
    candidates, seen = {}, set()
    for name, policy in sorted(baseline.items(), key=lambda item: (not item[1].get("array", True), item[0])):
        key = signature(policy)
        if key not in seen:
            candidates[name] = copy.deepcopy(policy)
            seen.add(key)
    if len(baseline) != 10 or len(candidates) != 8:
        raise ValueError("Expected ten unchanged baseline outputs representing eight effective settings")
    flare = plan["flare"]
    for reference, array, em, generations in itertools.product(
            flare["reference_panels"], flare["array"], flare["em"], flare["generations"]):
        policy = dict(family="flare", reference=reference, array=array, em=em, gen=generations)
        key = signature(policy)
        if key not in seen:
            name = f"flare_{reference}_{'array' if array else 'nonarray'}_{'em' if em else 'fixed'}_g{generations:03d}"
            candidates[name] = policy
            seen.add(key)
    gnomix = plan["gnomix"]
    for window, c2 in itertools.product(gnomix["window_sizes_cM"], gnomix["crf_l2"]):
        policy = dict(family="gnomix", reference="full", window_size_cM=window, crf_c2=c2)
        key = signature(policy)
        if key not in seen:
            name = f"gnomix_fast_{round(window * 100):03d}_l2_{c2:g}".replace(".", "p")
            candidates[name] = policy
            seen.add(key)
    if (len(candidates) != plan["expected_configurations"] or len(candidates) > plan["maximum_configurations"]
            or Counter(policy["family"] for policy in candidates.values()) !=
            Counter(flare=64, gnomix=32, phasefree=2)):
        raise ValueError("Expanded grid violates the frozen family/count contract")
    return candidates


def freeze(source_experiment, plan_path, output):
    output.mkdir()
    source_cohort = source_experiment / "freeze/cohort/frozen/cohort.json"
    original, plan = read(source_cohort), read(plan_path)
    if (original["truth_design"] != "synthetic_reference_haplotype_development"
            or plan["flare"]["verified_default_generations"] != 10
            or plan["gnomix"]["verified_default_crf_l2"] != 1):
        raise ValueError("Expansion requires the approved synthetic split and verified pinned defaults")
    cohort = copy.deepcopy(original)
    cohort["protocol"]["candidates"] = candidate_grid(original["protocol"]["candidates"], plan)
    baseline_candidates = [name for name in original["protocol"]["candidates"]
                           if name in cohort["protocol"]["candidates"]]
    aliases = {name: next(identity for identity in baseline_candidates
                         if signature(original["protocol"]["candidates"][identity]) == signature(policy))
               for name, policy in original["protocol"]["candidates"].items() if name not in baseline_candidates}
    cohort["candidate_expansion"] = dict(
        policy=plan, policy_sha256=digest(plan_path), source_experiment=str(source_experiment),
        source_cohort=str(source_cohort), source_cohort_sha256=digest(source_cohort),
        baseline_candidates=baseline_candidates, baseline_aliases=aliases,
        source_controller=(source_experiment / "controller.name").read_text().strip(),
        candidate_selection_uses_prediction_errors=False)
    write_json(output / "cohort.json", cohort)
    write_json(output / "candidate_plan.json", dict(
        configurations=len(cohort["protocol"]["candidates"]),
        family_counts=dict(Counter(row["family"] for row in cohort["protocol"]["candidates"].values())),
        candidates=cohort["protocol"]["candidates"], source_cohort_sha256=digest(source_cohort),
        baseline_aliases=aliases,
        policy_sha256=digest(plan_path), external_people_reserved=len(cohort["reserved_external_ids"])))
    print("Frozen98 distinct settings:64 FLARE,32 fast Gnomix,2 NNLS; ten original outputs provide eight unique baselines.", flush=True)


def verify_same_population(original, expanded):
    first, second = copy.deepcopy(original), copy.deepcopy(expanded)
    first["protocol"].pop("candidates")
    second["protocol"].pop("candidates")
    second.pop("candidate_expansion")
    if first != second:
        raise ValueError("Model-only expansion changed the split, preprocessing, truth or scoring contract")


def verify_source(cohort_path, source_scope, output):
    cohort = read(cohort_path)
    expansion = cohort["candidate_expansion"]
    original_path = Path(expansion["source_cohort"])
    if digest(original_path) != expansion["source_cohort_sha256"]:
        raise ValueError("Original frozen donor/candidate epoch changed")
    verify_same_population(read(original_path), cohort)
    manifest_path = source_scope / "published_verification.json"
    manifest = read(manifest_path)
    post = read(source_scope / "post_cleanup_verification.json")
    if post["manifest_sha256"] != digest(manifest_path):
        raise ValueError("Original scope lacks matching post-cleanup verification")
    require_hashes(manifest["published_sha256"])
    require_hashes(manifest["source_sha256"])
    audit_path = source_scope / "summary/summary/independent_verification.json"
    audit = read(audit_path)
    smoke = source_scope.name == "smoke"
    expected_status = ("smoke_passed_not_a_whole_genome_leaderboard" if smoke
                       else "complete_verified_development_comparison")
    if (audit["status"] != expected_status or audit["model_configurations"] != 10
            or audit["cohort_sha256"] != digest(original_path)
            or audit["chromosomes"] != ([22] if smoke else list(range(1, 23)))
            or audit["masks"] != (["lonni"] if smoke else cohort["protocol"]["masks"])
            or audit["n_evaluation_people"] != 15):
        raise ValueError("Reuse requires the complete, independently verified original scope")
    for alias, canonical in expansion["baseline_aliases"].items():
        for chromosome in audit["chromosomes"]:
            frames = [pd.read_csv(source_scope / f"models/{name}/chr{chromosome}/model/predictions.tsv",
                                 sep="\t", float_precision="round_trip").set_index(
                                     ["mask", "chromosome", "sample_id"]).sort_index()
                      for name in (alias, canonical)]
            if not frames[0].index.equals(frames[1].index):
                raise ValueError("Original array-alias prediction identities differ")
            np.testing.assert_array_equal(frames[0][list(SOURCES)].values, frames[1][list(SOURCES)].values)
    write_json(output, dict(status="verified_complete_original_scope_reusable",
                            source_scope=str(source_scope), source_manifest_sha256=digest(manifest_path),
                            source_cohort_sha256=digest(original_path), cohort_sha256=digest(cohort_path)))
    print(f"Verified reusable {source_scope.name} genotypes, panels, inheritance and ten baseline models.", flush=True)


def fine_subset(reference, cap):
    grouped = defaultdict(list)
    for row in reference:
        grouped[row["fine"]].append(row)
    return [row for fine in sorted(grouped)
            for row in sorted(grouped[fine], key=lambda row: row["sample_id"])[:cap]]


def capacity(cohort_path, source_scope, smoke, output):
    cohort = read(cohort_path)
    audit = read(smoke / "summary/summary/independent_verification.json")
    manifest = read(smoke / "published_verification.json")
    post = read(smoke / "post_cleanup_verification.json")
    if (audit["status"] != "smoke_passed_not_a_whole_genome_leaderboard"
            or audit["model_configurations"] != 98
            or post["manifest_sha256"] != digest(smoke / "published_verification.json")):
        raise ValueError("Capacity promotion requires the independently verified98-model smoke")
    require_hashes(manifest["published_sha256"])
    require_hashes(manifest["source_sha256"])
    counts = {str(chromosome): read(source_scope / f"panels/chr{chromosome}/panel/panel.json")["marker_count"]
              for chromosome in range(1, 23)}
    spans = cohort["chromosome_spans"]
    ratio = max(sum(counts.values()) / counts["22"],
                sum(end - start for start, end in spans.values()) / (spans["22"][1] - spans["22"][0]))
    measured = sum(path.stat().st_size for name in ("panels", "models", "reference_bases")
                   for path in (smoke / name).rglob("*") if path.is_file())
    required, available = int(4 * measured * ratio + 20 * 1024**3), shutil.disk_usage(smoke).free
    write_json(output, dict(measured_smoke_model_panel_and_base_bytes=measured, scaling_ratio=ratio,
                           required_bytes=required, available_bytes=available,
                           includes_two_masks_and_work_publish_copies=True,
                           scientific_inputs_reused_without_duplicate_preprocessing=True))
    if available < required:
        raise ValueError(f"Full98-model capacity gate failed: need {required} bytes; available {available}")
    print(f"Full98-model capacity gate passed: need {required} bytes; available {available}.", flush=True)


def panels(cohort_path, source_panel, gate, output):
    cohort, certificate = read(cohort_path), read(gate)
    if (certificate["status"] != "verified_complete_original_scope_reusable"
            or certificate["cohort_sha256"] != digest(cohort_path)):
        raise ValueError("Reusable source gate belongs to a different expansion")
    metadata = read(source_panel / "panel.json")
    expected = Path(certificate["source_scope"]) / f"panels/chr{metadata['chromosome']}/panel"
    if source_panel.resolve() != expected.resolve():
        raise ValueError("Source panel is not from the verified production scope")
    shutil.copytree(source_panel, output)
    reference = cohort["reference"]
    ids = [row["sample_id"] for row in reference]
    if (output / "reference.ids").read_text().splitlines() != ids:
        raise ValueError("Original panel's fitting identities differ from the expanded cohort")
    original = np.load(output / "reference.npz")["haplotypes"].reshape(len(ids), 2, -1)
    keys = [tuple(key) for key in metadata["keys"]]
    if original.shape[2] != len(keys):
        raise ValueError("Original reference marker axes differ")
    counts = {"full": len(ids), "finecap": metadata["finecap_people"]}
    for cap in (8, 32):
        name = f"finecap{cap}"
        selected = fine_subset(reference, cap)
        if {row["source"] for row in selected} != set(SOURCES):
            raise ValueError("Reference cap loses native-source coverage")
        identities = [row["sample_id"] for row in selected]
        (output / f"{name}.ids").write_text("".join(identity + "\n" for identity in identities))
        (output / f"{name}.panel").write_text("".join(f"{row['sample_id']}\t{row['source']}\n" for row in selected))
        target = output / f"{name}.vcf.gz"
        run(["bcftools", "view", "-S", output / f"{name}.ids", "-Oz", "-o", target,
             output / "reference.vcf.gz"], output / f"{name}.subset")
        run(["bcftools", "index", target], output / f"{name}.index")
        actual = read_call_matrix(target, keys, identities, output, phased=True)
        np.testing.assert_array_equal(actual, original[[ids.index(identity) for identity in identities]])
        counts[name] = len(identities)
    metadata["expanded_reference_people"] = counts
    metadata["source_panel_sha256"] = digest(source_panel / "panel.json")
    metadata["source_gate_sha256"] = digest(gate)
    metadata["cohort_sha256"] = digest(cohort_path)
    write_json(output / "panel.json", metadata)
    shutil.copy2(gate, output / "source_gate.json")


def fit_base(folder, window, output):
    output.mkdir()
    metadata = read(folder / "panel.json")
    keys = [tuple(key) for key in metadata["keys"]]
    cls, windowed = native.load_upstream()
    layout = native.window_layout(keys, native.map_points(metadata["map"], metadata["chromosome"]),
                                  metadata["span"],
                                  dict(window_size_cM=window, architecture="fast", smooth_size=1))
    data = np.load(folder / "reference.npz")
    x_path, y_path = output / "reference_X.npy", output / "reference_y.npy"
    np.save(x_path, data["haplotypes"])
    np.save(y_path, np.repeat(np.repeat(data["labels"], 2)[:, None], len(keys), axis=1).astype(np.int8))
    model = cls(C=len(keys), M=layout["M"], A=8, S=1, mode="fast",
                snp_pos=layout["positions"], snp_ref=np.asarray([key[2] for key in keys]),
                snp_alt=np.asarray([key[3] for key in keys]), population_order=np.asarray(SOURCES),
                missing_encoding=2, n_jobs=1, calibrate=False, context_ratio=.5,
                mode_filter=False, seed=20261022, verbose=False)
    model.base.base_multithread = False
    for estimator in model.base.models:
        estimator.set_params(random_state=20261022)
    model.base.train_windowed(windowed([str(x_path)], [str(y_path)], len(keys), layout["M"]), verbose=False)
    with (output / "base.pkl").open("wb") as handle:
        pickle.dump(dict(model=model, layout=layout, window_size_cM=window,
                         panel_sha256=digest(folder / "panel.json")), handle)
    write_json(output / "base.json", dict(window_size_cM=window, base_sha256=digest(output / "base.pkl"),
                                         panel_sha256=digest(folder / "panel.json"),
                                         fitting_role="reference_only", source_order=SOURCES))
    x_path.unlink()
    y_path.unlink()


def fit_smoother(folder, metadata, policy, mask, base, output):
    native.load_upstream()
    record = read(base / "base.json")
    if (record["base_sha256"] != digest(base / "base.pkl")
            or record["panel_sha256"] != digest(folder / "panel.json")
            or record["window_size_cM"] != policy["window_size_cM"]):
        raise ValueError("Shared Gnomix reference base belongs to a different panel/window")
    with (base / "base.pkl").open("rb") as handle:
        payload = pickle.load(handle)
    model, layout = payload["model"], payload["layout"]
    model.smooth.model.CRF.set_params(c2=policy["crf_c2"])
    training = np.load(folder / f"{mask}.training.npz")
    labels = training["labels"]
    if labels.shape != training["haplotypes"].shape or not np.isin(labels, range(8)).all():
        raise ValueError("Invalid training-only native-eight-source labels")
    model.smooth.train(model.base.predict_proba(training["haplotypes"]), native.window_labels(labels, layout))
    if list(map(int, model.smooth.model.classes_)) != list(range(8)):
        raise ValueError("Expanded Gnomix changed its class axis")
    path = output / f"{mask}.model.pkl"
    with path.open("wb") as handle:
        pickle.dump(dict(model=model, layout=layout, source_order=SOURCES), handle)
    with path.open("rb") as handle:
        restored = pickle.load(handle)["model"]
    evaluation = np.load(folder / f"{mask}.evaluation.npz")["haplotypes"]
    before = model.smooth.predict_proba(model.base.predict_proba(evaluation))
    after = restored.smooth.predict_proba(restored.base.predict_proba(evaluation))
    np.testing.assert_allclose(before, after, rtol=0, atol=1e-10)
    if (after.shape != (2 * len(metadata["evaluation"]), layout["W"], 8)
            or not np.isfinite(after).all() or (after < 0).any()
            or not np.allclose(after.sum(axis=2), 1, rtol=0, atol=1e-10)):
        raise ValueError("Expanded Gnomix returned malformed local posteriors")
    np.savez_compressed(output / f"{mask}.local_posteriors.npz", probabilities=after,
                        weights=layout["physical_weights"])
    values = np.einsum("nwa,w->na", after.reshape(-1, 2, layout["W"], 8).mean(axis=1),
                       layout["physical_weights"]) / (metadata["span"][1] - metadata["span"][0])
    return np.stack([validate_vector(value) for value in values])


def flare_command(folder, metadata, policy, mask, prefix):
    name = "reference" if policy["reference"] == "full" else policy["reference"]
    return ["java", "-Xmx6g", "-jar", PROJECT / "results_v3_global_precision_01/precision.jar",
            f"gt={folder / f'{mask}.evaluation.vcf.gz'}", f"ref={folder / (name + '.vcf.gz')}",
            f"ref-panel={folder / (policy['reference'] + '.panel')}", f"map={metadata['map']}",
            f"out={prefix}", "nthreads=1", "seed=20261022", f"array={str(policy['array']).lower()}",
            f"em={str(policy['em']).lower()}", f"gen={policy['gen']}", "update-p=false",
            "probs=true", "min-maf=0", "min-mac=0"]


def adopt_baseline(cohort_path, folder, candidate, output):
    cohort, panel = read(cohort_path), read(folder / "panel.json")
    expansion = cohort["candidate_expansion"]
    original = read(Path(expansion["source_cohort"]))
    if original["protocol"]["candidates"][candidate] != cohort["protocol"]["candidates"][candidate]:
        raise ValueError("Baseline reuse would change the fitted candidate policy")
    source_scope = Path(read(folder / "source_gate.json")["source_scope"])
    source = source_scope / f"models/{candidate}/chr{panel['chromosome']}/model"
    complete = read(source / "complete.json")
    if (complete["cohort_sha256"] != expansion["source_cohort_sha256"]
            or complete["panel_sha256"] != panel["source_panel_sha256"]
            or complete["predictions_sha256"] != digest(source / "predictions.tsv")):
        raise ValueError("Baseline model provenance does not match the independently verified source")
    shutil.copytree(source, output)
    shutil.copy2(source / "complete.json", output / "original_complete.json")
    complete.update(cohort_sha256=digest(cohort_path), panel_sha256=digest(folder / "panel.json"),
                    execution_mode="reused_independently_verified_identical_baseline",
                    original_complete_sha256=digest(source / "complete.json"),
                    original_cohort_sha256=expansion["source_cohort_sha256"],
                    original_source_scope=str(source_scope))
    write_json(output / "complete.json", complete)


def evaluate(cohort_path, folder, candidate, base, output):
    cohort, metadata = read(cohort_path), read(folder / "panel.json")
    if candidate in cohort["candidate_expansion"]["baseline_candidates"]:
        adopt_baseline(cohort_path, folder, candidate, output)
        return
    output.mkdir()
    policy = cohort["protocol"]["candidates"][candidate]
    keys = [tuple(key) for key in metadata["keys"]]
    predictions = []
    for mask in metadata["masks"]:
        started = time.monotonic()
        if policy["family"] == "flare":
            prefix = output / mask
            run(flare_command(folder, metadata, policy, mask, prefix), output / f"{mask}.flare")
            from screen_model_families import read_model
            values = physical_flare(str(prefix) + ".anc.vcf.gz", read_model(str(prefix) + ".model")["ancestries"],
                                    keys, metadata["span"], metadata["evaluation"])
        elif policy["family"] == "gnomix" and base is not None:
            values = fit_smoother(folder, metadata, policy, mask, base, output)
        else:
            raise ValueError("New candidates must be supported native FLARE/Gnomix models")
        for identity, value in zip(metadata["evaluation"], values):
            value = validate_vector(value)
            predictions.append(dict(candidate=candidate, mask=mask, chromosome=metadata["chromosome"],
                                    sample_id=identity, span_bp=metadata["span"][1] - metadata["span"][0],
                                    seconds=time.monotonic() - started,
                                    **{source: float(value[i]) for i, source in enumerate(SOURCES)}))
        print(f"{candidate}/{mask}/chr{metadata['chromosome']}: verified eight-source predictions.", flush=True)
    pd.DataFrame(predictions).to_csv(output / "predictions.tsv", sep="\t", index=False)
    write_json(output / "complete.json", dict(
        candidate=candidate, chromosome=metadata["chromosome"], rows=len(predictions),
        predictions_sha256=digest(output / "predictions.tsv"), panel_sha256=digest(folder / "panel.json"),
        cohort_sha256=digest(cohort_path), source_order=SOURCES, production_changed=False,
        execution_mode="new_native_fit", candidate_policy=policy))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("freeze", "verify-source", "capacity", "panels", "base", "evaluate"))
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--source-experiment", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--source-scope", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--gate", type=Path)
    parser.add_argument("--window", type=float)
    parser.add_argument("--candidate")
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.stage == "freeze":
            freeze(args.source_experiment, args.plan, args.output)
        elif args.stage == "verify-source":
            verify_source(args.cohort, args.source_scope, args.output)
        elif args.stage == "capacity":
            capacity(args.cohort, args.source_scope, args.folder, args.output)
        elif args.stage == "panels":
            panels(args.cohort, args.folder, args.gate, args.output)
        elif args.stage == "base":
            fit_base(args.folder, args.window, args.output)
        else:
            evaluate(args.cohort, args.folder, args.candidate, args.base, args.output)
