"""Build five declared-disjoint roles from retained, verified observed genotypes."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from pathlib import Path
import shutil

import numpy as np

import core


def components(groups):
    parent = {person: person for group in groups for person in group}

    def root(person):
        while parent[person] != person:
            parent[person] = parent[parent[person]]
            person = parent[person]
        return person

    for group in groups:
        for person in group[1:]:
            parent[root(person)] = root(group[0])
    result = defaultdict(list)
    for person in sorted(parent):
        result[root(person)].append(person)
    return sorted(result.values(), key=lambda people: (-len(people), people))


def family_id(donors):
    return "founders_" + hashlib.sha256(",".join(sorted(donors)).encode()).hexdigest()[:16]


def write_rows(path, records, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def plan(source, template, output):
    project = source.parent
    original = core.read(source / "experiment.json")
    completed = core.read(source / "completion.json")
    if completed["status"] != "complete_native_observed_only_pilot" or completed["imputation_performed"]:
        raise ValueError("Require the verified completed observed-only input source")
    config = core.read(template)
    core.validate_config(config)
    cohort = core.read(source / "tasks/lonni_chr22/cohort.json")
    if cohort["sources"] != config["sources"] or cohort["chromosome_spans"] != config["spans"]:
        raise ValueError("Template source axis or physical spans changed")
    legacy_path = project / "results_v3_phasefree_array_validation_01/experiment.json"
    donors_path = project / "results_v3_phasefree_rich_01/donors.json"
    legacy = core.read(legacy_path)
    audit_path = project / "results_v3_model_family_screen_05/input_audit.json"
    truth_audit = core.read(audit_path)
    family_donors = core.read(donors_path)
    paths = [template, source / "completion.json", source / "selection.json",
             legacy_path, donors_path, audit_path]
    reference_path = source / "inputs/reference_labels.tsv"
    pure_path = source / "inputs/pure_train_labels.tsv"
    reference, pure = core.rows(reference_path), core.rows(pure_path)
    paths += [reference_path, pure_path]
    pedigree = []
    batches = []
    for batch in cohort["pedigree_train"]["batches"]:
        metadata = source / "inputs" / Path(batch["samples_tsv"]).name
        selected = core.rows(metadata)
        pedigree += selected
        batches.append({"batch": batch["batch"], "samples": [r["sample_id"] for r in selected]})
        paths.append(metadata)
    families = sorted({row["family_id"] for row in pedigree})
    grouped = components([family_donors[family] for family in families])
    if len(grouped) != 3 or sorted(map(len, grouped)) != [4, 5, 55]:
        raise ValueError("Retained founder components changed; explicitly redesign the role split")
    donor_role = {donor: role for group, role in zip(grouped, ("train", "validation", "test"))
                  for donor in group}
    donor_family = {donor: family_id(group) for group in grouped for donor in group}
    ref_ids = {row["sample_id"] for row in reference}
    pure_ids = {row["sample_id"] for row in pure}
    if ref_ids & pure_ids or set(legacy["donor_ids"]) & (ref_ids | pure_ids):
        raise ValueError("Retained founder/reference identities overlap")
    if not set(donor_role) <= pure_ids:
        raise ValueError("A pedigree donor is absent from the retained pure genotype panel")
    holdouts = {}
    for source_group in config["sources"]:
        candidates = sorted(
            (row["sample_id"] for row in reference if row["source_group"] == source_group),
            key=lambda identity: hashlib.sha256(f"{config['seed']}:{identity}".encode()).hexdigest())
        if len(candidates) < 6:
            raise ValueError("Insufficient reference people for source-stratified holdouts")
        for identity, role in zip(candidates[:5],
                                  ("calibration", "calibration", "calibration", "validation", "test")):
            holdouts[identity] = role
    samples = []

    def add(identity, role, donors, source_group, fine, track):
        samples.append({
            "sample_id": identity, "role": role,
            "family_id": donor_family.get(identity, family_id(donors)),
            "donors": ",".join(sorted(donors)), "source": source_group,
            "fine": fine, "track": track,
        })

    for row in reference:
        identity = row["sample_id"]
        add(identity, holdouts.get(identity, "reference"), [identity],
            row["source_group"], row["fine"], "pure")
    for row in pure:
        identity = row["sample_id"]
        add(identity, donor_role.get(identity, "train"), [identity],
            row["source_group"], row["fine"], "pure")
    for row in pedigree:
        donors = family_donors[row["family_id"]]
        roles = {donor_role[donor] for donor in donors}
        if len(roles) != 1:
            raise ValueError("Connected founders were split across roles")
        identity = row["sample_id"]
        donor_family[identity] = donor_family[donors[0]]
        add(identity, roles.pop(), donors, "", "", "mixed")
    legacy_groups = components([row["ancestors"] for row in legacy["samples"]])
    for group in legacy_groups:
        for donor in group:
            donor_family[donor] = family_id(group)
    for row in legacy["samples"]:
        identity = row["sample_id"]
        donors = row["ancestors"]
        donor_family[identity] = donor_family[donors[0]]
        label = ""
        if row["track"] == "pure":
            if row["truth"].count(1.0) != 1 or sum(row["truth"]) != 1:
                raise ValueError("Pure source truth is not one-hot")
            label = config["sources"][row["truth"].index(1.0)]
        add(identity, "calibration", donors, label, "", "pure" if label else "mixed")
    samples.sort(key=lambda row: (core.ROLES.index(row["role"]), row["sample_id"]))
    core.validate_roles(samples, config["sources"])
    output.mkdir()
    write_rows(output / "samples.tsv", samples,
               ["sample_id", "role", "family_id", "donors", "source", "track", "fine"])
    all_ids = {row["sample_id"] for row in samples}
    mixed_ids = {row["sample_id"] for row in samples if row["track"] == "mixed"}
    truths = []
    truth_paths = [project / "results_v3_model_family_screen_05/inputs" /
                   f"pedigree_train_{batch['batch']}.segments.tsv" for batch in batches]
    truth_paths += [project / "results_v3_model_family_screen_05/inputs" /
                    f"pedigree_{batch}.segments.tsv"
                    for batch in sorted({row["batch"] for row in legacy["samples"] if row["batch"] != "pure"})]
    for path in truth_paths:
        expected = truth_audit["input_sha256"]["/results/inputs/" + path.name]
        if core.digest(path) != expected:
            raise ValueError(f"Canonical truth changed: {path}")
        if path.name.startswith("pedigree_train_") and expected != original["input_sha256"]["/results/inputs/" + path.name]:
            raise ValueError("Training truth certificates disagree")
        paths.append(path)
        for row in core.rows(path):
            if row["sample_id"] in mixed_ids:
                truths.append({
                    "sample_id": row["sample_id"], "chromosome": row["chrom"],
                    "haplotype": str(int(row["haplotype"]) - 1),
                    "start": row["start_bp"], "end": row["end_bp"],
                    "source": "EAS_SEA" if row["source_group"] in ("EAS", "SEA") else row["source_group"],
                })
    for row in samples:
        if row["role"] == "reference" or row["track"] != "pure":
            continue
        for chrom in config["chromosomes"]:
            start, end = config["spans"][str(chrom)]
            for hap in (0, 1):
                truths.append({"sample_id": row["sample_id"], "chromosome": str(chrom),
                               "haplotype": str(hap), "start": str(start), "end": str(end),
                               "source": row["source"]})
    write_rows(output / "segments.tsv", truths,
               ["sample_id", "chromosome", "haplotype", "start", "end", "source"])
    truth_config = {**config, "segments": str(output / "segments.tsv")}
    integrated, total = None, 0
    for chrom in config["chromosomes"]:
        targets, _ = core.truth(truth_config, samples, chrom)
        length = config["spans"][str(chrom)][1] - config["spans"][str(chrom)][0]
        integrated = length * targets if integrated is None else integrated + length * targets
        total += length
    nonreference = [row["sample_id"] for row in samples if row["role"] != "reference"]
    for row in legacy["samples"]:
        if not np.allclose(integrated[nonreference.index(row["sample_id"])] / total,
                           row["truth"], rtol=0, atol=1e-10):
            raise ValueError("Physical segment truth disagrees with the retained realized ancestry")
    expected_legacy = {row["sample_id"] for row in legacy["samples"]}
    if len(all_ids) != len(samples) or not expected_legacy <= all_ids:
        raise ValueError("Incomplete or duplicated retained people")
    for path in paths:
        relative = None
        if path.is_relative_to(source):
            relative = "/benchmark/" + str(path.relative_to(source))
        elif path == legacy_path:
            relative = "/legacy/experiment.json"
        elif path == donors_path:
            relative = "/workspace/" + str(path.relative_to(project))
        if relative in original["input_sha256"] and core.digest(path) != original["input_sha256"][relative]:
            raise ValueError(f"Historical scientific input changed: {path}")
    core.write(output / "plan.json", {
        "source": str(source), "template": config, "batches": batches,
        "holdouts": sorted(holdouts), "legacy": legacy["samples"],
        "samples": samples, "input_sha256": {str(path): core.digest(path) for path in paths},
        "role_counts": dict(Counter(row["role"] for row in samples)),
        "role_tracks": {role: dict(Counter(row["track"] for row in samples if row["role"] == role))
                        for role in core.ROLES},
        "training_founder_components": grouped,
        "selection_uses_model_performance": False, "old_model_outputs_adopted": False,
        "independent_test": False, "genetic_relatedness_certified": False,
        "limitations": [
            "Previously studied development cohort, not a fresh independent test.",
            "Validation and test each contain one small founder-connected pedigree group.",
            "Reference-stratified pure holdouts are removed before any new reference fitting.",
            "Observed-only; no passing imputed genotype source is available.",
        ],
    })
    print("Declared-disjoint role counts:", Counter(row["role"] for row in samples), flush=True)


def chromosome(plan_dir, chrom, output):
    plan = core.read(plan_dir / "plan.json")
    source = Path(plan["source"])
    samples = plan["samples"]
    output.mkdir()
    work = output.parent / (output.name + ".assembly")
    work.mkdir()
    hashes = {}
    certificates = {mask: core.read(source / f"tasks/{mask}_chr{chrom}/inputs_complete.json")
                    for mask in ("lonni", "frederik")}

    def certified(mask, path):
        key = "/benchmark/" + str(path.relative_to(source))
        certificate = certificates[mask]
        if certificate["status"] != "complete" or core.digest(path) != certificate["output_sha256"][key]:
            raise ValueError(f"Observed genotype certificate failed: {path}")
        hashes[str(path)] = certificate["output_sha256"][key]
        return path

    def subset(path, identities, name):
        selected = work / f"{name}.samples"
        selected.write_text("".join(identity + "\n" for identity in identities))
        target = work / f"{name}.vcf.gz"
        core.run(["bcftools", "view", "-S", selected, "-Oz", "-o", target, path], work / name)
        core.run(["bcftools", "index", target], work / f"{name}.index")
        return target

    def merge(paths, target, name):
        if len(paths) == 1:
            shutil.copyfile(paths[0], target)
        else:
            core.run(["bcftools", "merge", "--merge", "none", "-Oz", "-o", target, *paths], work / name)
        core.run(["bcftools", "index", target], work / f"{name}.index")

    def union(filename, identities, name):
        paths = [certified(mask, source / f"tasks/{mask}_chr{chrom}/{filename}.vcf.gz")
                 for mask in ("lonni", "frederik")]
        left, a = core.vcf(paths[0], chrom, identities, True)
        right, b = core.vcf(paths[1], chrom, identities, True)
        li, ri = {key: i for i, key in enumerate(left)}, {key: i for i, key in enumerate(right)}
        common = sorted(set(left) & set(right))
        if not common or not np.array_equal(a[:, :, [li[key] for key in common]],
                                             b[:, :, [ri[key] for key in common]]):
            raise ValueError(f"Genotype/phase disagreement between retained masks: {name}")
        keys = sorted(set(left) | set(right), key=lambda key: key[1])
        if len({key[1] for key in keys}) != len(keys):
            raise ValueError("Ambiguous alleles at a union marker")
        del a, b
        inputs = [subset(path, identities, f"{name}_{i}") for i, path in enumerate(paths)]
        target = work / f"{name}.vcf.gz"
        core.run(["bcftools", "concat", "-a", "-d", "exact", "-Oz", "-o", target, *inputs], work / name)
        core.run(["bcftools", "index", target], work / f"{name}.index")
        return target, keys

    ref_people = [row["sample_id"] for row in samples if row["role"] == "reference"]
    reference, keys = union("reference", ref_people, "reference_union")
    shutil.copyfile(reference, output / "reference.vcf.gz")
    shutil.copyfile(str(reference) + ".csi", output / "reference.vcf.gz.csi")
    train_people = {row["sample_id"] for row in samples if row["role"] == "train"}
    pure_ids = [row["sample_id"] for row in core.rows(source / "inputs/pure_train_labels.tsv")]
    groups = [("pure_train", pure_ids)] + [(batch["batch"], batch["samples"]) for batch in plan["batches"]]
    training = []
    for filename, identities in groups:
        chosen = sorted(set(identities) & train_people)
        if not chosen:
            continue
        path, panel = union(filename, chosen, filename + "_union")
        if panel != keys:
            raise ValueError("Training/reference exact union panels differ")
        training.append(path)
    merge(training, output / "training.vcf.gz", "training")
    query_people = sorted(row["sample_id"] for row in samples
                          if row["role"] in ("calibration", "validation", "test"))
    for mask in ("lonni", "frederik"):
        task = source / f"tasks/{mask}_chr{chrom}"
        parts = [subset(certified(mask, task / "reference.vcf.gz"), plan["holdouts"], mask + "_holdouts")]
        for filename, identities in groups:
            chosen = sorted(set(identities) & set(query_people))
            if chosen:
                parts.append(subset(certified(mask, task / f"{filename}.vcf.gz"), chosen,
                                    mask + "_" + filename))
        for row in plan["legacy"]:
            identity = row["sample_id"]
            parts.append(certified(mask, task / f"cases/{identity}/inherited.vcf.gz"))
        combined = work / f"query_{mask}.vcf.gz"
        merge(parts, combined, "merge_" + mask)
        querykeys, calls = core.vcf(combined, chrom, query_people, True)
        if not set(querykeys) <= set(keys):
            raise ValueError("Query contains keys outside the shared union reference")
        target = output / f"query_{mask}.vcf.gz"
        core.write_vcf(target, querykeys, query_people, calls, phased=False)
        observed_keys, observed = core.vcf(target, chrom, query_people, False)
        if observed_keys != querykeys or not np.array_equal(observed.sum(axis=1), calls.sum(axis=1)):
            raise ValueError("Harmonized adapter changed diploid genotypes")
    outputs = [path for path in output.iterdir() if path.is_file()]
    core.write(output / "audit.json", {
        "chromosome": chrom, "union_markers": len(keys), "old_model_outputs_adopted": False,
        "input_sha256": hashes, "output_sha256": {path.name: core.digest(path) for path in outputs},
        "mask_overlap_gt_phase_equal": True, "query_gt_unchanged": True,
    })
    print(f"chr{chrom}: union panel {len(keys)}; both harmonized masks verified", flush=True)


def finalize(plan_dir, bundles, destination, output):
    plan = core.read(plan_dir / "plan.json")
    config = plan["template"]
    expected = config["chromosomes"]
    audits = [core.read(path / "audit.json") for path in bundles]
    if sorted(record["chromosome"] for record in audits) != expected:
        raise ValueError("Incomplete chromosome bundle")
    for path, record in zip(bundles, audits):
        for name, value in record["output_sha256"].items():
            if core.digest(path / name) != value:
                raise ValueError(f"Changed prepared input: {path / name}")
    config["samples"] = str(destination / "metadata/plan/samples.tsv")
    config["segments"] = str(destination / "metadata/plan/segments.tsv")
    for kind in ("reference", "training"):
        config[kind] = {str(chrom): str(destination / f"chromosomes/chr{chrom}/bundle/{kind}.vcf.gz")
                        for chrom in expected}
    config["routes"] = [
        {"id": f"harmonized_{mask}", "stage": "harmonized",
         "vcfs": {str(chrom): str(destination / f"chromosomes/chr{chrom}/bundle/query_{mask}.vcf.gz")
                  for chrom in expected}} for mask in ("lonni", "frederik")]
    config["independent_test"] = False
    core.validate_config(config)
    output.mkdir()
    core.write(output / "full.json", config)
    smoke = {**config, "mode": "smoke", "chromosomes": [22]}
    chosen = {
        min((row["sample_id"] for row in plan["samples"] if row["role"] == role and row["track"] == track))
        for role in ("calibration", "validation", "test") for track in ("pure", "mixed")
    }
    smoke_samples = [row for row in plan["samples"]
                     if row["role"] in ("reference", "train") or row["sample_id"] in chosen]
    core.validate_roles(smoke_samples, config["sources"])
    write_rows(output / "smoke_samples.tsv", smoke_samples,
               ["sample_id", "role", "family_id", "donors", "source", "track", "fine"])
    smoke["samples"] = str(destination / "ready/smoke_samples.tsv")
    for kind in ("reference", "training", "maps", "spans"):
        smoke[kind] = {"22": config[kind]["22"]}
    smoke["routes"] = [{**route, "vcfs": {"22": route["vcfs"]["22"]}} for route in config["routes"]]
    core.validate_config(smoke)
    core.write(output / "smoke.json", smoke)
    core.write(output / "verification.json", {
        "status": "complete_verified_observed_inputs", "role_counts": plan["role_counts"],
        "role_tracks": plan["role_tracks"], "chromosomes": expected,
        "routes": [route["id"] for route in config["routes"]],
        "models": [row["id"] for row in config["models"]],
        "interpretations": [row["id"] for row in config["interpretations"]],
        "model_handover": False, "imputation_performed": False,
        "independent_test": False, "limitations": plan["limitations"], "chromosome_audits": audits,
        "smoke_query_samples": sorted(chosen), "smoke_is_mechanics_only": True,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("plan", "chromosome", "finalize"))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--bundles", nargs="+", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.step == "plan":
        plan(args.source, args.template, args.output)
    elif args.step == "chromosome":
        chromosome(args.plan, args.chromosome, args.output)
    else:
        finalize(args.plan, args.bundles, args.destination, args.output)


if __name__ == "__main__":
    main()
