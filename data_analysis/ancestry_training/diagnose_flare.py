"""Matched, non-promotional FLARE mechanism experiments on preserved real inputs."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from pathlib import Path
import shutil
import sys
import time
import zipfile

import numpy as np

import core
import models


PROJECT = next((parent / "FLARE/SNPster_v3_100model_pedigree"
                for parent in Path(__file__).resolve().parents
                if (parent / "FLARE/SNPster_v3_100model_pedigree").is_dir()), None)
if PROJECT is None:
    raise RuntimeError("Cannot locate the historical model helpers beside the ancestry-training workspace")
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from screen_model_families import read_model, write_model


def kinship_coefficients(reference, target):
    reference, target = np.asarray(reference), np.asarray(target)
    if (reference.ndim != 2 or target.shape != reference.shape[1:]
            or not np.isin(reference, [0, 1, 2]).all() or not np.isin(target, [0, 1, 2]).all()):
        raise ValueError("Kinship screening requires complete diploid SNP dosages on one exact axis")
    heterozygotes = (reference == 1).sum(axis=1) + (target == 1).sum()
    if (heterozygotes == 0).any():
        raise ValueError("Kinship coefficient is undefined without heterozygous markers")
    both = ((reference == 1) & (target == 1)).sum(axis=1)
    ibs0 = (np.abs(reference.astype(int) - target.astype(int)) == 2).sum(axis=1)
    return (both - 2 * ibs0) / heterozygotes, ibs0


def crossed_parameters(destination, donor, components):
    result = dict(destination)
    rows = [donor["ancestries"].index(name) for name in destination["ancestries"]]
    columns = [donor["panels"].index(name) for name in destination["panels"]]
    for component in components:
        if component in ("p", "theta"):
            result[component] = donor[component][np.ix_(rows, columns)]
        elif component == "rho":
            result[component] = donor[component][rows]
        else:
            raise ValueError(f"Unsupported crossed model component: {component}")
    return result


def native_vector(path, sources):
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1 or not set(sources) <= set(rows[0]):
        raise ValueError("Native global ancestry axes differ")
    return core.probability([float(rows[0][source]) for source in sources], len(sources), tolerance=.005)


def forensics(parent, protocol_path, output):
    protocol = core.read(protocol_path)
    if protocol.get("round", 1) not in (1, 2, 3, 4):
        raise ValueError("Unsupported diagnostic round")
    root = parent / "full/results"
    config = core.read(root / "frozen/frozen/config.json")
    manifest = core.read(root / "frozen/frozen/manifest.json")
    if config["mode"] != "full" or config["independent_test"]:
        raise ValueError("Require the preserved development run, never a fresh final test")
    output.mkdir()
    samples = core.rows(config["samples"])
    source = PROJECT / "results_v3_observed_dense_01"
    sources = config["sources"]
    cases = list(protocol["cases"])
    historical_reference_ids = {row["sample_id"] for row in core.rows(source / "inputs/reference_labels.tsv")}
    for source_group in protocol.get("control_sources", ["SSA", "EAS_SEA"]):
        choices = [row for row in samples if row["role"] == "train"
                   and row["track"] == "pure" and row["source"] == source_group]
        chosen = min(choices, key=lambda row: hashlib.sha256(
            f"{protocol['seed']}:{row['sample_id']}".encode()).hexdigest())
        cases.append(chosen["sample_id"])
    if len(set(cases)) != len(cases):
        raise ValueError("Duplicate diagnostic cases")
    evidence = core.read(root / "summary/summary/prediction_evidence.json")
    old = {(row["route"], row["model"], row["sample_id"]): row for row in evidence
           if row["interpreter"] == "native" and row["role"] == "validation"}
    hashes = {}

    def checked(path):
        path = Path(path)
        value = core.digest(path)
        if str(path) in manifest["input_sha256"] and manifest["input_sha256"][str(path)] != value:
            raise ValueError(f"Frozen scientific input changed: {path}")
        hashes[str(path)] = value
        return path

    for path in (config["samples"], config["segments"], protocol_path):
        checked(path)
    authority_path = checked("/srv/dependencies/ancestry/flare2_test/source_model/ancestry_role_audit.tsv")
    authority = {}
    for row in core.rows(authority_path):
        if row["role"] == "source":
            broad = "EAS_SEA" if row["source_group"] in ("EAS", "SEA") else row["source_group"]
            key = row["sample"], row["fine"]
            if key in authority and authority[key] != broad:
                raise ValueError("Contradictory authoritative reference labels")
            authority[key] = broad
    for row in samples:
        if row["sample_id"] in cases and authority.get((row["sample_id"], row["fine"])) != row["source"]:
            raise ValueError(f"Pure sample label differs from source authority: {row['sample_id']}")
    jar = Path(config["tools"]["flare_jar"])
    checked(jar)
    stock = Path("/opt/flare/flare.jar")
    shutil.copyfile(stock, output / "stock.jar")
    with zipfile.ZipFile(stock) as left, zipfile.ZipFile(jar) as right:
        changed = [name for name in left.namelist() if name.endswith(".class")
                   and left.read(name) != right.read(name)]
    if changed != ["admix/EstimatedGlobalAncProportions.class"]:
        raise ValueError(f"Stock/precision differences extend beyond export formatting: {changed}")
    datasets, audit = {}, []
    chrom = protocol["chromosome"]
    for mask in protocol["masks"]:
        route = f"harmonized_{mask}"
        folder = root / f"panels/{route}/chr{chrom}/panel"
        meta = core.read(checked(folder / "panel.json"))
        for row in meta["reference"]:
            if authority.get((row["sample_id"], row["fine"])) != row["source"]:
                raise ValueError(f"Reference label differs from source authority: {row['sample_id']}")
        data = np.load(checked(folder / "data.npz"))
        if core.digest(folder / "data.npz") != meta["data_sha256"]:
            raise ValueError("Published genotype tensor changed")
        keys = [tuple(key) for key in meta["keys"]]
        phased_path = checked(root / f"phased/{route}/chr{chrom}/phased/query.npy")
        phase_record = core.read(checked(phased_path.with_name("phase.json")))
        if core.digest(phased_path) != phase_record["query_sha256"]:
            raise ValueError("Published phase tensor changed")
        phased = np.load(phased_path)
        np.testing.assert_array_equal(phased.sum(axis=1), data["query"].sum(axis=1))
        task = source / f"tasks/{mask}_chr{chrom}"
        certificate = core.read(checked(task / "inputs_complete.json"))
        query_ids = [row["sample_id"] for row in meta["queries"]]
        for person in cases:
            row = next(row for row in samples if row["sample_id"] == person)
            if row["role"] not in ("train", "validation") or row["track"] != "pure":
                raise ValueError("Only declared training/validation pure controls may enter this diagnosis")
            case = output / f"{mask}_{person}"
            case.mkdir()
            original = task / ("reference.vcf.gz" if person in historical_reference_ids
                               else "pure_train.vcf.gz")
            original = checked(original)
            certificate_key = "/benchmark/" + str(original.relative_to(source))
            if core.digest(original) != certificate["output_sha256"][certificate_key]:
                raise ValueError("Inherited-phase input no longer matches its certificate")
            original_keys, inherited = core.vcf(original, chrom, [person], True)
            offsets = {key: index for index, key in enumerate(original_keys)}
            inherited = inherited[:, :, [offsets[key] for key in keys]]
            if row["role"] == "validation":
                current = phased[query_ids.index(person):query_ids.index(person) + 1]
                np.testing.assert_array_equal(current.sum(axis=1), inherited.sum(axis=1))
            else:
                current = inherited
            core.write_vcf(case / "target.vcf.gz", keys, [person], current)
            core.write_vcf(case / "inherited.vcf.gz", keys, [person], inherited)
            truth = np.eye(len(sources))[sources.index(row["source"])]
            reference_info = {}
            for reference, model_id in (("full", "flare_full_fixed_g010"),
                                        ("cap16", "flare_cap16_fixed_g010")):
                fitted = core.read(checked(root / f"models/{route}/chr{chrom}/{model_id}/fitted/model.json"))
                base = root / f"bases/{route}/chr{chrom}/{fitted['base_id']}/base"
                record = core.read(checked(base / "base.json"))
                for name, value in record["artifact_sha256"].items():
                    if core.digest(checked(base / name)) != value:
                        raise ValueError("Published reference base changed")
                selected = models.reference_subset(meta, {"reference_cap": 16} if reference == "cap16" else {})
                identities = [meta["reference"][index]["sample_id"] for index in selected]
                actual_keys, actual = core.vcf(base / "reference.vcf.gz", chrom, identities, True)
                if actual_keys != keys:
                    raise ValueError("Native reference exact allele axis differs")
                np.testing.assert_array_equal(actual, data["reference"][selected])
                labels = dict(line.split() for line in (base / "reference.panel").read_text().splitlines())
                expected_labels = {meta["reference"][index]["sample_id"]: meta["reference"][index]["source"]
                                   for index in selected}
                if labels != expected_labels or set(labels) & set(query_ids):
                    raise ValueError("Native reference labels or exclusions differ")
                ref = case / reference
                ref.mkdir()
                for name in ("reference.vcf.gz", "reference.vcf.gz.csi", "reference.panel"):
                    shutil.copyfile(base / name, ref / name)
                reversed_ids = identities[::-1]
                ids = ref / "reverse.samples"
                ids.write_text("".join(identity + "\n" for identity in reversed_ids))
                core.run(["bcftools", "view", "-S", ids, "-Oz", "-o", ref / "reverse.vcf.gz",
                          ref / "reference.vcf.gz"], ref / "reverse")
                core.run(["bcftools", "index", ref / "reverse.vcf.gz"], ref / "reverse.index")
                reverse_keys, reverse = core.vcf(ref / "reverse.vcf.gz", chrom, identities, True)
                if reverse_keys != keys:
                    raise ValueError("Reference permutation changed alleles")
                np.testing.assert_array_equal(reverse, actual)
                nearest = []
                for source_group in sources:
                    indices = [index for index, identity in enumerate(identities) if labels[identity] == source_group]
                    distances = np.mean(np.abs(actual[indices].sum(axis=1) - current.sum(axis=1)), axis=1)
                    nearest.append({"source": source_group, "minimum_dosage_L1_per_snp": float(distances.min())})
                phi, ibs0 = kinship_coefficients(actual.sum(axis=1), current[0].sum(axis=0))
                threshold = protocol.get("kinship_threshold", .177)
                related = [identity for identity, value in zip(identities, phi) if value > threshold]
                sham = []
                if protocol.get("round", 1) in (2, 3, 4):
                    keep = [identity for identity in identities if identity not in related]
                    if {labels[identity] for identity in keep} != set(sources):
                        raise ValueError("Genotype-only relatedness exclusion loses source coverage")
                    selections = {"unrelated": keep}
                    selected_references = protocol.get("references_by_case", {}).get(person, ["cap16", "full"])
                    if protocol.get("round", 1) == 4 and reference in selected_references:
                        fine = {row["sample_id"]: row["fine"] for row in meta["reference"]}
                        for population, count in Counter(fine[identity] for identity in related).items():
                            eligible = [identity for identity, value in zip(identities, phi)
                                        if fine[identity] == population and value < .0884 and identity not in related]
                            eligible.sort(key=lambda identity: hashlib.sha256(
                                f"{protocol['seed']}:{identity}".encode()).hexdigest())
                            if len(eligible) < count:
                                raise ValueError("Insufficient unrelated same-fine-population sham controls")
                            sham += eligible[:count]
                        if not related or len(sham) != len(related) or set(sham) & set(related):
                            raise ValueError("Sham exclusion must match real exclusions without removing relatives")
                        selections["sham"] = [identity for identity in identities if identity not in sham]
                    if person == "HG02954" and reference == "cap16":
                        selections.update({f"without_{identity}": [name for name in identities if name != identity]
                                           for identity in related})
                    for tag, retained in selections.items():
                        sample_file = ref / f"{tag}.samples"
                        sample_file.write_text("".join(identity + "\n" for identity in retained))
                        core.run(["bcftools", "view", "-S", sample_file, "-Oz", "-o", ref / f"{tag}.vcf.gz",
                                  ref / "reference.vcf.gz"], ref / tag)
                        core.run(["bcftools", "index", ref / f"{tag}.vcf.gz"], ref / f"{tag}.index")
                        selected_keys, selected_gt = core.vcf(ref / f"{tag}.vcf.gz", chrom, retained, True)
                        if selected_keys != keys:
                            raise ValueError("Kinship exclusion changed the marker allele axis")
                        np.testing.assert_array_equal(selected_gt, actual[[identities.index(identity) for identity in retained]])
                reference_info[reference] = {
                    "directory": str(ref.relative_to(output)), "reference_people": len(identities),
                    "counts": dict(Counter(labels.values())), "nearest_by_source": nearest,
                    "kinship_threshold": threshold, "related_ids": related,
                    "related_pairs": [{"reference": identity, "coefficient": float(value), "ibs0": int(zero)}
                                      for identity, value, zero in zip(identities, phi, ibs0) if value > threshold],
                    "sham_excluded_ids": sham,
                    "native_model": str(checked(root / f"predictions/{route}/chr{chrom}/{model_id}/prediction/HG02954.model")),
                }
            dataset_id = f"{mask}_{person}"
            datasets[dataset_id] = {
                "sample_id": person, "source": row["source"], "role": row["role"],
                "control_phase": "observed_phase_only" if row["role"] == "validation" else "inherited_training_phase",
                "mask": mask, "chromosome": chrom, "sources": sources,
                "keys": keys, "span": meta["span"], "references": reference_info,
                "map": str(checked(config["maps"][str(chrom)])), "jar": str(jar),
                "targets": {"observed": str((case / "target.vcf.gz").relative_to(output)),
                            "inherited": str((case / "inherited.vcf.gz").relative_to(output))},
                "truth": truth.tolist(),
                "expected_baselines": {
                    reference + "_" + mode: old.get((route, f"flare_{'full' if reference == 'full' else 'cap16'}_{mode}_g010", person))
                    for reference in ("full", "cap16") for mode in ("fixed", "em")},
            }
            audit.append({"dataset": dataset_id, "person": person, "role": row["role"],
                          "source": row["source"], "markers": len(keys),
                          "reference_labels_and_gt_match": True, "diploid_gt_unchanged": True})
    jobs = []

    def add(dataset_id, reference, name, **changes):
        jobs.append({"id": f"{dataset_id}_{reference}_{name}", "dataset": dataset_id,
                     "reference": reference, "contrast": name, **changes})

    for identity, dataset in datasets.items():
        references = protocol.get("references_by_case", {}).get(dataset["sample_id"], ["cap16", "full"])
        for reference in references:
            for em in (False, True):
                add(identity, reference, "baseline_" + ("em" if em else "fixed"), em=em)
        if protocol.get("round", 1) == 4:
            for reference in references:
                for filtered in ("unrelated", "sham"):
                    add(identity, reference, f"{filtered}_fixed", filtered_reference=filtered)
                    add(identity, reference, f"{filtered}_em", filtered_reference=filtered, em=True)
                    add(identity, reference, f"{filtered}_locked_fixed",
                        filtered_reference=filtered, parameter_model={})
        elif protocol.get("round", 1) == 3:
            for reference in references:
                add(identity, reference, "unrelated_fixed", filtered_reference="unrelated")
                add(identity, reference, "unrelated_em", filtered_reference="unrelated", em=True)
                add(identity, reference, "unrelated_locked_fixed", filtered_reference="unrelated", parameter_model={})
        elif protocol.get("round", 1) == 2:
            for reference in references:
                add(identity, reference, "explicit_fixed", parameter_model={})
                for name, components in [("rho", ["rho"]), ("p", ["p"]), ("all", ["p", "rho", "theta"])]:
                    add(identity, reference, f"crossed_{name}_fixed",
                        parameter_model={"crossed_components": components})
                add(identity, reference, "rho100_fixed", parameter_model={"equal_rho": 100.})
                add(identity, reference, "theta001_fixed", parameter_model={"equal_theta": .001})
                add(identity, reference, "unrelated_fixed", filtered_reference="unrelated")
                add(identity, reference, "unrelated_em", filtered_reference="unrelated", em=True)
                add(identity, reference, "unrelated_locked_fixed", filtered_reference="unrelated", parameter_model={})
                if dataset["sample_id"] == "HG02954" and reference == "cap16":
                    for relative in dataset["references"][reference]["related_ids"]:
                        add(identity, reference, f"without_{relative}_fixed",
                            filtered_reference=f"without_{relative}", parameter_model={})
        elif dataset["sample_id"] == "HG02954":
            for name, changes in [
                ("stock_fixed", {"stock": True}), ("stock_em", {"stock": True, "em": True}),
                ("inherited_fixed", {"target": "inherited"}),
                ("inherited_em", {"target": "inherited", "em": True}),
                ("reverse_fixed", {"reverse": True}), ("reverse_em", {"reverse": True, "em": True}),
                ("states500_fixed", {"params": {"states": 500}}),
                ("states500_em", {"params": {"states": 500}, "em": True}),
                ("ibs16_fixed", {"params": {"ibs-haps": 16}}),
                ("ibs16_em", {"params": {"ibs-haps": 16}, "em": True}),
                ("ibsstep05_fixed", {"params": {"ibs-step": .05}}),
                ("ibsstep10_fixed", {"params": {"ibs-step": .1}}),
                ("threads2_fixed", {"params": {"nthreads": 2}}),
                ("seed2_fixed", {"params": {"seed": protocol["seed"] + 1}}),
                ("seed2_em", {"params": {"seed": protocol["seed"] + 1}, "em": True}),
                ("updatep_em", {"params": {"update-p": True}, "em": True}),
            ]:
                add(identity, "cap16", name, **changes)
    core.write(output / "datasets.json", datasets)
    core.write(output / "jobs.json", jobs)
    core.write(output / "audit.json", {
        "parent": str(parent), "protocol": protocol, "input_sha256": hashes,
        "stock_changed_classes": changed, "checks": audit, "jobs": len(jobs),
        "test_people_used": False, "production_promotion": False,
        "reused_prior_evidence": "Native theta-crossing and diagonal-copying controls already tested historically; not presumed solutions.",
    })
    print(f"Prepared {len(jobs)} matched tests; source/genotype/phase/reference axes verified", flush=True)


def evaluate(prepared, identity, output):
    jobs = core.read(prepared / "jobs.json")
    job = next(row for row in jobs if row["id"] == identity)
    dataset = core.read(prepared / "datasets.json")[job["dataset"]]
    protocol = core.read(prepared / "audit.json")["protocol"]
    output.mkdir()
    ref = prepared / dataset["references"][job["reference"]]["directory"]
    target = prepared / dataset["targets"][job.get("target", "observed")]
    reference = ref / ("reverse.vcf.gz" if job.get("reverse") else "reference.vcf.gz")
    if job.get("filtered_reference"):
        reference = ref / f"{job['filtered_reference']}.vcf.gz"
    params = {"gen": 10, "em": job.get("em", False), "update-p": False, "states": 100,
              "panel-weight": .8, "nthreads": 1, "seed": protocol["seed"],
              "min-maf": 0, "min-mac": 0, "probs": True, "debug": True, **job.get("params", {})}
    command = ["java", core.java_heap(), "-jar", prepared / "stock.jar" if job.get("stock") else dataset["jar"],
               f"gt={target}", f"ref={reference}", f"ref-panel={ref / 'reference.panel'}",
               f"map={dataset['map']}", f"out={output / 'native'}"]
    command += [f"{key}={str(value).lower()}" for key, value in params.items()]
    if "parameter_model" in job:
        intervention = job["parameter_model"]
        native = read_model(dataset["references"][job["reference"]]["native_model"])
        if intervention.get("crossed_components"):
            donor_reference = "full" if job["reference"] == "cap16" else "cap16"
            donor = read_model(dataset["references"][donor_reference]["native_model"])
            native = crossed_parameters(native, donor, intervention["crossed_components"])
        if "equal_rho" in intervention:
            native["rho"] = np.full_like(native["rho"], intervention["equal_rho"])
        if "equal_theta" in intervention:
            native["theta"] = np.full_like(native["theta"], intervention["equal_theta"])
        write_model(output / "input.model", native)
        command += [f"model={output / 'input.model'}"]
    started = time.monotonic()
    core.run(command, output / "run")
    native_model = read_model(output / "native.model")
    values, correction = models.flare_posterior(
        output / "native.anc.vcf.gz", dataset["keys"], dataset["sources"], dataset["span"], dataset["sample_id"])
    values = values[[native_model["ancestries"].index(source) for source in dataset["sources"]]]
    global_values = native_vector(output / "native.global.anc.gz", dataset["sources"])
    expected = dataset["expected_baselines"].get(job["reference"] + ("_em" if params["em"] else "_fixed"))
    expected_error = None
    if job["contrast"].startswith("baseline_") and expected is not None:
        parent = Path(core.read(prepared / "audit.json")["parent"])
        route = f"harmonized_{dataset['mask']}"
        model_id = f"flare_{job['reference']}_{'em' if params['em'] else 'fixed'}_g010"
        original = parent / f"full/results/predictions/{route}/chr{dataset['chromosome']}/{model_id}/prediction/predictions.npz"
        meta = core.read(original.with_name("prediction.json"))
        index = [row["sample_id"] for row in meta["samples"]].index(dataset["sample_id"])
        expected_error = float(np.max(np.abs(values - np.load(original)["values"][index])))
        if expected_error > 1e-10:
            raise ValueError(f"Exact replay differs from the preserved baseline by {expected_error}")
    core.write(output / "result.json", {
        **job, "sample_id": dataset["sample_id"], "role": dataset["role"],
        "source": dataset["source"], "mask": dataset["mask"], "sources": dataset["sources"],
        "predictions": values.tolist(), "native_global": global_values.tolist(),
        "truth": dataset["truth"], "tv_pp": float(50 * np.abs(values - dataset["truth"]).sum()),
        "generation": native_model["generation"],
        "related_pairs": dataset["references"][job["reference"]]["related_pairs"],
        "sham_excluded_ids": dataset["references"][job["reference"]]["sham_excluded_ids"],
        "mu": dict(zip(native_model["ancestries"], native_model["mu"].tolist())),
        "p": native_model["p"].tolist(), "theta": native_model["theta"].tolist(),
        "rho": native_model["rho"].tolist(), "native_names": native_model["ancestries"],
        "seconds": time.monotonic() - started, "max_rounding_correction": correction,
        "baseline_max_fraction_delta": expected_error,
        "input_sha256": {str(path): core.digest(path) for path in (target, reference, ref / "reference.panel")},
        "output_sha256": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
    })
    print(f"{identity}: {dataset['source']} TV={50 * np.abs(values - dataset['truth']).sum():.4f} pp", flush=True)


def summarize(prepared, results, output):
    datasets = core.read(prepared / "datasets.json")
    expected = {row["id"] for row in core.read(prepared / "jobs.json")}
    records = [core.read(path / "result.json") for path in results]
    if {row["id"] for row in records} != expected or len(records) != len(expected):
        raise ValueError("Incomplete or duplicate diagnostic coverage")
    for path, record in zip(results, records):
        for name, value in record["output_sha256"].items():
            if core.digest(path / name) != value:
                raise ValueError("Diagnostic native output changed")
        measured = float(50 * np.abs(np.asarray(record["predictions"]) - record["truth"]).sum())
        if not np.isclose(measured, record["tv_pp"], atol=1e-10, rtol=0):
            raise ValueError("Diagnostic score does not reconstruct")
        dataset = datasets[record["dataset"]]
        keys, span = dataset["keys"], dataset["span"]
        boundaries = [span[0]] + [(left[1] + right[1]) / 2 for left, right in zip(keys, keys[1:])] + [span[1]]
        reconstructed = [0.] * len(dataset["sources"])
        count = 0
        with gzip.open(path / "native.anc.vcf.gz", "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if count >= len(keys) or [fields[0], int(fields[1]), fields[3], fields[4]] != keys[count]:
                    raise ValueError("Independent posterior audit: exact allele axis changed")
                formats, calls = fields[8].split(":"), fields[9].split(":")
                for haplotype in ("ANP1", "ANP2"):
                    values = [float(value) for value in calls[formats.index(haplotype)].split(",")]
                    if (len(values) != len(reconstructed) or not all(np.isfinite(values))
                            or any(value < 0 or value > 1 for value in values) or sum(values) <= 0
                            or sum(max(value - .005, 0) for value in values) > 1 + 1e-12
                            or sum(min(value + .005, 1) for value in values) < 1 - 1e-12):
                        raise ValueError("Independent posterior audit: invalid rounded probability mass")
                    factor = (boundaries[count + 1] - boundaries[count]) / (2 * (span[1] - span[0]) * sum(values))
                    for index, source in enumerate(dataset["sources"]):
                        reconstructed[index] += factor * values[record["native_names"].index(source)]
                count += 1
        if count != len(keys):
            raise ValueError("Independent posterior audit: incomplete marker coverage")
        np.testing.assert_allclose(reconstructed, record["predictions"], atol=1e-10, rtol=0)
    output.mkdir()
    core.write(output / "results.json", records)
    columns = ["id", "sample_id", "role", "source", "mask", "reference", "contrast", "tv_pp", "generation", "seconds"]
    with (output / "comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda row: row["id"]))
    core.write(output / "verification.json", {
        "verified": True, "calls": len(records), "test_people_used": False,
        "independent_native_posterior_reconstruction": True,
        "production_promotion": False,
        "files": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
    })
    print(f"Verified all {len(records)} diagnostic calls; no final-test model ranking", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("prepare", "evaluate", "summarize"))
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--identity")
    parser.add_argument("--results", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.step == "prepare":
        forensics(args.parent, args.protocol, args.output)
    elif args.step == "evaluate":
        evaluate(args.prepared, args.identity, args.output)
    else:
        summarize(args.prepared, args.results, args.output)
