"""Verified query-subset reuse without rebinding or modifying fitted models."""
import ast
import copy
import csv
import hashlib
import logging
from pathlib import Path
import shlex
import shutil

import numpy as np

import core
import models
from reuse import NATIVE_SUFFIXES, REQUIRED_NATIVE, inventory


def signature(panel, params):
    meta = core.read(panel / "panel.json")
    if core.digest(panel / "data.npz") != meta["data_sha256"]:
        raise ValueError("Projection reference panel changed")
    with np.load(panel / "data.npz") as data:
        reference = np.ascontiguousarray(data["reference"][models.reference_subset(meta, params)])
    selected = models.reference_subset(meta, params)
    chrom = str(meta["chromosome"])
    return {
        "reference_sha256": hashlib.sha256(reference.tobytes()).hexdigest(),
        "reference_shape": list(reference.shape), "reference_dtype": reference.dtype.str,
        "reference_people": [meta["reference"][index] for index in selected],
        "sources": meta["sources"], "keys": meta["keys"], "span": meta["span"],
        "route": {"id": meta["route"]["id"], "stage": meta["route"]["stage"],
                  "vcf": meta["route"]["vcfs"][chrom]},
        "params": params, "seed": meta["config"]["seed"],
        "map_sha256": core.digest(meta["config"]["maps"][chrom]),
        "tool_sha256": meta["tool_sha256"],
        "reference_filters": meta["config"]["reference_filters"],
        "target_policy": "one independently phased singleton per FLARE call; nthreads=1",
    }


def checked_files(record):
    folder = Path(record["directory"])
    for name, expected in record["files"].items():
        if Path(name).name != name or (folder / name).is_symlink() or core.digest(folder / name) != expected:
            raise ValueError(f"Projected source artifact changed or is work-dependent: {folder}/{name}")
    return folder


def prepare(record, frozen, route_id, chrom, output):
    asset = record["projections"].get(f"{route_id}/{chrom}")
    if asset is None:
        logging.info("No preprocessing projection planned for %s/chr%s; preparing fresh inputs", route_id, chrom)
        return False
    parent = checked_files(asset["panel"])
    old = core.read(parent / "panel.json")
    config = core.read(frozen / "config.json")
    frozen_inputs = core.read(frozen / "manifest.json")["input_sha256"]
    scoped = copy.deepcopy(config)
    scoped.update(chromosomes=[chrom], mode="development_screen")
    for key in ("reference", "training", "maps", "spans"):
        scoped[key] = {str(chrom): scoped[key][str(chrom)]}
    for route in scoped["routes"]:
        route["vcfs"] = {str(chrom): route["vcfs"][str(chrom)]}
    for path in core.input_paths(scoped):
        if core.digest(path) != frozen_inputs[path]:
            raise ValueError(f"Frozen projection input changed: {path}")
        if path != config["samples"] and asset["input_sha256"].get(path) != frozen_inputs[path]:
            raise ValueError(f"Projection source used a different scientific input: {path}")
    for key in ("sources", "seed", "reference_filters", "segments"):
        if config[key] != old["config"][key]:
            raise ValueError(f"Projection changes scientific policy: {key}")
    if (old["chromosome"] != chrom or old["span"] != config["spans"][str(chrom)]
            or old["route"]["id"] != route_id
            or old["route"]["stage"] != next(route["stage"] for route in config["routes"] if route["id"] == route_id)):
        raise ValueError("Projection route/chromosome/span identity differs")
    samples = core.rows(config["samples"])
    reference = [row for row in samples if row["role"] == "reference"]
    training = [row for row in samples if row["role"] == "train"]
    queries = [row for row in samples if row["role"] in ("calibration", "validation")]
    if reference != old["reference"] or training != old["training"]:
        raise ValueError("Projection must preserve complete, identical reference and training roles")
    indices = []
    for sample in queries:
        matches = [index for index, row in enumerate(old["queries"]) if row["sample_id"] == sample["sample_id"]]
        if len(matches) != 1 or old["queries"][matches[0]] != sample:
            raise ValueError("Projected query identity, labels or role ownership differs")
        indices.append(matches[0])
    if core.digest(parent / "data.npz") != old["data_sha256"]:
        raise ValueError("Parent panel genotype data changed")
    output.mkdir()
    with np.load(parent / "data.npz") as data:
        positions = np.asarray([key[1] for key in old["keys"]])
        expected, labels = core.truth(config, training + queries, chrom, positions)
        actual = np.concatenate([data["truth"][:len(training)], data["truth"][len(training) + np.asarray(indices)]])
        if (not np.array_equal(expected, actual)
                or not np.array_equal(labels[:len(training)], data["train_labels"])
                or not np.array_equal(data["keys"], np.asarray(old["keys"], dtype=str))):
            raise ValueError("Projected physical truth, marker axes or training labels differ")
        values = {name: data[name] for name in data.files}
        values.update(query=data["query"][indices], truth=actual)
        np.savez_compressed(output / "data.npz", **values)
    for name in asset["panel"]["files"]:
        if name not in ("data.npz", "panel.json", "execution_reuse.json"):
            shutil.copyfile(parent / name, output / name)
    meta = {**old, "config": config, "queries": queries,
            "route": next(route for route in config["routes"] if route["id"] == route_id),
            "base_plan": models.plan(config), "data_sha256": core.digest(output / "data.npz")}
    core.write(output / "panel.json", meta)
    core.write(output / "execution_reuse.json", {
        "status": "verified_query_subset", "source": asset["panel"],
        "query_indices": indices, "reference_training_unchanged": True,
        "fitted_models_rebound": False,
    })
    return True


def phase(record, panel, output):
    meta = core.read(panel / "panel.json")
    asset = record["projections"].get(f"{meta['route']['id']}/{meta['chromosome']}")
    if asset is None:
        logging.info("No phase projection planned for %s/chr%s; phasing fresh inputs",
                     meta["route"]["id"], meta["chromosome"])
        return False
    parent = checked_files(asset["phase"])
    old = core.read(Path(asset["panel"]["directory"]) / "panel.json")
    source = core.read(parent / "phase.json")
    requires_phase = meta["route"]["stage"] == "harmonized" and any(
        model["family"] in ("flare", "gnomix") for model in meta["config"]["models"])
    if (source["query_sha256"] != core.digest(parent / "query.npy")
            or source["phase_only"] != requires_phase
            or source["added_variants"] != 0 or source["changed_diploid_calls"] != 0):
        raise ValueError("Projected phase provenance or phase-only semantics differ")
    index = {row["sample_id"]: i for i, row in enumerate(old["queries"])}
    calls = np.load(parent / "query.npy", mmap_mode="r")[
        [index[row["sample_id"]] for row in meta["queries"]]]
    with np.load(panel / "data.npz") as data:
        if calls.shape != data["query"].shape or not np.array_equal(calls.sum(axis=1), data["query"].sum(axis=1)):
            raise ValueError("Projected phasing changed observed diploid genotypes")
    output.mkdir()
    for name in asset["phase"]["files"]:
        if name not in ("query.npy", "phase.json", "execution_reuse.json"):
            shutil.copyfile(parent / name, output / name)
    np.save(output / "query.npy", calls)
    core.write(output / "phase.json", {**source, "query_sha256": core.digest(output / "query.npy")})
    core.write(output / "execution_reuse.json", {
        "status": "verified_singleton_phase_subset", "source": asset["phase"],
        "sample_ids": [row["sample_id"] for row in meta["queries"]],
        "changed_diploid_calls": 0,
    })
    return True


def compatible(source):
    current = Path(__file__).resolve().parent
    verified = {}
    groups = {
        "core.py": ("prepare", "phase", "truth", "vcf", "write_vcf", "physical_weights"),
        "models.py": ("reference_subset", "base_policy", "plan", "fit_base", "fit_model",
                      "flare_command", "flare_posterior"),
    }
    for filename, names in groups.items():
        parsed = []
        for folder in (source, current):
            parsed.append({node.name: ast.dump(node, include_attributes=False)
                           for node in ast.parse((folder / filename).read_text()).body
                           if isinstance(node, ast.FunctionDef)})
        for name in names:
            if parsed[0][name] != parsed[1][name]:
                raise ValueError(f"Projection cannot reuse changed scientific producer: {filename}/{name}")
        verified[filename] = list(names)
    return verified


def capture(folder, meta, params, identity, base="base", output="prediction"):
    if (not all((folder / f"{identity}{suffix}").is_file() for suffix in REQUIRED_NATIVE)
            or "End Time" not in (folder / f"{identity}.log").read_text()):
        return None
    expected = models.flare_command(meta, params, Path(base), Path(output) / f"{identity}.vcf.gz",
                                    Path(output) / identity)
    actual = core.read(folder / f"{identity}.flare.command.json")
    if actual[:1] != expected[:1] or actual[2:] != expected[2:]:
        raise ValueError(f"Reusable FLARE command differs: {identity}")
    return {"directory": str(folder), "files": {
        f"{identity}{suffix}": core.digest(folder / f"{identity}{suffix}")
        for suffix in NATIVE_SUFFIXES if (folder / f"{identity}{suffix}").is_file()}}


def build_reuse(config, samples, source_run, output, *, stopped):
    root = source_run if (source_run / "source").is_dir() else source_run.parent
    if (root / "java-profile").read_text().strip() != "c1":
        raise ValueError("Staged reuse requires the same pinned C1 Java runtime")
    for line in (root / "source.sha256").read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        if core.digest(root / "source" / filename) != expected:
            raise ValueError(f"Frozen source changed: {filename}")
    producer_functions = compatible(root / "source")
    results = source_run / ("smoke/results" if stopped else "results")
    if stopped:
        transition = source_run / "reports/staged-scope-transition"
        controller = core.read(transition / "controller-after.json")
        tasks = core.read(transition / "owned-tasks-after.json")
        if controller["State"]["Running"] or any(task["State"]["Running"] for task in tasks):
            raise ValueError("Stop the superseded controller and its owned tasks before projection")
        environment = dict(value.split("=", 1) for value in controller["Config"]["Env"])
        if Path(environment["EXPERIMENT"]).resolve() != source_run.resolve() or any(
                not Path(task["Config"]["WorkingDir"]).is_relative_to(source_run / "smoke/work") for task in tasks):
            raise ValueError("Stop proof does not belong to the declared source run")
    else:
        from audit import verify_manifest
        verify_manifest(results)
    queries = [row for row in samples if row["role"] in ("calibration", "validation", "test")]
    identities = {row["sample_id"] for row in queries}
    source_inputs = core.read(results / "frozen/frozen/manifest.json")["input_sha256"]
    projections, native = {}, {}
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            key = f"{route['id']}/{chrom}"
            panel = results / f"panels/{route['id']}/chr{chrom}/panel"
            available = (panel / "panel.json").is_file()
            if available:
                meta = core.read(panel / "panel.json")
                phased = results / f"phased/{route['id']}/chr{chrom}/phased"
                if not (phased / "phase.json").is_file():
                    raise ValueError(f"Expected complete source phasing is missing: {phased}")
                phase_files = {path.name: core.digest(path) for path in phased.iterdir()
                               if path.is_file() and (path.name in ("query.npy", "phase.json")
                               or any(path.name.startswith(identity + ".") for identity in identities))}
                projections[key] = {
                    "panel": {"directory": str(panel), "files": inventory(panel)},
                    "phase": {"directory": str(phased), "files": phase_files},
                    "input_sha256": source_inputs,
                }
            for spec in config["models"]:
                if spec["family"] != "flare":
                    continue
                job = {"signature": signature(panel, spec["params"]) if available else None, "samples": {}}
                native[f"{key}/{spec['id']}"] = job
                if not available:
                    continue
                if not any(row == spec for row in meta["config"]["models"]):
                    raise ValueError("Projected native model parameters differ from the source config")
                prediction = results / f"predictions/{route['id']}/chr{chrom}/{spec['id']}/prediction"
                if (prediction / "prediction.json").is_file():
                    record = core.read(prediction / "prediction.json")
                    if record["config"] != meta["config"] or record["samples"] != meta["queries"]:
                        raise ValueError("Completed source predictions have inconsistent identities")
                    for identity in identities & {row["sample_id"] for row in meta["queries"]}:
                        found = capture(prediction, meta, spec["params"], identity)
                        if found is None:
                            raise ValueError("Completed prediction lacks required singleton evidence")
                        job["samples"][identity] = found
    if stopped:
        previous_path = source_run / "runtime/reuse-manifest.json"
        if core.digest(previous_path) != (source_run / "runtime/reuse-manifest.sha256").read_text().split()[0]:
            raise ValueError("Original frozen reuse manifest changed")
        previous = core.read(previous_path)
        if previous["config"] != core.read(source_run / "smoke.json"):
            raise ValueError("Previous frozen reuse manifest belongs to a different source config")
        for key, job in native.items():
            route, chromosome, identity = key.split("/")
            panel = results / f"panels/{route}/chr{chromosome}/panel"
            meta = core.read(panel / "panel.json")
            base_id = next(row["base_id"] for row in meta["base_plan"]["models"] if row["id"] == identity)
            from reuse import inference_inputs
            original_inputs = inference_inputs(
                panel, results / f"phased/{route}/chr{chromosome}/phased",
                results / f"bases/{route}/chr{chromosome}/{base_id}/base",
                results / f"models/{route}/chr{chromosome}/{identity}/fitted")
            if previous["native"][key]["input_sha256"] != original_inputs:
                raise ValueError("Original singleton manifest differs from the exact adopted source inputs")
            for identity, item in previous["native"][key]["samples"].items():
                if identity in identities:
                    job["samples"][identity] = item
        workdirs = [Path(task["Config"]["WorkingDir"]) for task in tasks]
        for trace in sorted((source_run / "reports").glob("*.smoke.trace.tsv")):
            with trace.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["name"].startswith("INFER (") and row["status"] == "COMPLETED":
                        matches = list((source_run / "smoke/work").glob(row["hash"] + "*"))
                        if len(matches) != 1:
                            raise ValueError("Source chunk work directory is missing or ambiguous")
                        workdirs.append(matches[0])
        checked_panels, checked_models = set(), set()
        parent_metadata = {key: core.read(Path(asset["panel"]["directory"]) / "panel.json")
                           for key, asset in projections.items()}
        for work in dict.fromkeys(workdirs):
            tokens = shlex.split("\n".join(line for line in (work / ".command.sh").read_text().splitlines()
                                           if not line.startswith("#")))
            if tokens[2] != "infer":
                continue
            args = dict(zip(tokens[3::2], tokens[4::2]))
            spec = core.read(work / args["--fitted"] / "model.json")
            if spec["family"] != "flare":
                continue
            if "bc523b63202a52dd2aad32ccadb7813c179c37c93fe89168bc3ae2c1c0172547" not in (
                    work / ".command.run").read_text():
                raise ValueError("Source singleton used a different Java task runtime")
            meta = core.read(work / args["--panel"] / "panel.json")
            key = f"{meta['route']['id']}/{meta['chromosome']}/{spec['id']}"
            if key not in native:
                continue
            params = next(row["params"] for row in config["models"] if row["id"] == spec["id"])
            parent_key = f"{meta['route']['id']}/{meta['chromosome']}"
            source_data = (work / args["--panel"] / "data.npz").resolve()
            if meta != parent_metadata[parent_key]:
                raise ValueError("Source chunk has different scientific inference inputs")
            if source_data not in checked_panels:
                if core.digest(source_data) != projections[parent_key]["panel"]["files"]["data.npz"]:
                    raise ValueError("Source chunk has different reference/training/query genotypes")
                checked_panels.add(source_data)
            fitted = (work / args["--fitted"] / "model.pkl").resolve()
            if fitted not in checked_models:
                saved = models.load(fitted)
                if (core.digest(fitted) != spec["sha256"] or saved["spec"]["params"] != params
                        or saved["panel_sha256"] != meta["data_sha256"]):
                    raise ValueError("Source chunk fitted parameters or model hash differs")
                checked_models.add(fitted)
            for identity in sorted(identities - native[key]["samples"].keys()):
                found = capture(work / args["--output"], meta, params, identity, args["--base"], args["--output"])
                if found is not None:
                    native[key]["samples"][identity] = found
    core.write(output, {
        "version": 2, "config": config, "java_profile": "c1", "source_run": str(source_run),
        "producer_functions": producer_functions, "projections": projections, "native": native,
        "reuse_policy": "unchanged complete reference/training; exact query-subset phase and native GT validation; refit models",
    })
    reused = sum(len(job["samples"]) for job in native.values())
    total = len(queries) * len(native)
    print(f"Projection plan: {len(projections)} prepared/phased panels; {reused}/{total} existing FLARE singletons; "
          f"{total - reused} planned new calls", flush=True)
