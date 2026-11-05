"""Preserve verified completed FLARE singletons before resuming interrupted chunks."""
import argparse
import copy
from pathlib import Path
import shlex
import shutil

import numpy as np

import core
import models
import projection


def checked_singleton(folder, meta, params, identity, expected):
    record = projection.capture(folder, meta, params, identity)
    if record is None:
        return None
    keys, calls = core.vcf(folder / f"{identity}.vcf.gz", meta["chromosome"], [identity], True)
    if keys != [tuple(key) for key in meta["keys"]] or not np.array_equal(calls[0], expected):
        raise ValueError(f"Recovery singleton genotype/key identity differs: {identity}")
    names = [line.split() for line in (folder / f"{identity}.model").read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")][0]
    components = params.get("components", meta["sources"])
    if len(names) != len(components) or set(names) != set(components):
        raise ValueError(f"Recovery singleton ancestry axis differs: {identity}")
    models.flare_posterior(folder / f"{identity}.anc.vcf.gz", meta["keys"], names, meta["span"], identity)
    return record


def snapshot(record, output):
    source = projection.checked_files(record)
    output.mkdir(parents=True)
    for name, expected in record["files"].items():
        shutil.copyfile(source / name, output / name)
        if core.digest(output / name) != expected:
            raise ValueError(f"Recovery snapshot changed while copying: {name}")
    return {**record, "directory": str(output)}


def build(root, stage, stop_receipt, output):
    root = root.resolve()
    if stage not in ("genome", "screen"):
        raise ValueError("Recovery supports only the declared development inference stage")
    stopped = core.read(stop_receipt)
    if (Path(stopped["root"]) != root or stopped["controller_after"]["State"]["Running"]
            or stopped["owned_workers_after"]):
        raise ValueError("Recovery requires exact proof that controller and owned workers stopped")
    if (root / "java-profile").read_text().strip() != "c1":
        raise ValueError("Recovery requires the unchanged pinned C1 runtime")
    for line in (root / "source.sha256").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if core.digest(root / "source" / name) != expected:
            raise ValueError("Frozen scientific source changed before recovery")
    config = core.read(root / "planning/design/config.json")
    core.validate_config(config)
    original = core.read(root / "planning/design/native-reuse.json")
    if original["version"] != 2 or original["config"] != config:
        raise ValueError("Recovery baseline manifest differs from the frozen scientific configuration")
    results = root / stage / "results"
    frozen = core.read(results / "frozen/frozen/manifest.json")
    if core.read(results / "frozen/frozen/config.json") != config:
        raise ValueError("Recovery stage config differs from planning")
    for path, expected in frozen["input_sha256"].items():
        if core.digest(path) != expected:
            raise ValueError(f"Frozen scientific input changed before recovery: {path}")
    output = output.resolve()
    if not output.is_relative_to(root / "reports") or output.exists():
        raise ValueError("Use a new dedicated recovery directory inside this run's reports")
    panels = {}
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            panel = results / f"panels/{route['id']}/chr{chrom}/panel"
            phased = results / f"phased/{route['id']}/chr{chrom}/phased"
            if not phased.is_dir():
                continue
            meta = core.read(panel / "panel.json")
            if meta["config"] != config or core.digest(panel / "data.npz") != meta["data_sha256"]:
                raise ValueError("Recovery panel/config changed")
            phase = core.read(phased / "phase.json")
            if core.digest(phased / "query.npy") != phase["query_sha256"]:
                raise ValueError("Recovery phased query changed")
            calls = np.load(phased / "query.npy")
            with np.load(panel / "data.npz") as data:
                if not np.array_equal(calls.sum(axis=1), data["query"].sum(axis=1)):
                    raise ValueError("Recovery phasing changed diploid genotypes")
            panels[(route["id"], chrom)] = (panel, meta, calls)
    manifest = copy.deepcopy(original)
    output.mkdir(parents=True)
    recovered = []
    skipped_unfinished = []
    for command_path in sorted((root / stage / "work").glob("*/*/.command.sh")):
        lines = [line for line in command_path.read_text().splitlines()
                 if "workflow.py' infer " in line or "workflow.py infer " in line]
        if not lines:
            continue
        if len(lines) != 1:
            raise ValueError("Ambiguous recovery inference task")
        words = shlex.split(lines[0])
        work = command_path.parent
        if work.resolve() != work:
            raise ValueError("Recovery refuses redirected task work")
        def argument(name):
            return words[words.index(name) + 1]
        folder = work / argument("--output")
        if not any("End Time" in path.read_text() for path in folder.glob("*.log")):
            skipped_unfinished.append(str(work))
            continue
        meta = core.read(work / argument("--panel") / "panel.json")
        if meta["config"] != config:
            raise ValueError("Recovery task belongs to another scientific config")
        fitted = work / argument("--fitted")
        fit = core.read(fitted / "model.json")
        model = next(row for row in config["models"] if row["id"] == fit["id"])
        if model["family"] != "flare":
            continue
        panel, expected_meta, calls = panels[(meta["route"]["id"], meta["chromosome"])]
        if meta != expected_meta or core.digest(fitted / "model.pkl") != fit["sha256"]:
            raise ValueError("Recovery task panel/fitted provenance differs")
        base = work / argument("--base")
        base_record = core.read(base / "base.json")
        fitted_model = models.load(fitted / "model.pkl")
        if (fitted_model["spec"]["params"] != model["params"] or fit["base_id"] != base_record["base_id"]
                or fitted_model["panel_sha256"] != expected_meta["data_sha256"]):
            raise ValueError("Recovery task native parameters/base identity differs")
        for name, expected in base_record["artifact_sha256"].items():
            if core.digest(base / name) != expected:
                raise ValueError("Recovery task reference artifact changed")
        signature = projection.signature(panel, model["params"])
        key = f"{meta['route']['id']}/{meta['chromosome']}/{model['id']}"
        job = manifest["native"][key]
        if job["signature"] is not None and job["signature"] != signature:
            raise ValueError("Recovery native signature differs from the planned job")
        job["signature"] = signature
        start, end = int(argument("--query-start")), int(argument("--query-end"))
        if not 0 <= start < end <= len(meta["queries"]):
            raise ValueError("Recovery task query range differs")
        for index in range(start, end):
            identity = meta["queries"][index]["sample_id"]
            record = checked_singleton(folder, meta, model["params"], identity, calls[index])
            if record is None:
                continue
            existing = job["samples"].get(identity)
            if existing is not None:
                projection.checked_files(existing)
                for suffix in (".vcf.gz", ".model", ".anc.vcf.gz"):
                    name = identity + suffix
                    if existing["files"][name] != record["files"][name]:
                        raise ValueError("Recovery duplicate singleton scientific outputs disagree")
                continue
            destination = output / "singletons" / key / identity
            job["samples"][identity] = snapshot(record, destination)
            recovered.append({"job": key, "sample": identity, "source": str(folder)})
    manifest["recovery"] = {"stop_receipt_sha256": core.digest(stop_receipt), "new_completed_singletons": len(recovered)}
    core.write(output / "native-reuse.json", manifest)
    core.write(output / "recovery.json", {
        "new_completed_singletons": len(recovered),
        "total_reusable_singletons": sum(len(job["samples"]) for job in manifest["native"].values()),
        "recovered": recovered, "work_caches_deleted": False, "scientific_inputs_and_sources_reverified": True,
        "tasks_without_completed_singletons": skipped_unfinished,
        "manifest_sha256": core.digest(output / "native-reuse.json"),
    })
    print(f"Recovery preserved {len(recovered)} additional completed singletons; scientific config and frozen sources unchanged")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--stage", choices=("genome", "screen"), required=True)
    parser.add_argument("--stop-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.root, args.stage, args.stop_receipt, args.output)
