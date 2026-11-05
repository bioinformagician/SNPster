"""Explicit, hash-verified adoption from a stopped, scientifically identical run."""
import argparse
import ast
import csv
import logging
from pathlib import Path
import shlex
import shutil

import numpy as np

from core import digest, read, vcf, write


NATIVE_SUFFIXES = (
    ".vcf.gz", ".vcf.gz.csi", ".vcf.compress.command.json", ".vcf.compress.log",
    ".vcf.index.command.json", ".vcf.index.log", ".flare.command.json", ".flare.log",
    ".model", ".anc.vcf.gz", ".global.anc.gz", ".log",
)
REQUIRED_NATIVE = (".vcf.gz", ".vcf.gz.csi", ".flare.command.json", ".model", ".anc.vcf.gz",
                   ".global.anc.gz", ".log")


def inventory(folder):
    paths = sorted(folder.rglob("*"))
    if any(path.is_symlink() for path in paths):
        raise ValueError(f"Reusable output contains a work-dependent symlink: {folder}")
    return {str(path.relative_to(folder)): digest(path) for path in paths if path.is_file()}


def checked_copy(record, output):
    folder = Path(record["directory"])
    if inventory(folder) != record["files"]:
        raise ValueError(f"Reusable output changed: {folder}")
    shutil.copytree(folder, output)
    if inventory(output) != record["files"]:
        raise ValueError("Adopted output differs from its verified source")


def manifest(path, config):
    record = read(path)
    if record["version"] not in (1, 2) or record["config"] != config:
        raise ValueError("Reuse manifest differs from this exact scientific scope/configuration")
    return record


def inference_inputs(panel, phased, base, fitted):
    return {
        "panel.json": digest(panel / "panel.json"), "data.npz": digest(panel / "data.npz"),
        "phase.json": digest(phased / "phase.json"), "query.npy": digest(phased / "query.npy"),
        "model.json": digest(fitted / "model.json"), "model.pkl": digest(fitted / "model.pkl"),
        "base.json": digest(base / "base.json"),
    }


def adopt_stage(args):
    if args.reuse_manifest is None or args.stage not in ("prepare", "phase", "base", "fit", "infer"):
        return False
    if args.stage == "prepare":
        config = read(args.frozen / "config.json")
        route, chrom = args.route, args.chromosome
    else:
        meta = read(args.panel / "panel.json")
        config, route, chrom = meta["config"], meta["route"]["id"], meta["chromosome"]
    record = manifest(args.reuse_manifest, config)
    if record["version"] == 2:
        import projection
        if args.stage == "prepare":
            return projection.prepare(record, args.frozen, route, chrom, args.output)
        if args.stage == "phase":
            return projection.phase(record, args.panel, args.output)
        logging.info("Projection preserves preprocessing/native evidence, not fitted objects; executing %s", args.stage)
        return False
    key = f"{args.stage}/{route}/{chrom}"
    if args.stage == "base":
        key += f"/{args.base_id}"
    elif args.stage == "fit":
        key += f"/{args.model}"
    elif args.stage == "infer":
        fitted = read(args.fitted / "model.json")
        if fitted["family"] == "flare":
            return False
        key += f"/{fitted['id']}"
    if key not in record["assets"]:
        logging.info("No planned reusable asset for %s; executing it normally", key)
        return False
    asset = record["assets"][key]
    checked_copy(asset, args.output)
    write(args.output / "execution_reuse.json", {
        "status": "adopted_not_refitted", "asset": key, "source": asset,
        "manifest_sha256": digest(args.reuse_manifest), "source_run": record["source_run"],
    })
    logging.info("Adopted hash-verified completed output: %s", key)
    return True


def native_job(path, meta, spec, inputs, *, panel=None, params=None):
    if path is None:
        return None
    record = manifest(path, meta["config"])
    key = f"{meta['route']['id']}/{meta['chromosome']}/{spec['id']}"
    if key not in record["native"]:
        raise ValueError("Reuse manifest is missing the declared singleton model/mask job")
    job = record["native"][key]
    if record["version"] == 2:
        from projection import signature
        if job["signature"] is None:
            if job["samples"]:
                raise ValueError("Planned fresh inference unexpectedly contains reusable samples")
            logging.info("Executing planned fresh singleton job: %s", key)
            return None
        if panel is None or params is None or signature(panel, params) != job["signature"]:
            raise ValueError("Projected singleton reference, parameters, marker axes or software differs")
        return job
    if job["input_sha256"] != inputs:
        raise ValueError("Reusable singleton inference has different panel, phase, base or fitted inputs")
    return job


def copy_native(record, output, sample_id, meta, expected):
    source = Path(record["directory"])
    for name, expected_hash in record["files"].items():
        if Path(name).name != name or not name.startswith(sample_id + "."):
            raise ValueError("Invalid reusable singleton artifact name")
        if digest(source / name) != expected_hash:
            raise ValueError(f"Reusable singleton artifact changed: {name}")
        shutil.copyfile(source / name, output / name)
    keys, calls = vcf(output / f"{sample_id}.vcf.gz", meta["chromosome"], [sample_id], True)
    if keys != [tuple(key) for key in meta["keys"]] or not np.array_equal(calls[0], expected):
        raise ValueError("Reusable singleton target differs from the exact phased query")


def producer_compatibility(old_source, new_source):
    for name in ("core.py", "evaluation.py", "audit.py"):
        if digest(old_source / name) != digest(new_source / name):
            raise ValueError(f"Cannot adopt outputs after scientific producer changes: {name}")
    functions = ("base_policy", "plan", "reference_subset", "layout", "gnomix_classes", "pca_features",
                 "transform", "fit_base", "fit_model", "regression", "flare_posterior", "dump", "load")
    old = {node.name: ast.dump(node, include_attributes=False)
           for node in ast.parse((old_source / "models.py").read_text()).body if isinstance(node, ast.FunctionDef)}
    new = {node.name: ast.dump(node, include_attributes=False)
           for node in ast.parse((new_source / "models.py").read_text()).body if isinstance(node, ast.FunctionDef)}
    if any(old[name] != new[name] for name in functions):
        raise ValueError("Cannot adopt outputs after reference fitting, training or posterior-integration changes")
    return list(functions)


def build(old_run, stop_record, trace, output):
    from models import flare_command
    stopped = read(stop_record)
    if stopped["controller"]["State"]["Running"] or stopped["running_task_workdirs"]:
        raise ValueError("Stop the source controller and its owned tasks before freezing reuse evidence")
    environment = dict(value.split("=", 1) for value in stopped["controller"]["Config"]["Env"])
    if Path(environment["EXPERIMENT"]).resolve() != old_run.resolve():
        raise ValueError("Stopped controller belongs to a different run")
    new_source, old_source = Path(__file__).resolve().parent, old_run / "source"
    compatible = producer_compatibility(old_source, new_source)
    source_hashes = {}
    for line in (old_run / "source.sha256").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        path = old_source / name
        if digest(path) != expected:
            raise ValueError(f"Frozen source changed: {path}")
        source_hashes[name] = expected
    config = read(old_run / "smoke.json")
    results = old_run / "smoke/results"
    assets, native, workdirs = {}, {}, []
    for route in config["routes"]:
        for chrom in config["chromosomes"]:
            root = f"{route['id']}/{chrom}"
            panel = results / f"panels/{route['id']}/chr{chrom}/panel"
            meta = read(panel / "panel.json")
            if meta["config"] != config:
                raise ValueError("Source panel does not match the original frozen smoke")
            planned = {
                f"prepare/{root}": panel,
                f"phase/{root}": results / f"phased/{route['id']}/chr{chrom}/phased",
            }
            for base_id in meta["base_plan"]["bases"]:
                planned[f"base/{root}/{base_id}"] = results / f"bases/{route['id']}/chr{chrom}/{base_id}/base"
            for model in config["models"]:
                planned[f"fit/{root}/{model['id']}"] = results / f"models/{route['id']}/chr{chrom}/{model['id']}/fitted"
                if model["family"] != "flare":
                    planned[f"infer/{root}/{model['id']}"] = results / f"predictions/{route['id']}/chr{chrom}/{model['id']}/prediction"
            for key, folder in planned.items():
                assets[key] = {"directory": str(folder), "files": inventory(folder)}
                if not assets[key]["files"]:
                    raise ValueError(f"Missing expected reusable output: {folder}")
    with trace.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["name"].startswith("INFER (") and row["status"] == "COMPLETED":
                matches = list((old_run / "smoke/work").glob(row["hash"] + "*"))
                if len(matches) != 1:
                    raise ValueError("Completed source inference work directory is ambiguous")
                workdirs.append(matches[0])
    workdirs += [Path(row["Config"]["WorkingDir"]) for row in stopped["previously_running_tasks"]]
    skipped = []
    for workdir in workdirs:
        command = shlex.split("\n".join(line for line in (workdir / ".command.sh").read_text().splitlines()
                                       if not line.startswith("#")))
        if command[2] != "infer":
            raise ValueError("Source task is not singleton inference")
        arguments = dict(zip(command[3::2], command[4::2]))
        panel, phased, base, fitted = [workdir / arguments[name] for name in ("--panel", "--phased", "--base", "--fitted")]
        meta, spec = read(panel / "panel.json"), read(fitted / "model.json")
        if spec["family"] != "flare":
            continue
        if (old_run / "java-profile").read_text().strip() != "c1" or "bc523b63202a52dd2aad32ccadb7813c179c37c93fe89168bc3ae2c1c0172547" not in (workdir / ".command.run").read_text():
            raise ValueError("Reuse requires the same pinned C1 Java profile")
        if meta["config"] != config:
            raise ValueError("Source singleton belongs to a different scientific scope")
        key = f"{meta['route']['id']}/{meta['chromosome']}/{spec['id']}"
        if key in native:
            raise ValueError("Ambiguous duplicate singleton source jobs")
        folder = workdir / arguments["--output"]
        job = {"input_sha256": inference_inputs(panel, phased, base, fitted), "samples": {}}
        for sample in meta["queries"]:
            identity = sample["sample_id"]
            if (not all((folder / f"{identity}{suffix}").is_file() for suffix in REQUIRED_NATIVE)
                    or "End Time" not in (folder / f"{identity}.log").read_text()):
                skipped.append({"job": key, "sample_id": identity, "reason": "singleton_not_finished_before_stop"})
                continue
            expected_command = flare_command(
                meta, next(model for model in config["models"] if model["id"] == spec["id"])["params"],
                Path(arguments["--base"]), Path(arguments["--output"]) / f"{identity}.vcf.gz",
                Path(arguments["--output"]) / identity)
            actual_command = read(folder / f"{identity}.flare.command.json")
            if actual_command[:1] != expected_command[:1] or actual_command[2:] != expected_command[2:]:
                raise ValueError(f"Source singleton flags differ: {key}/{identity}")
            files = {f"{identity}{suffix}": digest(folder / f"{identity}{suffix}") for suffix in NATIVE_SUFFIXES
                     if (folder / f"{identity}{suffix}").is_file()}
            job["samples"][identity] = {"directory": str(folder), "files": files}
        native[key] = job
    expected_native = len(config["routes"]) * len(config["chromosomes"]) * sum(model["family"] == "flare" for model in config["models"])
    if len(native) != expected_native:
        raise ValueError("Incomplete source model/mask job coverage")
    result = {
        "version": 1, "config": config, "source_run": str(old_run), "java_profile": "c1",
        "source_sha256": source_hashes, "compatible_producer_functions": compatible,
        "stop_record_sha256": digest(stop_record), "trace_sha256": digest(trace),
        "assets": assets, "native": native, "not_finished": skipped,
    }
    write(output, result)
    print(f"Frozen {len(assets)} completed reusable stages and "
          f"{sum(len(job['samples']) for job in native.values())} completed singleton predictions; "
          f"{len(skipped)} unfinished singletons will execute normally", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-run", required=True, type=Path)
    parser.add_argument("--stop-record", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    build(args.old_run, args.stop_record, args.trace, args.output)
