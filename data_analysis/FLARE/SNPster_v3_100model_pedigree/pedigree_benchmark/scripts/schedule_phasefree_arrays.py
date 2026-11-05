#!/usr/bin/env python3
"""Run the frozen array validation and verification inside one durable container."""

import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import model_family_gnomix as adapter
import validate_phasefree_arrays as validation


def archive_incomplete_preparation(out, chromosome):
    directory = out / f"chr{chromosome}"
    if not directory.exists() or (directory / "prepared.json").exists():
        return None
    if any(directory.glob("*/completion.json")) or (directory / "completion.json").exists():
        raise ValueError(f"Unprepared chromosome contains completed work: {directory}")
    attempt = 1
    archive = directory.with_name(f"{directory.name}.incomplete_preparation_{attempt}")
    while archive.exists():
        attempt += 1
        archive = directory.with_name(f"{directory.name}.incomplete_preparation_{attempt}")
    directory.rename(archive)
    return dict(chromosome=chromosome, archive=str(archive),
                reason="Incomplete preparation preserved; no completed case checkpoints present.")


def run_batch(commands, out, state):
    workers = []
    with ExitStack() as stack:
        try:
            for label, command in commands:
                log = stack.enter_context((out / f"durable_{label}.log").open("a", buffering=1))
                log.write(f"\nStarting {time.strftime('%Y-%m-%dT%H:%M:%S%z')}: {command}\n")
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                workers.append((label, process))
            results = {label: process.wait() for label, process in workers}
        finally:
            for _, process in workers:
                if process.poll() is None:
                    process.terminate()
            for _, process in workers:
                if process.poll() is None:
                    process.wait()
    failures = {label: code for label, code in results.items() if code}
    if failures:
        state.update(status="failed", failed_jobs=failures,
                     logs=[str(out / f"durable_{label}.log") for label in failures])
        adapter.write_json(out / "durable_scheduler.json", state)
        raise RuntimeError(f"Frozen array jobs failed: {failures}; inspect the recorded logs")


def schedule(out, recover_preparation=False):
    config = validation.configuration(out)
    lock = out / "scheduler.lock"
    lock.mkdir()
    state = dict(
        status="starting", scheduler_pid=os.getpid(),
        experiment_sha256=adapter.digest(out / "experiment.json"),
        model_sha256=config["model_sha256"], implementation_sha256=adapter.digest(__file__),
        parallel_chromosomes=2, inherited_completed_chromosomes=[],
        recovered_preparations=[], completed_chromosomes=[],
        started_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    adapter.write_json(lock / "owner.json", state)
    try:
        experiment_hash = state["experiment_sha256"]
        for chromosome in range(1, 23):
            if validation.verified(out / f"chr{chromosome}/completion.json", experiment_hash):
                state["inherited_completed_chromosomes"].append(chromosome)
                state["completed_chromosomes"].append(chromosome)
            elif recover_preparation:
                archived = archive_incomplete_preparation(out, chromosome)
                if archived:
                    state["recovered_preparations"].append(archived)
        pending = [chromosome for chromosome in range(1, 23)
                   if chromosome not in state["completed_chromosomes"]]
        script = Path(validation.__file__)
        for index in range(0, len(pending), 2):
            batch = pending[index:index + 2]
            state.update(status="running", active_chromosomes=batch)
            adapter.write_json(out / "durable_scheduler.json", state)
            print(f"Starting chromosomes {batch}; completed {state['completed_chromosomes']}", flush=True)
            run_batch([
                (f"chr{chromosome}", [
                    sys.executable, str(script), "--stage", "chromosome",
                    "--chromosome", str(chromosome), "--output", str(out), "--retry-incomplete"])
                for chromosome in batch
            ], out, state)
            for chromosome in batch:
                if not validation.verified(out / f"chr{chromosome}/completion.json", experiment_hash):
                    raise ValueError(f"Worker returned without completed chromosome {chromosome}")
                state["completed_chromosomes"].append(chromosome)
            state["completed_chromosomes"].sort()
            adapter.write_json(out / "durable_scheduler.json", state)
        state.update(status="scoring", active_chromosomes=[])
        adapter.write_json(out / "durable_scheduler.json", state)
        run_batch([("summarize", [sys.executable, str(script), "--stage", "summarize",
                                 "--output", str(out)])], out, state)
        state.update(status="verifying")
        adapter.write_json(out / "durable_scheduler.json", state)
        run_batch([("verify", [sys.executable, str(script.with_name("verify_phasefree_arrays.py")),
                              "--output", str(out)])], out, state)
        verification = json.loads((out / "verification.json").read_text())
        if (verification["status"] != "complete_verified_seen_donor_stress_test"
                or verification["completion_sha256"] != adapter.digest(out / "completion.json")):
            raise ValueError("Final verification is absent or incompatible")
        state.update(status="complete_verified", finished_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                     completion_sha256=adapter.digest(out / "completion.json"),
                     verification_sha256=adapter.digest(out / "verification.json"))
        adapter.write_json(out / "durable_scheduler.json", state)
        print("Frozen array stress tests and independent verification completed.", flush=True)
    finally:
        (lock / "owner.json").unlink()
        lock.rmdir()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("/validation"))
    parser.add_argument("--recover-preparation", action="store_true")
    args = parser.parse_args()
    schedule(args.output, args.recover_preparation)
