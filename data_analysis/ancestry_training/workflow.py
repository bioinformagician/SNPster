#!/usr/bin/env python3
"""Small task entry points; Nextflow owns scheduling and concurrency."""
import argparse
import logging
from pathlib import Path
import shutil

import core
import evaluation
import models
from reuse import adopt_stage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=(
        "freeze", "prepare", "phase", "base", "fit", "infer", "collect-inference",
        "aggregate", "interpret", "select", "verify"))
    for name in ("config", "frozen", "panel", "phased", "base", "fitted", "aggregate", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--reuse-manifest", type=Path)
    for name in ("route", "base-id", "model", "interpreter"):
        parser.add_argument("--" + name)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--query-start", type=int, default=0)
    parser.add_argument("--query-end", type=int)
    parser.add_argument("--results", type=Path, nargs="+")
    args = parser.parse_args()
    if adopt_stage(args):
        return
    if args.stage == "freeze":
        core.freeze(args.config, args.output)
        core.write(args.output / "plan.json", models.plan(core.read(args.output / "config.json")))
        if args.reuse_manifest:
            shutil.copyfile(args.reuse_manifest, args.output / "execution_reuse.json")
    elif args.stage == "prepare":
        core.prepare(args.frozen, args.route, args.chromosome, args.output)
    elif args.stage == "phase":
        core.phase(args.panel, args.output)
    elif args.stage == "base":
        models.fit_base(args.panel, args.base_id, args.output)
    elif args.stage == "fit":
        models.fit_model(args.panel, args.base, args.model, args.output)
    elif args.stage == "infer":
        models.infer(args.panel, args.phased, args.base, args.fitted, args.output,
                     query_start=args.query_start, query_end=args.query_end, reuse_manifest=args.reuse_manifest)
    elif args.stage == "collect-inference":
        models.collect_inference(args.panel, args.results, args.output)
    elif args.stage == "aggregate":
        evaluation.aggregate(args.results, args.output)
    elif args.stage == "interpret":
        evaluation.interpret(args.aggregate, args.interpreter, args.output)
    elif args.stage == "select":
        evaluation.select(args.results, args.output)
    else:
        core.write(args.output, evaluation.verify(args.aggregate))
    logging.info("%s complete", args.stage)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
