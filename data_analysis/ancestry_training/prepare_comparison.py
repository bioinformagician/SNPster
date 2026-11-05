"""Derive a focused comparison or an unchanged-model chr22 smoke from a validated config."""
import argparse
from collections import Counter
import copy
from pathlib import Path

import core


def derive(config, stage):
    core.validate_config(config)
    if config["mode"] != "full":
        raise ValueError("Comparison preparation requires a full-autosome input config")
    result = copy.deepcopy(config)
    if stage == "smoke":
        result.update(mode="smoke", chromosomes=[22], independent_test=False)
        for key in ("reference", "training", "maps", "spans"):
            result[key] = {"22": result[key]["22"]}
        result["routes"] = [{**route, "vcfs": {"22": route["vcfs"]["22"]}} for route in result["routes"]]
    elif stage == "focused16":
        templates = {model["id"]: model for model in config["models"]}
        required = {
            "flare_full_fixed_g010", "flare_full_em_g010",
            "flare_cap16_fixed_g010", "flare_cap16_em_g010",
            "gnomix_fast_w050_c2_0p1", "gnomix_fast_w050_c2_1",
            "frequency_control", "pca_ridge_control",
        }
        if not required <= templates.keys():
            raise ValueError(f"Missing reviewed model templates: {sorted(required - templates.keys())}")
        if {head["kind"] for head in config["interpretations"]} != {"native", "simplex", "ridge", "rf", "xgboost"}:
            raise ValueError("The focused comparison requires the five reviewed output interpretations")
        candidates = []
        for reference in ("full", "cap16"):
            for generations in (2, 5, 10, 20):
                model = copy.deepcopy(templates[f"flare_{reference}_fixed_g010"])
                model["id"] = f"flare_{reference}_fixed_g{generations:03d}"
                model["params"]["gen"] = generations
                candidates.append(model)
            candidates.append(copy.deepcopy(templates[f"flare_{reference}_em_g010"]))
        for window in (.25, .5):
            for strength in ("0p1", "1"):
                model = copy.deepcopy(templates[f"gnomix_fast_w050_c2_{strength}"])
                model["id"] = f"gnomix_fast_w{round(window * 100):03d}_c2_{strength}"
                model["params"]["window_size_cM"] = window
                candidates.append(model)
        candidates.extend(copy.deepcopy(templates[key]) for key in ("frequency_control", "pca_ridge_control"))
        result["models"] = candidates
    else:
        raise ValueError(f"Unknown comparison stage: {stage}")
    core.validate_config(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("focused16", "smoke"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.config.resolve()
    result = derive(core.read(source), args.stage)
    provenance = {
        "stage": args.stage, "source_config": str(source), "source_sha256": core.digest(source),
        "model_counts": dict(Counter(model["family"] for model in result["models"])),
        "interpretations": [head["id"] for head in result["interpretations"]],
        "chromosomes": result["chromosomes"], "routes": [route["id"] for route in result["routes"]],
    }
    metadata = args.output.with_suffix(".provenance.json")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path, value in ((args.output, result), (metadata, provenance)):
        if path.exists():
            if core.read(path) != value:
                raise ValueError(f"Refusing to change a previously generated comparison: {path}")
        else:
            core.write(path, value)
    print(f"{args.stage}: {len(result['models'])} settings; "
          f"{len(result['interpretations'])} interpretations; "
          f"{len(result['routes'])} routes; {len(result['chromosomes'])} chromosomes")


if __name__ == "__main__":
    main()
