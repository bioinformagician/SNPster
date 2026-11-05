#!/usr/bin/env python3
"""Independent durable-output reconstruction; usable after deleting Nextflow work."""
import argparse
import gzip
from pathlib import Path

import numpy as np

from core import digest, probability, read, rounded_flare_probability, rows, write
from evaluation import verify


def audit(root):
    frozen = root / "frozen/frozen"
    config = read(frozen / "config.json")
    provenance = read(frozen / "manifest.json")
    for path, expected in provenance["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Scientific input/native dependency changed: {path}")
    summary = verify(root / "summary/summary")
    sources = config["sources"]
    segments = rows(config["segments"])
    denominator = sum(config["spans"][str(chrom)][1] - config["spans"][str(chrom)][0]
                      for chrom in config["chromosomes"])
    truth = {}
    for sample in rows(config["samples"]):
        if sample["role"] not in ("calibration", "validation", "test"):
            continue
        value = np.zeros(len(sources))
        for chrom in config["chromosomes"]:
            span = config["spans"][str(chrom)]
            for hap in (0, 1):
                selected = sorted((int(row["start"]), int(row["end"]), row["source"])
                                  for row in segments if row["sample_id"] == sample["sample_id"]
                                  and row["chromosome"] == str(chrom) and row["haplotype"] == str(hap))
                ends = [row[1] for row in selected]
                starts = [row[0] for row in selected]
                if not starts or starts != [span[0]] + ends[:-1] or ends[-1] != span[1]:
                    raise ValueError("Independent physical truth coverage failed")
                for start, end, source in selected:
                    value[sources.index(source)] += (end - start) / (2 * denominator)
        truth[sample["sample_id"]] = probability(value, len(sources))
    checked_local, checked_aggregates = 0, 0
    for route in config["routes"]:
        for model in config["models"]:
            weighted = []
            for chrom in config["chromosomes"]:
                folder = root / "predictions" / route["id"] / f"chr{chrom}" / model["id"] / "prediction"
                meta = read(folder / "prediction.json")
                data = np.load(folder / "predictions.npz")
                panel = root / "panels" / route["id"] / f"chr{chrom}" / "panel"
                panelmeta = read(panel / "panel.json")
                positions = np.asarray([key[1] for key in panelmeta["keys"]], dtype=float)
                span = config["spans"][str(chrom)]
                # Reconstruct the physical integration separately from the inference adapter.
                edges = [float(span[0])] + [(positions[i] + positions[i + 1]) / 2
                                            for i in range(len(positions) - 1)] + [float(span[1])]
                cells = np.diff(edges)
                if not np.all(cells > 0) or cells.sum() != span[1] - span[0]:
                    raise ValueError("Independent midpoint coverage failed")
                if model["family"] == "flare":
                    for index, sample in enumerate(meta["samples"]):
                        native = [line.split() for line in (folder / f"{sample['sample_id']}.model").read_text().splitlines()
                                  if line.strip() and not line.lstrip().startswith("#")][0]
                        integrated, found = np.zeros(len(native)), []
                        header = None
                        with gzip.open(folder / f"{sample['sample_id']}.anc.vcf.gz", "rt") as handle:
                            for line in handle:
                                if line.startswith("#CHROM"):
                                    header = line.rstrip().split("\t")[9:]
                                if line.startswith("#"):
                                    continue
                                fields = line.rstrip().split("\t")
                                found.append([fields[0], int(fields[1]), fields[3], fields[4]])
                                formats = dict(zip(fields[8].split(":"), fields[9].split(":")))
                                for tag in ("ANP1", "ANP2"):
                                    value = rounded_flare_probability(list(map(float, formats[tag].split(","))), len(native))
                                    integrated += cells[len(found) - 1] * value / value.sum() / 2
                        if found != panelmeta["keys"] or header != [sample["sample_id"]]:
                            raise ValueError("Independent native key/sample audit failed")
                        integrated = integrated[[native.index(name) for name in meta["components"]]] / cells.sum()
                        np.testing.assert_allclose(integrated, data["values"][index], rtol=0, atol=1e-10)
                        checked_local += 1
                elif model["family"] == "gnomix":
                    local = np.load(folder / "local_posteriors.npz")
                    windows = [cells[a:b].sum() for a, b in zip(local["starts"], local["ends"])]
                    np.testing.assert_array_equal(windows, local["weights"])
                    probabilities = local["probabilities"].reshape(len(meta["samples"]), 2, len(windows), len(sources))
                    reconstructed = (probabilities * np.asarray(windows)[None, None, :, None]).sum(axis=(1, 2)) / (2 * cells.sum())
                    np.testing.assert_allclose(reconstructed, data["values"], rtol=0, atol=1e-10)
                    checked_local += len(meta["samples"])
                weighted.append(data["values"] * cells.sum())
            aggregate = root / "aggregated" / route["id"] / model["id"] / "aggregate"
            saved = np.load(aggregate / "aggregate.npz")
            np.testing.assert_allclose(np.sum(weighted, axis=0) / denominator, saved["values"], rtol=0, atol=1e-10)
            meta = read(aggregate / "aggregate.json")
            if [row["sample_id"] for row in meta["samples"]] != list(truth):
                raise ValueError("Independent declared-person coverage failed")
            np.testing.assert_allclose(np.stack([truth[row["sample_id"]] for row in meta["samples"]]),
                                       saved["truth"], rtol=0, atol=1e-10)
            if model["params"].get("reference_label") == "fine":
                reference_rows = rows(config["samples"])
                populations = {}
                for row in reference_rows:
                    if row["role"] == "reference":
                        if row["fine"] in populations and populations[row["fine"]] != row["source"]:
                            raise ValueError("Fine-population mapping is ambiguous")
                        populations[row["fine"]] = row["source"]
                if set(meta["components"]) != set(populations):
                    raise ValueError("Native fine-population components lack exact reference coverage")
                mapping = np.asarray([[float(populations[name] == source) for source in sources]
                                      for name in meta["components"]])
                np.testing.assert_array_equal(mapping, model["params"]["component_map"])
                interpreted = root / "interpreters" / route["id"] / model["id"] / "native/interpreted"
                np.testing.assert_allclose(np.load(interpreted / "interpreted.npz")["values"],
                                           saved["values"] @ mapping, rtol=0, atol=1e-10)
            checked_aggregates += 1
    for item in read(root / "summary/summary/prediction_evidence.json"):
        np.testing.assert_allclose(item["truth"], truth[item["sample_id"]], rtol=0, atol=1e-10)
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Published result must survive work-cache deletion: {path}")
        if path.is_file() and path.name != "independent_verification.json":
            files[str(path.relative_to(root))] = digest(path)
    return {**summary, "independent_local_predictions_checked": checked_local,
            "independent_aggregates_checked": checked_aggregates, "files": files,
            "production_promotion": False, "genetic_relatedness_independently_certified": False}


def verify_manifest(root):
    record = read(root / "independent_verification.json")
    for name, expected in record["files"].items():
        if digest(root / name) != expected:
            raise ValueError(f"Durable artifact changed: {name}")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--manifest-only", action="store_true",
                        help="Verify durable hashes without requiring original scientific inputs.")
    args = parser.parse_args()
    result = verify_manifest(args.results) if args.manifest_only else audit(args.results)
    if not args.manifest_only:
        write(args.results / "independent_verification.json", result)
    print(f"Durable verification passed: {result['independent_aggregates_checked']} aggregates, "
          f"{result['independent_local_predictions_checked']} local predictions; scope={result['scope']}")
