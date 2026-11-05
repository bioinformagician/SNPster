"""Independently check the reusable bundle's genotypes, ancestry and split constraints."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import tempfile

import numpy as np

import core


def chromosome(master, bundle, samples_path, chrom, output, masks=("lonni", "frederik")):
    samples = core.rows(samples_path)
    identities = [row["sample_id"] for row in samples]
    keys, values = core.vcf(master / "all.phased.vcf.gz", chrom, identities, True)
    index = {key: offset for offset, key in enumerate(keys)}
    people = {identity: offset for offset, identity in enumerate(identities)}
    checked = {}
    for name, selected, phased in [
        ("reference", [r["sample_id"] for r in samples if r["role"] == "reference"], True),
        ("training", [r["sample_id"] for r in samples if r["role"] == "train"], True),
        *[(f"query_{mask}", [r["sample_id"] for r in samples if r["role"] in
                             ("calibration", "validation", "test")], False) for mask in masks],
    ]:
        path = bundle / f"{name}.vcf.gz"
        actual_keys, actual = core.vcf(path, chrom, selected, phased)
        if phased and actual_keys != keys:
            raise ValueError("Reference/training lost union markers")
        expected = values[[people[identity] for identity in selected]][:, :, [index[key] for key in actual_keys]]
        if phased:
            np.testing.assert_array_equal(actual, expected)
        else:
            np.testing.assert_array_equal(actual.sum(axis=1), expected.sum(axis=1))
            sites = core.rows(master / f"{name.removeprefix('query_')}.sites.tsv")
            if actual_keys != [(r["chromosome"], int(r["position"]), r["ref"], r["alt"]) for r in sites]:
                raise ValueError("Mask allele keys changed")
        checked[name] = {"people": len(selected), "markers": len(actual_keys), "sha256": core.digest(path)}
    core.write(output, {"verified": True, "chromosome": chrom, "exact_phased_gt_or_diploid_gt_preserved": True,
                        "sources": checked, "samples_sha256": core.digest(samples_path)})


def verify(dataset, output, interfaces=False):
    descriptor = core.read(dataset / "configs/dataset.json")
    for name, expected in descriptor["files_sha256"].items():
        path = (dataset / name).resolve()
        if not path.is_relative_to(dataset.resolve()) or core.digest(path) != expected:
            raise ValueError(f"Changed or unsafe dataset path: {name}")
    samples = core.rows(dataset / descriptor["samples"])
    core.validate_roles(samples, descriptor["sources"])
    manifest = core.read(dataset / descriptor["individuals"])
    nodes = manifest["individuals"]
    roles = {row["sample_id"]: row["role"] for row in samples}
    owners = defaultdict(set)
    for node in nodes.values():
        if node["component_role"]:
            owners[node["kinship_component"]].add(node["component_role"])
        for parent in filter(None, (node["father_id"], node["mother_id"], node["donor_id"])):
            if parent not in nodes or nodes[parent]["kinship_component"] != node["kinship_component"]:
                raise ValueError("Cross-component parent/founder relationship")
    if any(len(values) != 1 for values in owners.values()):
        raise ValueError("Biological family component crosses roles")
    screen = core.read(dataset / "kinship/screen/screen.json")
    pairs = core.rows(dataset / "kinship/screen/related.kin0")
    if core.digest(dataset / "kinship/screen/related.kin0") != screen["kinship_table_sha256"]:
        raise ValueError("Kinship results changed")
    for row in pairs:
        left, right = row["IID1"], row["IID2"]
        if float(row["KINSHIP"]) < .0884 or int(row["NSNP"]) < 10000:
            raise ValueError("Invalid genome-wide kinship evidence")
        if nodes[left]["kinship_component"] != nodes[right]["kinship_component"]:
            raise ValueError("Genotype-related pair not connected")
        if left in roles and right in roles and roles[left] != roles[right]:
            raise ValueError("Cross-role genotype relative retained")
    # Published pedigree rows are checked directly, not through the producer's graph.
    public_path = dataset.parent / "inputs/1000genomes.ped"
    if not public_path.is_file():
        raise ValueError("Retain the dataset run's inputs/1000genomes.ped for independent parentage verification")
    families = defaultdict(set)
    for row in core.rows(public_path):
        if row["Individual ID"] in nodes:
            node = nodes[row["Individual ID"]]
            families[row["Family ID"]].add(node["kinship_component"])
            for key, column in (("father_id", "Paternal ID"), ("mother_id", "Maternal ID")):
                expected = None if row[column] in ("", "0") else row[column]
                if node[key] != expected:
                    raise ValueError("Published immediate parental identity was lost or changed")
    if any(len(group) > 1 for group in families.values()):
        raise ValueError("Published biological family was split")
    by_homolog = defaultdict(list)
    for row in core.rows(dataset / descriptor["segments"]):
        by_homolog[(row["sample_id"], row["chromosome"], row["haplotype"])].append(row)
    sources = descriptor["sources"]
    length = 2 * sum(end - start for start, end in descriptor["spans"].values())
    checked_truth = 0
    for identity, node in nodes.items():
        if node["realized_ancestry"] is None:
            if node["genotype_available"]:
                raise ValueError("Missing candidate realized ancestry")
            continue
        counts = Counter()
        for chrom, (start, end) in descriptor["spans"].items():
            for hap in ("0", "1"):
                cursor = start
                for row in sorted(by_homolog[(identity, chrom, hap)], key=lambda row: int(row["start"])):
                    a, b, source = int(row["start"]), int(row["end"]), row["source"]
                    if a != cursor or b <= a or b > end or source not in sources:
                        raise ValueError("Invalid ancestry segment tiling")
                    counts[source] += b - a
                    cursor = b
                if cursor != end:
                    raise ValueError("Ancestry does not cover both homologs/all autosomes")
        np.testing.assert_allclose([counts[source] / length for source in sources],
                                   [node["realized_ancestry"][source] for source in sources],
                                   rtol=0, atol=1e-12)
        if node["kind"] == "simulated_person":
            parents = [nodes[node[key]]["expected_ancestry"] for key in ("father_id", "mother_id")]
            np.testing.assert_allclose([node["expected_ancestry"][s] for s in sources],
                                       [(parents[0][s] + parents[1][s]) / 2 for s in sources],
                                       rtol=0, atol=1e-12)
        checked_truth += 1
    audits = [core.read(dataset / f"verification/chromosomes/chr{chrom}.json") for chrom in range(1, 23)]
    if any(not row["verified"] or row["chromosome"] != chrom for chrom, row in enumerate(audits, 1)):
        raise ValueError("Missing independent chromosome genotype check")
    if interfaces:
        with tempfile.TemporaryDirectory(prefix="snpster-kinship-interface-") as temporary:
            folder = Path(temporary)
            core.freeze(dataset / "configs/smoke.json", folder / "frozen")
            for mask in ("lonni", "frederik"):
                core.prepare(folder / "frozen", f"harmonized_{mask}", 22, folder / mask)
    core.write(output, {
        "verified": True, "genotype_chromosomes_independently_checked": 22,
        "ancestry_individuals_independently_checked": checked_truth,
        "genome_wide_related_pairs_checked": len(pairs),
        "cross_role_known_or_flagged_relationships": 0,
        "pipeline_freeze_and_both_mask_prepare_passed": interfaces,
        "role_counts": dict(Counter(roles.values())), "independent_test": False,
        "scientific_certification": False, "production_promotion": False,
        "dataset_descriptor_sha256": core.digest(dataset / "configs/dataset.json"),
    })
    print(f"Independent reusable dataset audit passed: {len(samples)} people, {checked_truth} ancestry records", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("chromosome", "verify"))
    for name in ("master", "bundle", "samples", "dataset", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--interfaces", action="store_true")
    args = parser.parse_args()
    if args.step == "chromosome":
        chromosome(args.master, args.bundle, args.samples, args.chromosome, args.output)
    else:
        verify(args.dataset, args.output, args.interfaces)


if __name__ == "__main__":
    main()
