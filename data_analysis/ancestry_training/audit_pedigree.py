"""Independently reconstruct founder haplotypes, ancestry and held-out role ownership."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import tempfile

import numpy as np

import core
from audit_kinship import chromosome as audit_exports
from prepare_kinship import ancestry_totals


def chromosome(master, bundle, plan_dir, chrom, output, source_master=None):
    data = core.read(plan_dir / "plan.json")
    samples = core.rows(plan_dir / "samples.tsv")
    people = [row["sample_id"] for row in samples]
    keys, calls = core.vcf(master / "all.phased.vcf.gz", chrom, people, True)
    index = {identity: offset for offset, identity in enumerate(people)}
    nodes = data["nodes"]
    real = [identity for identity in people if nodes[identity]["kind"] == "real_person"]
    source = (source_master / "all.phased.vcf.gz" if source_master is not None
              else Path(data["parent"]) / f"genotypes/all/chr{chrom}/master/all.phased.vcf.gz")
    original_keys, original = core.vcf(source, chrom, real, True)
    if keys != original_keys:
        raise ValueError("Source allele keys changed")
    np.testing.assert_array_equal(calls[[index[identity] for identity in real]], original)
    # Parse raw simulator evidence separately from the producer's BP parser.
    parsed = defaultdict(dict)
    for line in (master / "pedsim.bp").read_text().splitlines():
        fields = line.split()
        identity, hap = data["aliases"][fields[0]], int(fields[2][-1])
        cursor, pieces = None, []
        for token in fields[3:]:
            if "|" in token:
                chromosome_id, start = token.split("|")
                if chromosome_id.removeprefix("chr") != str(chrom):
                    raise ValueError("Unexpected breakpoint chromosome")
                cursor = int(start)
            else:
                donor_hap, end = map(int, token.split(":"))
                if cursor is None or end < cursor:
                    raise ValueError("Invalid raw breakpoint interval")
                pieces.append((cursor, end + 1, donor_hap))
                cursor = end + 1
        if hap in parsed[identity]:
            raise ValueError("Duplicate breakpoint homolog")
        parsed[identity][hap] = pieces
    mapping = {}
    founders = {identity for people_in_role in data["selected_founders"].values() for identity in people_in_role}
    if set(parsed) != set(data["aliases"].values()):
        raise ValueError("Incomplete raw breakpoint evidence")
    for identity in founders:
        for hap in (0, 1):
            ids = {piece[2] for piece in parsed[identity][hap]}
            if len(ids) != 1 or next(iter(ids)) in mapping:
                raise ValueError("Ambiguous founder homolog")
            mapping[ids.pop()] = (identity, hap)
    by_hap = defaultdict(list)
    for row in core.rows(master / "inheritance.tsv"):
        by_hap[(row["sample_id"], int(row["haplotype"]))].append(row)
    positions = np.array([key[1] for key in keys])
    start, end = data["spans"][str(chrom)]
    checked = 0
    for identity in people:
        node = nodes[identity]
        for hap in (0, 1):
            if node["kind"] == "real_person":
                expected = [(start, end, identity, hap)]
            else:
                expected = [(max(start, a), min(end, b), *mapping[hid])
                            for a, b, hid in parsed[identity][hap] if max(start, a) < min(end, b)]
            rows = by_hap[(identity, hap)]
            actual = [(int(r["start"]), int(r["end"]), r["donor_id"], int(r["donor_haplotype"]))
                      for r in rows]
            if actual != expected:
                raise ValueError("Published inheritance differs from raw founder breakpoints")
            cursor = start
            for row, (a, b, donor, donor_hap) in zip(rows, actual):
                if a != cursor or b <= a or row["source"] != nodes[donor]["source"]:
                    raise ValueError("Invalid ancestry/founder tiling")
                cursor = b
                left, right = np.searchsorted(positions, [a, b])
                np.testing.assert_array_equal(calls[index[identity], hap, left:right],
                                              calls[index[donor], donor_hap, left:right])
                checked += int(right - left)
            if cursor != end:
                raise ValueError("Incomplete physical span")
        if node["kind"] == "simulated_person":
            for hap, parent in enumerate((node["father_id"], node["mother_id"])):
                inherited = calls[index[identity], hap]
                if not np.all((inherited == calls[index[parent], 0]) | (inherited == calls[index[parent], 1])):
                    raise ValueError("Non-Mendelian transmitted allele")
    segments = core.rows(master / "segments.tsv")
    if len(segments) != sum(len(rows) for rows in by_hap.values()):
        raise ValueError("Missing or extra ancestry intervals")
    for row, inheritance in zip(segments, core.rows(master / "inheritance.tsv")):
        if any(row[key] != inheritance[key] for key in row):
            raise ValueError("Source truth differs from inherited donor labels")
    audit_exports(master, bundle, plan_dir / "samples.tsv", chrom, output,
                  masks=data.get("masks", ("lonni", "frederik")))
    result = core.read(output)
    result.update(founder_haplotype_alleles_independently_checked=checked,
                  original_real_people_checked=len(real), mendelian_alleles_verified=True,
                  breakpoint_sha256=core.digest(master / "pedsim.bp"))
    core.write(output, result)


def verify(dataset, output):
    descriptor = core.read(dataset / "configs/dataset.json")
    excluded = {"configs/dataset.json", "verification/independent.json"}
    actual_files = {str(path.relative_to(dataset)) for path in dataset.rglob("*")
                    if path.is_file()} - excluded
    if actual_files != set(descriptor["files_sha256"]):
        raise ValueError("Dataset hash inventory is incomplete or contains extra files")
    for relative, expected in descriptor["files_sha256"].items():
        path = (dataset / relative).resolve()
        if not path.is_relative_to(dataset.resolve()) or core.digest(path) != expected:
            raise ValueError(f"Changed or unsafe release file: {relative}")
    parent = Path(descriptor["parent"])
    if core.digest(parent / "configs/dataset.json") != descriptor["parent_descriptor_sha256"]:
        raise ValueError("Kinship source descriptor changed")
    samples = core.rows(dataset / descriptor["samples"])
    core.validate_roles(samples, descriptor["sources"])
    nodes = core.read(dataset / descriptor["individuals"])["individuals"]
    roles = {row["sample_id"]: row["role"] for row in samples}
    owners, original_owners = defaultdict(set), defaultdict(set)
    for identity, node in nodes.items():
        if node["component_role"]:
            owners[node["kinship_component"]].add(node["component_role"])
        if node["source_kinship_component"] and node["component_role"]:
            original_owners[node["source_kinship_component"]].add(node["component_role"])
        if identity in roles:
            if node["role"] != roles[identity]:
                raise ValueError("Metadata role disagreement")
            if node["kind"] == "simulated_person":
                father, mother = (nodes[node[key]] for key in ("father_id", "mother_id"))
                if (father["role"] != node["role"] or mother["role"] != node["role"]
                        or father["generation"] != node["generation"] - 1
                        or mother["generation"] != node["generation"] - 1
                        or set(father["founder_donors"]) & set(mother["founder_donors"])):
                    raise ValueError("Cross-role, inbred or nonadjacent-generation pedigree")
                for source in descriptor["sources"]:
                    expected = (father["expected_ancestry"][source] + mother["expected_ancestry"][source]) / 2
                    if node["expected_ancestry"][source] != expected:
                        raise ValueError("Mendelian expected ancestry differs")
    if any(len(group) != 1 for group in list(owners.values()) + list(original_owners.values())):
        raise ValueError("Known family/kinship component crosses roles")
    pairs = core.rows(parent / "kinship/screen/related.kin0")
    source_descriptor = core.read(parent / "configs/dataset.json")
    source_nodes = core.read(parent / source_descriptor["individuals"])["individuals"]
    for row in pairs:
        left, right = (source_nodes[row[key]]["kinship_component"] for key in ("IID1", "IID2"))
        if left != right:
            raise ValueError("Source kinship edge was not included in its screened component")
        if len(original_owners[left] | original_owners[right]) > 1:
            raise ValueError("Flagged genome-wide relative crosses roles")
    # Compare real parent identities to the retained public pedigree, not a producer-derived table.
    for row in core.rows(parent.parent / "inputs/1000genomes.ped"):
        if row["Individual ID"] in nodes:
            for key, column in (("father_id", "Paternal ID"), ("mother_id", "Maternal ID")):
                expected = None if row[column] in ("", "0") else row[column]
                if nodes[row["Individual ID"]][key] != expected:
                    raise ValueError("Published biological parent changed")
    values = ancestry_totals(core.rows(dataset / descriptor["segments"]), descriptor["sources"], descriptor["spans"])
    if set(values) != set(roles):
        raise ValueError("Truth sample axis differs from the active cohort")
    for identity, vector in values.items():
        np.testing.assert_allclose(list(vector.values()),
                                   [nodes[identity]["realized_ancestry"][s] for s in vector], rtol=0, atol=1e-12)
    for role in ("train", "calibration", "validation", "test"):
        counts = Counter(nodes[identity]["generation"] for identity in roles if roles[identity] == role)
        if dict(counts) != {g: counts[0] // 2**g for g in range(4)} or counts[0] % 8:
            raise ValueError("Halving pedigree counts differ")
        for generation in (1, 2, 3):
            descendants = [identity for identity in roles if roles[identity] == role
                           and nodes[identity]["generation"] == generation]
            if any(sum(v > 1e-12 for v in values[identity].values()) < 2 for identity in descendants):
                raise ValueError("A purported mixed descendant is source-pure")
    audits = [core.read(dataset / f"verification/chromosomes/chr{chrom}.json") for chrom in range(1, 23)]
    if any(not item["verified"] or not item["mendelian_alleles_verified"]
           or item["chromosome"] != chrom for chrom, item in enumerate(audits, 1)):
        raise ValueError("Missing independent full-genotype chromosome audit")
    with tempfile.TemporaryDirectory(prefix="snpster-pedigree-interface-") as temporary:
        root = Path(temporary)
        core.freeze(dataset / "configs/smoke.json", root / "frozen")
        for mask in ("lonni", "frederik"):
            core.prepare(root / "frozen", f"harmonized_{mask}", 22, root / mask)
    core.write(output, {
        "verified": True, "independent_test": False, "models_fitted": False,
        "genotype_chromosomes_independently_checked": 22, "cross_role_known_or_flagged_relationships": 0,
        "pipeline_freeze_and_both_mask_prepare_passed": True, "ancestry_people_checked": len(values),
        "role_counts": dict(Counter(roles.values())), "role_generation_counts": descriptor["role_generation_counts"],
        "family_components_per_role": {role: len({nodes[identity]["kinship_component"]
                                                 for identity in roles if roles[identity] == role})
                                       for role in core.ROLES},
        "dataset_descriptor_sha256": core.digest(dataset / "configs/dataset.json"),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("chromosome", "verify"))
    for name in ("master", "bundle", "plan", "dataset", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    args = parser.parse_args()
    if args.step == "chromosome":
        chromosome(args.master, args.bundle, args.plan, args.chromosome, args.output)
    else:
        verify(args.dataset, args.output)


if __name__ == "__main__":
    main()
