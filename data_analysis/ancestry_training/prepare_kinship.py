"""Build a portable observed-genotype dataset with pedigree/kinship-component roles."""
import argparse
from collections import Counter, defaultdict
import copy
import hashlib
from pathlib import Path
import subprocess

import numpy as np

import core
from prepare_observed import components, write_rows


SAMPLE_FIELDS = ["sample_id", "role", "family_id", "donors", "source", "track", "fine"]
SEGMENT_FIELDS = ["sample_id", "chromosome", "haplotype", "start", "end", "source"]
MASKS = ("lonni", "frederik")


def collapse(source):
    return "EAS_SEA" if source in ("EAS", "SEA") else source


def pedigree_closure(identities, nodes):
    selected, visiting = set(), set()

    def visit(identity):
        if identity in visiting:
            raise ValueError(f"Cyclic synthetic pedigree: {identity}")
        if identity in selected:
            return
        if identity not in nodes:
            raise ValueError(f"Missing synthetic parent: {identity}")
        visiting.add(identity)
        node = nodes[identity]
        parents = [node["father_id"], node["mother_id"]]
        if bool(parents[0]) != bool(parents[1]):
            raise ValueError(f"Incomplete synthetic parents: {identity}")
        for parent in filter(None, parents):
            visit(parent)
        visiting.remove(identity)
        selected.add(identity)

    for identity in identities:
        visit(identity)
    return selected


def onehot(source, sources):
    if source not in sources:
        raise ValueError(f"Unrecognized ancestry label: {source}")
    return {name: float(name == source) for name in sources}


def expected_vector(row, sources):
    values = dict.fromkeys(sources, 0.)
    for source in ("EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"):
        values[collapse(source)] += float(row[f"expected_{source}"])
    core.probability(list(values.values()), len(sources))
    return values


def catalog(parent, pedigree, authority, output):
    old = parent / "prepare/results/metadata/plan"
    plan = core.read(old / "plan.json")
    project = Path(plan["source"]).parent
    config = core.read(parent / "prepare/results/ready/full.json")
    config["seed"] = 20261025
    verification = core.read(parent / "full/results/independent_verification.json")
    if not verification["verified"] or verification["scope"] != "full":
        raise ValueError("Require the independently verified completed parent experiment")
    frozen = core.read(parent / "full/results/frozen/frozen/manifest.json")
    for path in (old / "samples.tsv", old / "segments.tsv"):
        if core.digest(path) != frozen["input_sha256"][str(path)]:
            raise ValueError(f"Changed parent scientific metadata: {path}")
    paths = [old / "plan.json", old / "samples.tsv", old / "segments.tsv",
             parent / "prepare/results/ready/full.json", pedigree, authority]
    labels = {}
    for row in core.rows(authority):
        if row["role"] == "source":
            value = (collapse(row["source_group"]), row["fine"])
            if row["sample"] in labels and labels[row["sample"]] != value:
                raise ValueError("Conflicting authoritative source labels")
            labels[row["sample"]] = value
    simulated = {}
    for folder in ("results_v3_pedigree_training_screen_01/cohort/pedigree/plans/plans",
                   "results_v3_pedigree_pilot_10/pedigree/plans/plans"):
        path = project / folder / "expected_ancestry.tsv"
        paths.append(path)
        for row in core.rows(path):
            if row["sample_id"] in simulated:
                raise ValueError("Duplicate synthetic pedigree nodes")
            simulated[row["sample_id"]] = row
    candidates = {row["sample_id"]: dict(row) for row in plan["samples"]}
    closure = pedigree_closure([identity for identity in candidates if identity in simulated], simulated)
    founder_cache = {}

    def founder_donors(identity):
        if identity not in founder_cache:
            row = simulated[identity]
            founder_cache[identity] = ({row["founder_donor"]} if row["founder_donor"] else
                                       founder_donors(row["father_id"]) | founder_donors(row["mother_id"]))
        return founder_cache[identity]

    for identity in closure & candidates.keys():
        candidates[identity]["donors"] = ",".join(sorted(founder_donors(identity)))
    donors = {simulated[identity]["founder_donor"] for identity in closure
              if simulated[identity]["founder_donor"]}
    for identity in sorted(donors - candidates.keys()):
        source, fine = labels[identity]
        candidates[identity] = {"sample_id": identity, "role": "support", "family_id": "",
                                "donors": identity, "source": source, "fine": fine, "track": "pure"}
    nodes = {}
    for identity, row in candidates.items():
        if identity in simulated:
            continue
        source, fine = labels[identity]
        if source != row["source"] or (row["fine"] and row["fine"] != fine):
            raise ValueError(f"Pure source authority disagreement: {identity}")
        row["fine"] = fine
        nodes[identity] = {
            "sample_id": identity, "kind": "real_person", "father_id": None, "mother_id": None,
            "donor_id": None, "source": source, "fine": fine,
            "expected_ancestry": onehot(source, config["sources"]),
            "ancestry_basis": "operational_population_label_not_measured_individual_admixture",
            "genotype_available": True,
        }
    for identity in sorted(closure):
        row = simulated[identity]
        expected = expected_vector(row, config["sources"])
        nodes[identity] = {
            "sample_id": identity, "kind": "simulated_founder_alias" if row["founder_donor"] else "simulated_person",
            "father_id": row["father_id"] or None, "mother_id": row["mother_id"] or None,
            "donor_id": row["founder_donor"] or None, "source": collapse(row["founder_source"]),
            "fine": labels[row["founder_donor"]][1] if row["founder_donor"] else None,
            "expected_ancestry": expected, "ancestry_basis": "pedigree_mendelian_expectation",
            "genotype_available": identity in candidates,
            "simulation_family_id": row["family_id"], "generation": int(row["generation"]),
            "founder_donors": sorted(founder_donors(identity)),
            "parentage_evidence": "original_simulation_pedigree_plan",
        }
        if not row["founder_donor"]:
            parents = [expected_vector(simulated[row[key]], config["sources"])
                       for key in ("father_id", "mother_id")]
            if not np.allclose(list(expected.values()),
                               [(parents[0][s] + parents[1][s]) / 2 for s in config["sources"]],
                               rtol=0, atol=1e-12):
                raise ValueError(f"Parent expectation mismatch: {identity}")
    public = core.rows(pedigree)
    public_ids = {row["Individual ID"]: row for row in public}
    if len(public_ids) != len(public):
        raise ValueError("Duplicated public pedigree identities")
    active_families = {public_ids[identity]["Family ID"] for identity in nodes if identity in public_ids}
    family_groups = defaultdict(list)
    for row in public:
        if row["Family ID"] not in active_families:
            continue
        identity = row["Individual ID"]
        if identity not in nodes:
            label = labels.get(identity)
            nodes[identity] = {
                "sample_id": identity, "kind": "real_person", "donor_id": None,
                "source": label[0] if label else None, "fine": label[1] if label else None,
                "expected_ancestry": onehot(label[0], config["sources"]) if label else None,
                "ancestry_basis": "operational_population_label" if label else "unknown",
                "genotype_available": False,
            }
        nodes[identity].update({
            "father_id": None if row["Paternal ID"] in ("0", "") else row["Paternal ID"],
            "mother_id": None if row["Maternal ID"] in ("0", "") else row["Maternal ID"],
            "public_family_id": row["Family ID"], "parentage_evidence": "1000_genomes_public_pedigree",
        })
        family_groups[row["Family ID"]].append(identity)
    groups = list(family_groups.values())
    for identity, node in list(nodes.items()):
        for other in filter(None, (node.get("father_id"), node.get("mother_id"), node.get("donor_id"))):
            groups.append([identity, other])
            if other not in nodes:
                label = labels.get(other)
                nodes[other] = {
                    "sample_id": other, "kind": "real_person", "father_id": None, "mother_id": None,
                    "donor_id": None, "source": label[0] if label else None,
                    "fine": label[1] if label else None, "genotype_available": False,
                    "expected_ancestry": onehot(label[0], config["sources"]) if label else None,
                    "ancestry_basis": "operational_population_label" if label else "unknown",
                }
    segments = []
    truth_audit_path = project / "results_v3_model_family_screen_05/input_audit.json"
    truth_audit = core.read(truth_audit_path)
    paths.append(truth_audit_path)
    segment_paths = sorted(path for path in plan["input_sha256"] if path.endswith(".segments.tsv"))
    for filename in segment_paths:
        path = Path(filename)
        if core.digest(path) != truth_audit["input_sha256"]["/results/inputs/" + path.name]:
            raise ValueError(f"Changed physical ancestry truth: {path}")
        paths.append(path)
        for row in core.rows(path):
            if row["sample_id"] in nodes:
                segments.append(dict(zip(SEGMENT_FIELDS, (
                    row["sample_id"], row["chrom"], str(int(row["haplotype"]) - 1),
                    row["start_bp"], row["end_bp"], collapse(row["source_group"])))))
    for identity, node in nodes.items():
        if node["kind"] in ("real_person", "simulated_founder_alias") and node["source"]:
            for chrom in config["chromosomes"]:
                start, end = config["spans"][str(chrom)]
                for hap in (0, 1):
                    segments.append(dict(zip(SEGMENT_FIELDS, (identity, str(chrom), str(hap),
                                                             str(start), str(end), node["source"]))))
    output.mkdir()
    write_rows(output / "segments.tsv", segments, SEGMENT_FIELDS)
    realized = ancestry_totals(segments, config["sources"], config["spans"])
    for identity, node in nodes.items():
        node["realized_ancestry"] = realized.get(identity)
        if node["genotype_available"] and node["realized_ancestry"] is None:
            raise ValueError(f"Missing ancestry truth for genotype candidate: {identity}")
    for row in plan["legacy"]:
        if not np.allclose(list(realized[row["sample_id"]].values()), row["truth"], rtol=0, atol=1e-10):
            raise ValueError("Realized truth differs from certified legacy ancestry")
    core.write(output / "catalog.json", {
        "version": 1, "parent": str(parent), "source": plan["source"], "template": config,
        "candidates": candidates, "individuals": nodes, "pedigree_groups": groups,
        "extra_donors": sorted(donors - {row["sample_id"] for row in plan["samples"]}),
        "batches": plan["batches"], "legacy": plan["legacy"],
        "input_sha256": {str(path.resolve()): core.digest(path) for path in paths},
        "selection_uses_model_predictions": False,
    })
    print(f"Catalog: {len(candidates)} genotype candidates; {len(nodes)} individuals including parents", flush=True)


def ancestry_totals(segments, sources, spans):
    grouped = defaultdict(list)
    for row in segments:
        grouped[(row["sample_id"], row["chromosome"], row["haplotype"])].append(
            (int(row["start"]), int(row["end"]), row["source"]))
    total = 2 * sum(end - start for start, end in spans.values())
    result = {}
    for identity in sorted({key[0] for key in grouped}):
        counts = dict.fromkeys(sources, 0)
        for chrom, span in spans.items():
            for hap in ("0", "1"):
                cursor = span[0]
                for start, end, source in sorted(grouped[(identity, chrom, hap)]):
                    if start != cursor or not start < end <= span[1] or source not in sources:
                        raise ValueError(f"Invalid ancestry tiling: {identity}/{chrom}/{hap}")
                    counts[source] += end - start
                    cursor = end
                if cursor != span[1]:
                    raise ValueError(f"Incomplete ancestry tiling: {identity}/{chrom}/{hap}")
        result[identity] = {source: counts[source] / total for source in sources}
    return result


def chromosome(catalog_dir, chrom, dense_root, output):
    data = core.read(catalog_dir / "catalog.json")
    source = Path(data["source"])
    identities = sorted(data["candidates"])
    assembly = output.parent / "assembly"
    assembly.mkdir()
    output.mkdir()
    hashes = {}

    def index(path, name):
        core.run(["bcftools", "index", path], assembly / (name + ".index"))

    paths, panels = [], {}
    for mask in MASKS:
        task = source / f"tasks/{mask}_chr{chrom}"
        certificate = core.read(task / "inputs_complete.json")
        if certificate["status"] != "complete":
            raise ValueError("Uncertified observed input")
        parts = [task / "reference.vcf.gz", task / "pure_train.vcf.gz"]
        parts += [task / f"{batch['batch']}.vcf.gz" for batch in data["batches"]]
        parts += [task / f"cases/{row['sample_id']}/inherited.vcf.gz" for row in data["legacy"]]
        for path in parts:
            expected = certificate["output_sha256"]["/benchmark/" + str(path.relative_to(source))]
            if core.digest(path) != expected:
                raise ValueError(f"Observed input hash failed: {path}")
            hashes[str(path)] = expected
        merged = assembly / f"{mask}.vcf.gz"
        core.run(["bcftools", "merge", "--merge", "none", "-Oz", "-o", merged, *parts],
                 assembly / (mask + ".merge"))
        index(merged, mask)
        existing = sorted(set(identities) - set(data["extra_donors"]))
        keys, calls = core.vcf(merged, chrom, existing, True)
        panels[mask] = keys
        del calls
        paths.append(merged)
    legacy = core.read(Path(data["parent"]) / "prepare/results/metadata/plan/plan.json")
    experiment = core.read(Path(legacy["source"]).parent / "results_v3_phasefree_array_validation_01/experiment.json")
    dense = dense_root / f"hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"
    expected = experiment["reference_sha256"]["/imputer/beagle_references/" + dense.name]
    if core.digest(dense) != expected:
        raise ValueError(f"Dense founder resource hash failed: {dense}")
    hashes[str(dense)] = expected
    union_keys = sorted(set(panels["lonni"]) | set(panels["frederik"]), key=lambda key: key[1])
    union = assembly / "union.vcf.gz"
    core.run(["bcftools", "concat", "-a", "-d", "exact", "-Oz", "-o", union, *paths],
             assembly / "union")
    index(union, "union")
    regions = assembly / "regions.tsv"
    regions.write_text("".join(f"{key[0]}\t{key[1]}\t{key[1]}\n" for key in union_keys))
    raw_extra, extra = assembly / "extra.raw.vcf.gz", assembly / "extra.vcf.gz"
    core.run(["bcftools", "view", "-R", regions, "-s", ",".join(data["extra_donors"]),
              "-m2", "-M2", "-v", "snps", "-Oz", "-o", raw_extra, dense], assembly / "extra")
    index(raw_extra, "extra.raw")
    core.run(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz", "-o", extra,
              raw_extra, union], assembly / "extra.exact")
    index(extra, "extra")
    extra_keys, extra_calls = core.vcf(extra, chrom, data["extra_donors"], True)
    if extra_keys != union_keys:
        raise ValueError("External founder genotypes differ from exact union allele keys")
    del extra_calls
    existing = sorted(set(identities) - set(data["extra_donors"]))
    left_keys, left = core.vcf(paths[0], chrom, existing, True)
    right_keys, right = core.vcf(paths[1], chrom, existing, True)
    li, ri = {key: i for i, key in enumerate(left_keys)}, {key: i for i, key in enumerate(right_keys)}
    common = sorted(set(left_keys) & set(right_keys))
    if not common or not np.array_equal(left[:, :, [li[key] for key in common]],
                                         right[:, :, [ri[key] for key in common]]):
        raise ValueError("Cross-mask inherited GT or phase disagreement")
    del left, right
    master = output / "all.phased.vcf.gz"
    core.run(["bcftools", "merge", "--merge", "none", "-Oz", "-o", master, union, extra],
             assembly / "master")
    index(master, "master")
    actual_keys, calls = core.vcf(master, chrom, identities, True)
    if actual_keys != union_keys:
        raise ValueError("Master genotypes lost exact observed keys")
    del calls
    for mask in MASKS:
        write_rows(output / f"{mask}.sites.tsv", [
            dict(zip(["chromosome", "position", "ref", "alt"], key)) for key in panels[mask]],
            ["chromosome", "position", "ref", "alt"])
    prefix = output / "screen"
    plink = ["plink2", "--threads", "1", "--memory", "6144"]
    core.run(plink + ["--vcf", master, "--double-id", "--set-all-var-ids", "@:#:$r:$a",
                     "--maf", ".05", "--geno", "0", "--make-pgen", "--out", prefix],
             assembly / "plink-import")
    core.run(plink + ["--pfile", prefix, "--indep-pairwise", "200kb", "1", ".2",
                     "--out", output / "ld"], assembly / "plink-ld")
    core.run(plink + ["--pfile", prefix, "--extract", output / "ld.prune.in",
                     "--make-pgen", "--out", output / "pruned"], assembly / "plink-pruned")
    for suffix in (".pgen", ".pvar", ".psam"):
        prefix.with_suffix(suffix).unlink()
    core.write(output / "audit.json", {
        "chromosome": chrom, "people": len(identities), "markers": len(union_keys),
        "pruned_markers": len((output / "ld.prune.in").read_text().splitlines()),
        "input_sha256": hashes, "cross_mask_gt_phase_equal": True, "missing_gt_allowed": False,
        "output_sha256": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
    })
    print(f"chr{chrom}: {len(identities)} complete phased people, {len(union_keys)} observed union SNPs", flush=True)


def kinship(bundles, output):
    output.mkdir()
    records = sorted((core.read(path / "audit.json")["chromosome"], path) for path in bundles)
    if [chrom for chrom, _ in records] != list(range(1, 23)):
        raise ValueError("Genome-wide kinship requires all 22 autosomes")
    listing = output / "merge.list"
    listing.write_text("".join(f"{path.resolve() / 'pruned'}\n" for _, path in records))
    command = ["plink2", "--threads", "1", "--memory", "6144"]
    core.run(command + ["--pmerge-list", listing, "--make-pgen", "--out", output / "autosomes"],
             output / "merge")
    variants = sum(1 for line in (output / "autosomes.pvar").read_text().splitlines()
                   if not line.startswith("#"))
    if variants < 10000:
        raise ValueError(f"Insufficient genome-wide LD-pruned markers: {variants}")
    core.run(command + ["--pfile", output / "autosomes", "--make-king-table", "counts",
                        "--king-table-filter", ".0884", "--out", output / "related"],
             output / "king")
    core.write(output / "screen.json", {
        "method": "PLINK2 KING-robust between-family all-pairs",
        "threshold": .0884, "scope": "second_degree_or_closer_screen_not_absolute_unrelatedness",
        "chromosomes": list(range(1, 23)), "markers": variants,
        "maf": .05, "ld_pruning": "200kb step 1 r2 0.2, pooled candidates",
        "tool_version": subprocess.check_output(["plink2", "--version"], text=True).strip(),
        "kinship_table_sha256": core.digest(output / "related.kin0"),
    })
    print(f"Genome-wide KING completed on {variants} LD-pruned SNPs", flush=True)


def choose_roles(groups, candidates, seed, sources):
    role_by_component = {}
    eligible = [(index, [person for person in group if person in candidates])
                for index, group in enumerate(groups)]
    eligible = [(index, people) for index, people in eligible if people]

    def rank(index):
        return hashlib.sha256(f"{seed}:{','.join(groups[index])}".encode()).hexdigest()

    mixed = [(index, people) for index, people in eligible
             if any(candidates[person]["track"] == "mixed" for person in people)]
    if len(mixed) < 4:
        raise ValueError("Fewer than four independent mixed-family components; new simulations are required")
    mixed.sort(key=lambda item: (-sum(candidates[p]["track"] == "mixed" for p in item[1]), rank(item[0])))
    role_by_component[mixed[0][0]] = "train"
    for (index, _), role in zip(mixed[1:4], ("calibration", "validation", "test")):
        role_by_component[index] = role
    for index, _ in mixed[4:]:
        counts = Counter(role_by_component.values())
        role_by_component[index] = min(("calibration", "validation", "test"),
                                       key=lambda role: (counts[role], role))
    for index, people in eligible:
        if index not in role_by_component and any(candidates[p]["role"] == "train" for p in people):
            role_by_component[index] = "train"
    for role in ("calibration", "validation", "test"):
        for source in sources:
            def count():
                return sum(candidates[p]["source"] == source and candidates[p]["track"] == "pure"
                           for index, people in eligible if role_by_component.get(index) == role for p in people)
            available = sorted((index for index, people in eligible if index not in role_by_component
                                and any(candidates[p]["source"] == source for p in people)), key=rank)
            while count() < 3 and available:
                role_by_component[available.pop(0)] = role
            if count() < 3:
                raise ValueError(f"Insufficient separate pure family components for {role}/{source}")
    for index, _ in eligible:
        role_by_component.setdefault(index, "reference")
    return {person: role_by_component[index] for index, people in eligible for person in people}


def split(catalog_dir, screen_dir, output):
    data = core.read(catalog_dir / "catalog.json")
    screen = core.read(screen_dir / "screen.json")
    pairs = core.rows(screen_dir / "related.kin0")
    groups = list(data["pedigree_groups"]) + [[identity] for identity in data["individuals"]]
    for row in pairs:
        if float(row["KINSHIP"]) < screen["threshold"] or int(row["NSNP"]) < 10000:
            raise ValueError("Invalid or underpowered genome-wide kinship pair")
        groups.append([row["IID1"], row["IID2"]])
    grouped = components(groups)
    candidates, sources = data["candidates"], data["template"]["sources"]
    roles = choose_roles(grouped, candidates, data["template"]["seed"], sources)
    samples, excluded = [], []
    family = {}
    for group in grouped:
        component_id = "kin_" + hashlib.sha256(",".join(group).encode()).hexdigest()[:16]
        family.update(dict.fromkeys(group, component_id))
        members = sorted(person for person in group if person in candidates)
        if not members:
            continue
        if len({roles[person] for person in members}) != 1:
            raise ValueError("Family component split across roles")
        # A reference is a donor panel, not another pedigree evaluation role.
        if roles[members[0]] == "reference":
            members.sort(key=lambda person: hashlib.sha256(
                f"{data['template']['seed']}:{person}".encode()).hexdigest())
            for person in members[1:]:
                excluded.append({"sample_id": person, "reason": "one_reference_representative_per_kinship_component",
                                 "retained_representative": members[0], "family_id": component_id})
            members = members[:1]
        for person in members:
            row = dict(candidates[person])
            row.update(role=roles[person], family_id=component_id)
            samples.append({key: row[key] for key in SAMPLE_FIELDS})
    samples.sort(key=lambda row: (core.ROLES.index(row["role"]), row["sample_id"]))
    core.validate_roles(samples, sources)
    for role in core.ROLES:
        tracks = Counter(row["track"] for row in samples if row["role"] == role)
        if role != "reference" and not tracks["mixed"] and role != "train":
            raise ValueError(f"No mixed people in {role}")
    reference_counts = Counter(row["source"] for row in samples if row["role"] == "reference")
    if min(reference_counts[source] for source in sources) < 5:
        raise ValueError("Insufficient unrelated reference representation")
    active = {row["sample_id"]: row for row in samples}
    for identity, node in data["individuals"].items():
        node["kinship_component"] = family[identity]
        component_roles = {roles[p] for p in grouped[next(i for i, group in enumerate(grouped) if identity in group)]
                           if p in roles}
        node["component_role"] = next(iter(component_roles)) if component_roles else None
        node["role"] = active[identity]["role"] if identity in active else None
        node["included_in_pipeline"] = identity in active
        node["genotype_location"] = "genotypes/all/chr{chromosome}/master/all.phased.vcf.gz" if node["genotype_available"] else None
        node["parent_ancestry"] = {
            key: ({"sample_id": node[key], "expected_ancestry": data["individuals"][node[key]]["expected_ancestry"],
                   "realized_ancestry": data["individuals"][node[key]]["realized_ancestry"]}
                  if node.get(key) else None) for key in ("father_id", "mother_id")}
    output.mkdir()
    write_rows(output / "samples.tsv", samples, SAMPLE_FIELDS)
    write_rows(output / "exclusions.tsv", excluded,
               ["sample_id", "reason", "retained_representative", "family_id"])
    core.write(output / "individuals.json", {
        "schema_version": 1, "sources": sources, "fraction_units": "0_to_1",
        "coordinates": "1_based_physical_half_open_intervals; haplotypes 0/1",
        "unknown_values": "null means unknown or unavailable, never zero ancestry",
        "individuals": data["individuals"],
    })
    flat = []
    for identity, node in sorted(data["individuals"].items()):
        row = {key: node.get(key) for key in ("sample_id", "kind", "father_id", "mother_id", "donor_id",
                                             "source", "fine", "role", "component_role", "kinship_component",
                                             "genotype_available", "included_in_pipeline", "ancestry_basis")}
        for kind in ("expected", "realized"):
            for source in sources:
                values = node[f"{kind}_ancestry"]
                row[f"{kind}_{source}"] = values[source] if values is not None else None
        flat.append(row)
    write_rows(output / "individuals.tsv", flat, list(flat[0]))
    core.write(output / "split.json", {
        "roles": roles, "component_members": grouped, "role_counts": dict(Counter(r["role"] for r in samples)),
        "role_tracks": {role: dict(Counter(r["track"] for r in samples if r["role"] == role)) for role in core.ROLES},
        "reference_source_counts": dict(reference_counts), "excluded": excluded, "screen": screen,
        "catalog_sha256": core.digest(catalog_dir / "catalog.json"),
        "selection_uses_model_predictions": False, "independent_test": False,
    })
    print("Kinship-component role counts:", Counter(row["role"] for row in samples), flush=True)


def release(master_dir, split_dir, chrom, output, masks=MASKS):
    samples = core.rows(split_dir / "samples.tsv")
    output.mkdir()
    for kind, roles in (("reference", ["reference"]), ("training", ["train"])):
        people = [row["sample_id"] for row in samples if row["role"] in roles]
        subset_vcf(master_dir / "all.phased.vcf.gz", output / f"{kind}.vcf.gz", people, output / kind)
    query_ids = [row["sample_id"] for row in samples if row["role"] in ("calibration", "validation", "test")]
    for mask in masks:
        sites = core.rows(master_dir / f"{mask}.sites.tsv")
        regions = output / f"{mask}.regions"
        regions.write_text("".join(f"{row['chromosome']}\t{row['position']}\t{row['position']}\n" for row in sites))
        target = output / f"query_{mask}.vcf.gz"
        subset_vcf(master_dir / "all.phased.vcf.gz", target, query_ids, output / mask, regions)
        keys, calls = core.vcf(target, chrom, query_ids, True)
        if [(k[0], str(k[1]), k[2], k[3]) for k in keys] != [
                (r["chromosome"], r["position"], r["ref"], r["alt"]) for r in sites]:
            raise ValueError("Query mask changed its exact allele keys")
        # Use a distinct temporary path: write_vcf refuses accidental indexed-file reuse.
        temporary = output / f"unphased_{mask}.vcf.gz"
        core.write_vcf(temporary, keys, query_ids, calls, phased=False)
        target.unlink()
        Path(str(target) + ".csi").unlink()
        temporary.rename(target)
        Path(str(temporary) + ".csi").rename(str(target) + ".csi")
        regions.unlink()
    core.write(output / "audit.json", {
        "chromosome": chrom, "master_sha256": core.digest(master_dir / "all.phased.vcf.gz"),
        "samples_sha256": core.digest(split_dir / "samples.tsv"),
        "output_sha256": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
    })


def subset_vcf(source, target, people, log, regions=None):
    sample_file = target.with_suffix(".samples")
    sample_file.write_text("".join(person + "\n" for person in people))
    command = ["bcftools", "view", "-S", sample_file, "-Oz", "-o", target]
    if regions:
        command += ["-R", regions]
    core.run(command + [source], log)
    core.run(["bcftools", "index", target], Path(str(log) + ".index"))
    sample_file.unlink()


def finalize(catalog_dir, split_dir, bundles, destination, output):
    data = core.read(catalog_dir / "catalog.json")
    split_data = core.read(split_dir / "split.json")
    if sorted(core.read(path / "audit.json")["chromosome"] for path in bundles) != list(range(1, 23)):
        raise ValueError("Incomplete pipeline release")
    config = copy.deepcopy(data["template"])
    config["seed"] = 20261025
    config["samples"] = str(destination / "metadata/split/samples.tsv")
    config["segments"] = str(destination / "catalog/catalog/segments.tsv")
    for kind in ("reference", "training"):
        config[kind] = {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/{kind}.vcf.gz")
                        for chrom in config["chromosomes"]}
    config["routes"] = [
        {"id": f"harmonized_{mask}", "stage": "harmonized",
         "vcfs": {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/query_{mask}.vcf.gz")
                  for chrom in config["chromosomes"]}} for mask in MASKS]
    config["independent_test"] = False
    core.validate_config(config)
    output.mkdir()
    core.write(output / "sweep.json", config)
    compact = copy.deepcopy(config)
    compact["models"] = [model for model in config["models"] if model["id"] in (
        "flare_full_em_g010", "flare_cap16_em_g010", "gnomix_fast_w050_c2_1",
        "frequency_control", "pca_ridge_control")]
    compact["interpretations"] = [head for head in config["interpretations"] if head["kind"] == "native"]
    if len(compact["models"]) != 5 or len(compact["interpretations"]) != 1:
        raise ValueError("Compact comparison template changed")
    core.validate_config(compact)
    core.write(output / "compact.json", compact)
    smoke = copy.deepcopy(compact)
    smoke["mode"], smoke["chromosomes"] = "smoke", [22]
    for kind in ("reference", "training", "maps", "spans"):
        smoke[kind] = {"22": smoke[kind]["22"]}
    smoke["routes"] = [{**route, "vcfs": {"22": route["vcfs"]["22"]}} for route in smoke["routes"]]
    core.validate_config(smoke)
    core.write(output / "smoke.json", smoke)
    files = {str(path.relative_to(destination)): core.digest(path)
             for path in destination.rglob("*") if path.is_file()}
    files.update({f"configs/{path.name}": core.digest(path) for path in output.iterdir() if path.is_file()})
    core.write(output / "dataset.json", {
        "schema_version": 1, "dataset_id": "kinship_observed_v1_20261025",
        "genome_build": "GRCh38", "sources": config["sources"], "chromosomes": config["chromosomes"],
        "spans": config["spans"], "fraction_units": "0_to_1", "imputation_performed": False,
        "samples": "metadata/split/samples.tsv", "individuals": "metadata/split/individuals.json",
        "individuals_table": "metadata/split/individuals.tsv", "segments": "catalog/catalog/segments.tsv",
        "genotypes": {
            "all_phased": "genotypes/all/chr{chromosome}/master/all.phased.vcf.gz",
            "reference": "genotypes/pipeline/chr{chromosome}/bundle/reference.vcf.gz",
            "training": "genotypes/pipeline/chr{chromosome}/bundle/training.vcf.gz",
            "harmonized": "genotypes/pipeline/chr{chromosome}/bundle/query_{mask}.vcf.gz",
            "masks": list(MASKS), "missing_gt_policy": "forbidden; no filling/DS replacement",
        },
        "configs": {"sweep": "configs/sweep.json", "compact": "configs/compact.json", "smoke": "configs/smoke.json"},
        "role_counts": split_data["role_counts"], "role_tracks": split_data["role_tracks"],
        "kinship_screen": split_data["screen"], "independent_test": False,
        "production_promotion": False, "scientific_certification": False,
        "limitations": [
            "Previously studied development data; rebuilt roles are not a fresh independent test.",
            "KING-robust second-degree threshold is a screening boundary, not proof of zero relatedness.",
            "Admixture, population structure and inbreeding affect relationship estimators.",
            "Published parental information is incomplete outside available pedigrees; unknown is null.",
            "Pure ancestry vectors are operational population labels, not individual biological admixture measurements.",
            "No model fitting or sweep was run while preparing this dataset.",
        ],
        "portable_paths": "Dataset paths are relative to the release root; pipeline configs use local absolute tool/data paths.",
        "files_sha256": files,
    })
    print("Reusable dataset assembled; independent published-output audit is still required", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("catalog", "chromosome", "kinship", "split", "release", "finalize"))
    for name in ("parent", "pedigree", "authority", "catalog", "screen", "split", "master", "dense-root", "destination", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--bundles", nargs="+", type=Path)
    args = parser.parse_args()
    if args.step == "catalog":
        catalog(args.parent, args.pedigree, args.authority, args.output)
    elif args.step == "chromosome":
        chromosome(args.catalog, args.chromosome, args.dense_root, args.output)
    elif args.step == "kinship":
        kinship(args.bundles, args.output)
    elif args.step == "split":
        split(args.catalog, args.screen, args.output)
    elif args.step == "release":
        release(args.master, args.split, args.chromosome, args.output)
    else:
        finalize(args.catalog, args.split, args.bundles, args.destination, args.output)


if __name__ == "__main__":
    main()
