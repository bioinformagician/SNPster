"""Create persistent, disjoint founder trees with pinned Ped-sim meioses."""
import argparse
from collections import Counter, defaultdict
import copy
import hashlib
from pathlib import Path
import shutil

import numpy as np

import core
from prepare_kinship import MASKS, SAMPLE_FIELDS, SEGMENT_FIELDS, ancestry_totals, subset_vcf
from prepare_observed import components, write_rows


def rank(seed, value):
    return hashlib.sha256(f"{seed}:{value}".encode()).hexdigest()


def allocate(nodes, specification, sources):
    if (set(specification) != {"version", "seed", "generations", "minimum_reference_per_source", "founders"}
            or specification["version"] != 1 or specification["generations"] != 3
            or type(specification["seed"]) is not int or not 0 <= specification["seed"] < 2**31
            or type(specification["minimum_reference_per_source"]) is not int
            or specification["minimum_reference_per_source"] < 5
            or set(specification["founders"]) != set(core.ROLES) - {"reference"}):
        raise ValueError("Require a version-1, three-generation, five-role founder design")
    pools = defaultdict(list)
    for identity, node in nodes.items():
        if node["kind"] == "real_person" and node["genotype_available"] and node["source"] in sources:
            pools[node["kinship_component"]].append(identity)
    available, used, selected = defaultdict(list), set(), {}
    for component, people in pools.items():
        labels = {nodes[identity]["source"] for identity in people}
        if len(labels) == 1:
            available[next(iter(labels))].append(component)
    for source in sources:
        available[source].sort(key=lambda value: rank(specification["seed"], value))
    for role in ("train", "calibration", "validation", "test"):
        quota = specification["founders"][role]
        if (set(quota) != set(sources) or any(type(n) is not int or n < 1 for n in quota.values())
                or sum(quota.values()) % 8):
            raise ValueError("Each role needs all sources and a founder count divisible by eight")
        selected[role] = []
        for source in sources:
            choices = [component for component in available[source] if component not in used]
            if len(choices) < quota[source]:
                raise ValueError(f"Not enough independent {source} founders for {role}")
            for component in choices[:quota[source]]:
                used.add(component)
                selected[role].append(min(pools[component],
                                          key=lambda identity: rank(specification["seed"], identity)))
    reference = []
    counts = Counter()
    # Reserve scarce sources first when a legacy connected component contains several labels.
    priority = sorted(sources, key=lambda source: (len(available[source]), source))
    for component in sorted(set(pools) - used):
        person = min(pools[component], key=lambda identity: (
            priority.index(nodes[identity]["source"]), rank(specification["seed"], identity)))
        reference.append(person)
        counts[nodes[person]["source"]] += 1
    if any(counts[source] < specification["minimum_reference_per_source"] for source in sources):
        raise ValueError(f"Insufficient independent reference founders after allocation: {dict(counts)}")
    return selected, reference


def pair_generation(people, nodes, seed):
    remaining = sorted(people, key=lambda identity: rank(seed, identity))
    pairs = []
    while remaining:
        left = remaining.pop(0)
        a = nodes[left]["expected_ancestry"]
        choices = remaining
        if nodes[left]["source"]:
            counts = Counter(nodes[identity]["source"] for identity in remaining)
            choices = [identity for identity in remaining if nodes[identity]["source"] != nodes[left]["source"]
                       and max(count - (source == nodes[identity]["source"]) for source, count in counts.items())
                       <= (len(remaining) - 1) // 2]
            if not choices:
                raise ValueError("Founder design cannot provide different-source first-generation parents")
        right = max(choices, key=lambda identity: (
            sum(abs(a[source] - nodes[identity]["expected_ancestry"][source]) for source in a),
            rank(seed, identity)))
        remaining.remove(right)
        if nodes[left]["source"] and nodes[left]["source"] == nodes[right]["source"]:
            raise ValueError("Founder design cannot provide different-source first-generation parents")
        pairs.append((left, right))
    return pairs


def build_trees(selected, original, sources, seed):
    nodes = copy.deepcopy(original)
    trees = {}
    for role, founders in selected.items():
        people = list(founders)
        for identity in founders:
            nodes[identity].update(role=role, generation=0, founder_donors=[identity])
        for generation in range(1, 4):
            children = []
            for index, (father, mother) in enumerate(pair_generation(people, nodes, seed + generation), 1):
                identity = f"ped_{role}_g{generation}_{index:04d}"
                if identity in nodes:
                    raise ValueError("Synthetic identity collision")
                parents = [nodes[father], nodes[mother]]
                donors = sorted(set(parents[0]["founder_donors"]) | set(parents[1]["founder_donors"]))
                if len(donors) != 2**generation:
                    raise ValueError("Founder reuse or inbreeding in generated tree")
                nodes[identity] = {
                    "sample_id": identity, "kind": "simulated_person", "source": "",
                    "fine": "", "father_id": father, "mother_id": mother, "donor_id": None,
                    "role": role, "generation": generation, "founder_donors": donors,
                    "expected_ancestry": {source: sum(p["expected_ancestry"][source] for p in parents) / 2
                                          for source in sources},
                    "ancestry_basis": "pedigree_mendelian_expectation",
                    "parentage_evidence": "frozen_three_generation_plan",
                }
                children.append(identity)
            people = children
        trees[role] = people
    return nodes, trees


def write_simulation_plan(output, nodes, trees):
    definitions, assignments, aliases = [], [], {}
    for role, roots in trees.items():
        for number, root in enumerate(roots, 1):
            name = f"pd{role}{number:04d}"
            levels = {3: [root]}
            for generation in (2, 1, 0):
                levels[generation] = [nodes[child][key] for child in levels[generation + 1]
                                      for key in ("father_id", "mother_id")]
            for generation in range(4):
                for branch, identity in enumerate(levels[generation], 1):
                    alias = f"{name}1_g{generation + 1}-b{branch}-i1"
                    aliases[alias] = identity
                    if generation == 0:
                        assignments.append(f"{alias}\t{identity}\n")
            definitions.append(f"def {name} 1 4")
            for generation in range(4):
                tokens = []
                for branch, identity in enumerate(levels[generation], 1):
                    sex = "M" if branch % 2 else "F"
                    tokens.append(f"{branch}s{sex}")
                    if generation:
                        tokens.append(f"{branch}:{2 * branch - 1}_{2 * branch}")
                definitions.append(f"{generation + 1} 1 {len(levels[generation])} " + " ".join(tokens))
    (output / "pedigrees.def").write_text("\n".join(definitions) + "\n")
    (output / "set_founders.tsv").write_text("".join(assignments))
    return aliases


def plan(parent, design, output):
    descriptor = core.read(parent / "configs/dataset.json")
    verification = core.read(parent / "verification/independent.json")
    if not verification["verified"] or verification["cross_role_known_or_flagged_relationships"] != 0:
        raise ValueError("Require independently verified kinship-screened source genotypes")
    for relative, expected in descriptor["files_sha256"].items():
        if core.digest(parent / relative) != expected:
            raise ValueError(f"Source dataset changed: {relative}")
    original = core.read(parent / descriptor["individuals"])["individuals"]
    original = {identity: node for identity, node in original.items() if node["kind"] == "real_person"}
    spec = core.read(design)
    selected, reference = allocate(original, spec, descriptor["sources"])
    nodes, trees = build_trees(selected, original, descriptor["sources"], spec["seed"])
    active = set(reference) | {identity for people in selected.values() for identity in people}
    active |= {identity for identity, node in nodes.items() if node["kind"] == "simulated_person"}
    for identity in reference:
        nodes[identity].update(role="reference", generation=0, founder_donors=[identity])
    groups = defaultdict(list)
    for identity, node in original.items():
        groups[node["kinship_component"]].append(identity)
    links = list(groups.values()) + [[identity] for identity in nodes]
    links += [[identity, node["father_id"], node["mother_id"]] for identity, node in nodes.items()
              if node["kind"] == "simulated_person"]
    grouped = components(links)
    for group in grouped:
        family = "pedkin_" + rank(spec["seed"], ",".join(group))[:16]
        roles = {nodes[identity]["role"] for identity in group if identity in active}
        if len(roles) > 1:
            raise ValueError("Biological/simulated connected component crosses roles")
        for identity in group:
            node = nodes[identity]
            node["source_kinship_component"] = node.get("kinship_component")
            node.update(kinship_component=family, component_role=next(iter(roles)) if roles else None,
                        role=node["role"] if identity in active else None,
                        included_in_pipeline=identity in active, genotype_available=identity in active)
    samples = []
    for identity in sorted(active, key=lambda value: (core.ROLES.index(nodes[value]["role"]), value)):
        node = nodes[identity]
        samples.append({
            "sample_id": identity, "role": node["role"], "family_id": node["kinship_component"],
            "donors": ",".join(node["founder_donors"]), "source": node["source"], "fine": node["fine"] or "",
            "track": "pure" if node["generation"] == 0 else f"mixed_g{node['generation']}",
        })
    core.validate_roles(samples, descriptor["sources"])
    output.mkdir()
    write_rows(output / "samples.tsv", samples, SAMPLE_FIELDS)
    aliases = write_simulation_plan(output, nodes, trees)
    core.write(output / "plan.json", {
        "specification": spec, "parent": str(parent), "parent_descriptor_sha256": core.digest(parent / "configs/dataset.json"),
        "sources": descriptor["sources"], "spans": descriptor["spans"], "nodes": nodes, "aliases": aliases,
        "selected_founders": selected, "reference": reference,
        "role_generation_counts": {role: dict(Counter(str(nodes[r["sample_id"]]["generation"])
                                                        for r in samples if r["role"] == role))
                                   for role in core.ROLES},
        "terminal_trees": trees, "independent_test": False,
        "simulation": "Ped-sim --pois; sex-averaged GRCh38 map; zero genotyping error/missingness; phased founder GT",
    })
    print("Founder-disjoint role/generation counts:", core.read(output / "plan.json")["role_generation_counts"], flush=True)


def bp_origins(path, chrom, aliases):
    from realized_ancestry import load_bp
    parsed = load_bp(path, chrom)
    if set(parsed) != set(aliases):
        raise ValueError("Ped-sim breakpoint sample axis changed")
    return {aliases[identity]: {int(hap[-1]): [(a, b + 1, hid) for a, b, hid in segments]
                               for hap, segments in homologs.items()}
            for identity, homologs in parsed.items()}


def simulate(plan_dir, chrom, output, source_master=None):
    from convert_plink_map import convert
    from verify_founder_assignments import run as verify_assignments
    data = core.read(plan_dir / "plan.json")
    parent = Path(data["parent"])
    config = data["template"] if "template" in data else core.read(parent / "configs/compact.json")
    origin = source_master if source_master is not None else parent / f"genotypes/all/chr{chrom}/master"
    output.mkdir()
    scratch = output.parent / "simulation"
    scratch.mkdir()
    convert(config["maps"][str(chrom)], scratch / "maps")
    founders = sorted(identity for people in data["selected_founders"].values() for identity in people)
    subset_vcf(origin / "all.phased.vcf.gz", scratch / "founders.vcf.gz", founders, scratch / "founders")
    seed = (data["specification"]["seed"] + chrom * 15427) % 2000000000
    core.run(["ped-sim", "-d", plan_dir / "pedigrees.def", "-m", scratch / f"maps/chr{chrom}.simmap",
              "-i", scratch / "founders.vcf.gz", "-o", scratch / "sim",
              "--set_founders", plan_dir / "set_founders.tsv", "--founder_ids", "--bp_all",
              "--keep_phase", "--nogz", "--err_rate", "0", "--miss_rate", "0", "--seed", seed, "--pois"],
             scratch / "pedsim")
    verify_assignments(plan_dir / "set_founders.tsv", scratch / "sim.ids")
    origins = bp_origins(scratch / "sim.bp", chrom, data["aliases"])
    haplotypes = {}
    for identity in founders:
        for hap in (0, 1):
            segments = origins[identity][hap]
            ids = {hid for _, _, hid in segments}
            if len(ids) != 1:
                raise ValueError("Founder haplotype changes identity")
            hid = ids.pop()
            if hid in haplotypes:
                raise ValueError("Founder haplotype ID reused")
            haplotypes[hid] = (identity, hap)
    identities = [row["sample_id"] for row in core.rows(plan_dir / "samples.tsv")]
    real = [identity for identity in identities if data["nodes"][identity]["kind"] == "real_person"]
    keys, real_calls = core.vcf(origin / "all.phased.vcf.gz", chrom, real, True)
    simulated_ids = [alias for alias, identity in data["aliases"].items()
                     if data["nodes"][identity]["kind"] == "simulated_person"]
    simkeys, simcalls = core.vcf(scratch / "sim.vcf", chrom, simulated_ids, True)
    if keys != simkeys:
        raise ValueError("Ped-sim changed exact observed allele keys")
    simindex = {data["aliases"][alias]: index for index, alias in enumerate(simulated_ids)}
    realindex = {identity: index for index, identity in enumerate(real)}
    calls = np.empty((len(identities), 2, len(keys)), dtype=np.int8)
    for index, identity in enumerate(identities):
        calls[index] = real_calls[realindex[identity]] if identity in realindex else simcalls[simindex[identity]]
    core.write_vcf(output / "all.phased.vcf.gz", keys, identities, calls)
    segments, inheritance = [], []
    start, end = data["spans"][str(chrom)]
    for identity in identities:
        node = data["nodes"][identity]
        for hap in (0, 1):
            if node["kind"] == "real_person":
                pieces = [(start, end, identity, hap)]
            else:
                pieces = [(max(start, a), min(end, b), *haplotypes[hid])
                          for a, b, hid in origins[identity][hap] if max(start, a) < min(end, b)]
            cursor = start
            for a, b, donor, donor_hap in pieces:
                if a != cursor or b <= a or donor not in node["founder_donors"]:
                    raise ValueError("Invalid founder-origin tiling")
                cursor = b
                segments.append({"sample_id": identity, "chromosome": str(chrom), "haplotype": str(hap),
                                 "start": a, "end": b, "source": data["nodes"][donor]["source"]})
                inheritance.append({**segments[-1], "donor_id": donor, "donor_haplotype": donor_hap})
            if cursor != end:
                raise ValueError("Ped-sim breakpoint coverage does not include complete benchmark span")
    write_rows(output / "segments.tsv", segments, SEGMENT_FIELDS)
    write_rows(output / "inheritance.tsv", inheritance, SEGMENT_FIELDS + ["donor_id", "donor_haplotype"])
    for mask in data.get("masks", MASKS):
        shutil.copyfile(origin / f"{mask}.sites.tsv", output / f"{mask}.sites.tsv")
    for suffix in ("bp", "ids"):
        shutil.copyfile(scratch / f"sim.{suffix}", output / f"pedsim.{suffix}")
    shutil.copyfile(Path(str(scratch / "pedsim") + ".log"), output / "pedsim.log")
    core.write(output / "audit.json", {
        "chromosome": chrom, "seed": seed, "plan_sha256": core.digest(plan_dir / "plan.json"),
        "source_vcf_sha256": core.digest(origin / "all.phased.vcf.gz"),
        "map_sha256": core.digest(config["maps"][str(chrom)]),
        "people": len(identities), "markers": len(keys),
        "output_sha256": {path.name: core.digest(path) for path in output.iterdir() if path.is_file()},
    })


def finalize(plan_dir, bundles, destination, output):
    data = core.read(plan_dir / "plan.json")
    ordered = sorted(bundles, key=lambda path: core.read(path / "audit.json")["chromosome"])
    if [core.read(path / "audit.json")["chromosome"] for path in ordered] != list(range(1, 23)):
        raise ValueError("All 22 autosomes are required")
    output.mkdir()
    metadata = output / "metadata"
    metadata.mkdir()
    segments = [row for bundle in ordered for row in core.rows(bundle / "segments.tsv")]
    write_rows(metadata / "segments.tsv", segments, SEGMENT_FIELDS)
    samples = core.rows(plan_dir / "samples.tsv")
    shutil.copyfile(plan_dir / "samples.tsv", metadata / "samples.tsv")
    values = ancestry_totals(segments, data["sources"], data["spans"])
    nodes = data["nodes"]
    for identity, node in nodes.items():
        node["realized_ancestry"] = values.get(identity)
    for node in nodes.values():
        node["parent_ancestry"] = {
            key: ({"sample_id": node[key], "expected_ancestry": nodes[node[key]]["expected_ancestry"],
                   "realized_ancestry": nodes[node[key]]["realized_ancestry"]}
                  if node.get(key) and node[key] in nodes else None) for key in ("father_id", "mother_id")}
    core.write(metadata / "individuals.json", {
        "schema_version": 1, "sources": data["sources"], "fraction_units": "0_to_1",
        "coordinates": "1_based_physical_half_open_intervals; haplotypes 0/1",
        "individuals": nodes,
    })
    flat = []
    for identity, node in sorted(nodes.items()):
        row = {key: node.get(key) for key in (
            "sample_id", "role", "generation", "kind", "source", "fine", "father_id", "mother_id",
            "kinship_component", "source_kinship_component", "genotype_available")}
        row["founder_donors"] = ",".join(node.get("founder_donors", []))
        for kind in ("expected", "realized"):
            for source in data["sources"]:
                vector = node.get(f"{kind}_ancestry")
                row[f"{kind}_{source}"] = vector[source] if vector else None
        flat.append(row)
    write_rows(metadata / "individuals.tsv", flat, list(flat[0]))
    config = core.read(Path(data["parent"]) / "configs/sweep.json")
    config.update(seed=data["specification"]["seed"], independent_test=False,
                  samples=str(destination / "metadata/samples.tsv"), segments=str(destination / "metadata/segments.tsv"))
    for kind in ("reference", "training"):
        config[kind] = {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/{kind}.vcf.gz")
                        for chrom in range(1, 23)}
    config["routes"] = [
        {"id": f"harmonized_{mask}", "stage": "harmonized",
         "vcfs": {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/query_{mask}.vcf.gz")
                  for chrom in range(1, 23)}} for mask in MASKS]
    configs = output / "configs"
    configs.mkdir()
    core.validate_config(config)
    core.write(configs / "sweep.json", config)
    compact = copy.deepcopy(config)
    template = core.read(Path(data["parent"]) / "configs/compact.json")
    compact.update(models=template["models"], interpretations=template["interpretations"])
    core.validate_config(compact)
    core.write(configs / "compact.json", compact)
    smoke = copy.deepcopy(compact)
    smoke.update(mode="smoke", chromosomes=[22])
    for key in ("reference", "training", "maps", "spans"):
        smoke[key] = {"22": smoke[key]["22"]}
    smoke["routes"] = [{**route, "vcfs": {"22": route["vcfs"]["22"]}} for route in smoke["routes"]]
    core.validate_config(smoke)
    core.write(configs / "smoke.json", smoke)
    files = {str(path.relative_to(destination)): core.digest(path)
             for path in destination.rglob("*") if path.is_file()}
    files.update({str(path.relative_to(output)): core.digest(path)
                  for path in output.rglob("*") if path.is_file()})
    core.write(configs / "dataset.json", {
        "schema_version": 1, "dataset_id": destination.parent.name, "genome_build": "GRCh38",
        "sources": data["sources"], "spans": data["spans"], "chromosomes": list(range(1, 23)),
        "samples": "metadata/samples.tsv", "individuals": "metadata/individuals.json",
        "individuals_table": "metadata/individuals.tsv", "segments": "metadata/segments.tsv",
        "role_generation_counts": data["role_generation_counts"],
        "role_counts": dict(Counter(row["role"] for row in samples)),
        "genotypes": {
            "all_phased": "genotypes/all/chr{chromosome}/master/all.phased.vcf.gz",
            "reference": "genotypes/pipeline/chr{chromosome}/bundle/reference.vcf.gz",
            "training": "genotypes/pipeline/chr{chromosome}/bundle/training.vcf.gz",
            "harmonized": "genotypes/pipeline/chr{chromosome}/bundle/query_{mask}.vcf.gz",
            "inheritance": "genotypes/all/chr{chromosome}/master/inheritance.tsv",
        },
        "configs": {kind: f"configs/{kind}.json" for kind in ("sweep", "compact", "smoke")},
        "simulation": data["simulation"], "parent": data["parent"],
        "parent_descriptor_sha256": data["parent_descriptor_sha256"],
        "independent_test": False, "imputation_performed": False, "scientific_certification": False,
        "limitations": [
            "Reused development donors: not an untouched external test.",
            "Pure labels are population labels, not measured individual ancestry.",
            "Generations are pedigree depth, not a guaranteed number of ancestry sources.",
            "Three-generation halving trees have 8 founders per terminal family; relatives are correlated.",
            "This design does not certify low-fraction ancestry below nominal 1/8.",
            "Oceanian and Native American founder diversity is limited.",
            "KING screening is not proof of absolute unrelatedness; inherited components are conservatively retained.",
            "Phased simulated union SNP VCFs are observed-panel truth, not dense whole-genome or validated imputation.",
        ],
        "files_sha256": files,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("plan", "simulate", "finalize"))
    for name in ("parent", "design", "plan", "destination", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--bundles", nargs="+", type=Path)
    args = parser.parse_args()
    if args.step == "plan":
        plan(args.parent, args.design, args.output)
    elif args.step == "simulate":
        simulate(args.plan, args.chromosome, args.output)
    else:
        finalize(args.plan, args.bundles, args.destination, args.output)


if __name__ == "__main__":
    main()
