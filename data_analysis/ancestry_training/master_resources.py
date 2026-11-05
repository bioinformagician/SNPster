"""Consolidate every local reference resource without inventing observed calls or ancestry."""
import argparse
from collections import Counter, defaultdict
import copy
import csv
import gzip
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

import core
from prepare_kinship import SAMPLE_FIELDS, SEGMENT_FIELDS, ancestry_totals, collapse, onehot, release
from prepare_observed import components, write_rows
from prepare_pedigree import build_trees, rank, simulate, write_simulation_plan


RESOURCE_IDS = ("hgdp1kg", "gse77508", "gse73996", "kazakh", "sgdp")
SITE_FIELDS = ["chromosome", "position", "ref", "alt"]
PRIMARY_ROLES = ("train", "calibration", "validation", "test")


def canonical(identity):
    result = re.sub(r"\s+", "_", identity.strip())
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", result):
        raise ValueError(f"Unsupported source identity: {identity!r}")
    return result


def identities(path):
    values = subprocess.check_output(["bcftools", "query", "-l", str(path)], text=True).splitlines()
    if not values or len(values) != len(set(values)):
        raise ValueError(f"Empty or duplicate sample axis: {path}")
    return values


def sites(path, chrom, regions=None):
    command = ["bcftools", "query"]
    command += ["-R", str(regions)] if regions else ["-r", str(chrom)]
    command += ["-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(path)]
    result = []
    for line in subprocess.check_output(command, text=True).splitlines():
        chromosome, pos, ref, alt = line.split("\t")
        if chromosome == str(chrom) and len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and ref != alt:
            result.append((chromosome, int(pos), ref, alt))
    if len(result) != len(set(result)):
        raise ValueError(f"Duplicate exact source variants: {path}")
    return result


def read_gt(path, chrom, people, regions=None):
    command = ["bcftools", "query", "-s", ",".join(people)]
    command += ["-R", str(regions)] if regions else ["-r", str(chrom)]
    command += ["-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)]
    keys, values = [], []
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t", 4)
            if (fields[0] != str(chrom) or len(fields[2]) != 1 or len(fields[3]) != 1
                    or fields[2] not in "ACGT" or fields[3] not in "ACGT" or fields[2] == fields[3]):
                continue
            try:
                alleles = core.decode_gt(fields[4], len(people), allow_missing=True)
            except ValueError as error:
                raise ValueError(f"{error}: {path}:{fields[1]}") from error
            keys.append((fields[0], int(fields[1]), fields[2], fields[3]))
            values.append(alleles)
        if process.wait():
            raise RuntimeError(f"Source GT query failed: {path}")
    if not keys:
        raise ValueError(f"No qualified source variants: {path}")
    return keys, np.asarray(values, dtype=np.int8).transpose(1, 2, 0)


def write_gt(path, keys, people, calls, phased=False):
    if (calls.shape != (len(people), 2, len(keys)) or not calls.size or calls.dtype.kind not in "biu"
            or calls.min() < -1 or calls.max() > 1):
        raise ValueError("Invalid archive GT axes or alleles")
    raw = path.with_suffix("")
    sep = "|" if phased else "/"
    coding = {-1: ".", 0: "0", 1: "1"}
    encoded = np.asarray([f"{coding[a]}{sep}{coding[b]}\t".encode("ascii")
                          for a in (-1, 0, 1) for b in (-1, 0, 1)], dtype="S4")
    with raw.open("wb") as handle:
        handle.write((f'##fileformat=VCFv4.2\n##contig=<ID={keys[0][0]}>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Source calls; missing is not filled">\n'
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(people) + "\n").encode("utf-8"))
        for offset, (chrom, pos, ref, alt) in enumerate(keys):
            values = encoded[(calls[:, 0, offset] + 1) * 3 + calls[:, 1, offset] + 1].tobytes()
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t".encode("ascii"))
            handle.write(values[:-1] + b"\n")
    core.run(["bcftools", "view", "-Oz", "-o", path, raw], path.with_suffix(".compress"))
    core.run(["bcftools", "index", path], path.with_suffix(".index"))
    raw.unlink()


def source_files(spec, chrom):
    root, project = Path(spec["reference_root"]), Path(spec["project"])
    qc = core.read(project / spec["public_qc"] / "preparation/plan.json")
    return {
        "hgdp1kg": [root / f"hgdp_1kg_references/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"],
        "gse77508": [root / f"gse77508/GSE77508.chr{chrom}.phased.vcf.gz"],
        "gse73996": [root / f"harmonized/gse73996/selected/gse73996.chr{chrom}.grch38.vcf.gz"],
        "kazakh": [root / f"harmonized/kazakh/vcf/kazakh.chr{chrom}.grch38.vcf.gz"],
        "sgdp": [Path(qc["jobs"][str(chrom)]["old"]), Path(qc["jobs"][str(chrom)]["new"])],
    }


def inventory(specification, output):
    specification = specification.resolve()
    spec = core.read(specification)
    sources = spec["ancestry_axis"]
    if spec["version"] != 1 or len(set(sources)) != 7 or set(spec["panels"]) != {"observed", "gse77508", "gse73996", "sgdp"}:
        raise ValueError("Unsupported full-resource specification")
    root, project = Path(spec["reference_root"]), Path(spec["project"])
    authority = root / "flare2_test/source_model/ancestry_role_audit.tsv"
    authority_rows = core.rows(authority)
    labels = {row["sample"]: row for row in authority_rows}
    if len(labels) != len(authority_rows):
        raise ValueError("Duplicate or conflicting authoritative sample labels")
    old = core.read(Path(spec["observed_dataset"]) / "metadata/split/individuals.json")["individuals"]
    public = core.read(project / spec["public_qc"] / "preparation/plan.json")
    if set(public["sources"]) != set(public["public_ids"]):
        raise ValueError("Public observed donor source labels are incomplete")
    template = core.read(Path(spec["template_dataset"]) / "configs/sweep.json")
    compact = core.read(Path(spec["template_dataset"]) / "configs/compact.json")
    design_path = specification.resolve().parent / spec["founder_design"]
    design = core.read(design_path)
    resources, nodes = {}, {}
    files = {str(path): core.digest(path) for path in (
        specification, design_path, authority, Path(spec["public_pedigree"]),
        project / spec["public_qc"] / "preparation/plan.json",
        project / spec["public_qc"] / "summary/genetic_screen.json",
        Path(spec["template_dataset"]) / "configs/sweep.json",
        Path(spec["template_dataset"]) / "configs/compact.json")}
    for filename in set(template["maps"].values()) | {template["tools"]["beagle_jar"]}:
        files[str(Path(filename).resolve())] = core.digest(filename)
    for resource, paths in source_files(spec, 22).items():
        mappings, people = [], []
        for number, path in enumerate(paths):
            raw = identities(path)
            mapped = {identity: canonical(identity) for identity in raw}
            if len(set(mapped.values())) != len(mapped):
                raise ValueError(f"Identity sanitization collision: {resource}")
            mappings.append(mapped)
            selected = [identity for identity in raw
                        if not (resource == "sgdp" and number == 1 and identity == "S_Chane-1.DG")]
            for identity in selected:
                person = mapped[identity]
                if person in nodes:
                    raise ValueError(f"Unresolved duplicate resource identity: {person}")
                if resource == "sgdp":
                    source, fine, eligibility, affinity = collapse(public["sources"][person]), "SGDP_" + person.rsplit("-", 1)[0], "source", ""
                else:
                    if person not in labels:
                        raise ValueError(f"Resource person has no reviewed label: {person}")
                    label = labels[person]
                    source, fine, eligibility, affinity = collapse(label["source_group"]), label["fine"], label["role"], label["affinity_group"]
                if eligibility == "source" and source not in sources:
                    raise ValueError("Source label outside configured ancestry axis")
                prior = old.get(person, {})
                nodes[person] = {
                    "sample_id": person, "raw_sample_id": identity, "resource": resource, "source_file_index": number,
                    "kind": "real_person", "source": source if eligibility == "source" else "",
                    "fine": fine, "label_role": eligibility, "affinity_group": affinity,
                    "father_id": prior.get("father_id"), "mother_id": prior.get("mother_id"), "donor_id": None,
                    "public_family_id": prior.get("public_family_id"),
                    "ancestry_basis": "operational_population_label_not_measured_individual_admixture" if eligibility == "source" else "unknown_additive_ancestry",
                    "expected_ancestry": onehot(source, sources) if eligibility == "source" else None,
                    "genotype_available": True, "generation": 0, "founder_donors": [person],
                    "in_fixed_production_imputation_reference": person in set(public["fixed_ids"]),
                }
                people.append(person)
        resources[resource] = {"people": people, "raw_to_canonical": mappings,
                               "count": len(people), "input_files": {
                                   str(chrom): [str(path) for path in source_files(spec, chrom)[resource]]
                                   for chrom in range(1, 23)}}
    if set(labels) != {identity for identity, node in nodes.items() if node["resource"] != "sgdp"}:
        raise ValueError("Full legacy reference inventory was not covered exactly")
    public_families = defaultdict(list)
    for row in core.rows(spec["public_pedigree"]):
        person = row["Individual ID"]
        if person in nodes:
            nodes[person].update(father_id=None if row["Paternal ID"] in ("0", "") else row["Paternal ID"],
                                 mother_id=None if row["Maternal ID"] in ("0", "") else row["Maternal ID"],
                                 public_family_id=row["Family ID"])
            public_families[row["Family ID"]].append(person)
    output.mkdir()
    core.write(output / "inventory.json", {
        "specification": spec, "resources": resources, "nodes": nodes, "template": template, "compact": compact,
        "founder_design": design, "public_family_groups": list(public_families.values()), "input_sha256": files,
        "public_source_expected_sha256": {str(chrom): public["jobs"][str(chrom)]["input_sha256"]
                                           for chrom in range(1, 23)},
        "coordinate_policy": "GRCh38; numeric autosomes; half-open physical ancestry intervals",
        "missingness_policy": "Missing and unobserved GT stay missing; no DS/GP substitution or genotype filling",
        "gse77508_provenance": "Provided phased VCF; original per-call pre-phase missingness is not independently certified",
        "full_source_policy": "Original full-resource VCF paths and SHA256 are retained; derived archive views preserve target-panel GT and missingness",
    })
    print("COMPLETE RESOURCE CENSUS:", {r: value["count"] for r, value in resources.items()},
          "total sample identities", len(nodes), flush=True)


def project_chromosome(inventory_dir, chrom, output):
    data = core.read(inventory_dir / "inventory.json")
    spec, resources = data["specification"], data["resources"]
    output.mkdir()
    scratch = output.parent / "projection_work"
    scratch.mkdir()
    backbone = Path(spec["observed_dataset"]) / f"genotypes/all/chr{chrom}/master/all.phased.vcf.gz"
    desired = set(sites(backbone, chrom))
    for resource in ("gse77508", "gse73996", "kazakh"):
        desired.update(sites(Path(resources[resource]["input_files"][str(chrom)][0]), chrom))
    regions = scratch / "positions.tsv"
    regions.write_text("".join(f"{chrom}\t{pos}\t{pos}\n" for pos in sorted({key[1] for key in desired})))
    canonical_keys = [key for key in sites(Path(resources["hgdp1kg"]["input_files"][str(chrom)][0]), chrom, regions)
                      if key in desired]
    counts = Counter(key[1] for key in canonical_keys)
    canonical_keys = [key for key in canonical_keys if counts[key[1]] == 1]
    allowed = set(canonical_keys)
    print(f"chr{chrom}: {len(canonical_keys)} exact canonical model-union SNPs; source GT projection beginning", flush=True)
    regions.write_text("".join(f"{chrom}\t{key[1]}\t{key[1]}\n" for key in canonical_keys))
    blocks, reports, all_ids = {}, {}, []
    for resource in RESOURCE_IDS:
        info = resources[resource]
        ids = info["people"]
        calls = np.full((len(ids), 2, len(canonical_keys)), -1, dtype=np.int8)
        row_index = {identity: offset for offset, identity in enumerate(ids)}
        column_index = {key: offset for offset, key in enumerate(canonical_keys)}
        present = set()
        hashes = {}
        for number, path_string in enumerate(info["input_files"][str(chrom)]):
            path = Path(path_string)
            raw = identities(path)
            mapping = info["raw_to_canonical"][number]
            if set(raw) != set(mapping):
                raise ValueError(f"Source sample axis differs on chr{chrom}: {resource}")
            selected = [identity for identity in raw if mapping[identity] in row_index
                        and not (resource == "sgdp" and number == 1 and identity == "S_Chane-1.DG")]
            keys, values = read_gt(path, chrom, selected, regions)
            selected_offsets = [(index, column_index[key]) for index, key in enumerate(keys) if key in allowed]
            for old_column, new_column in selected_offsets:
                calls[[row_index[mapping[identity]] for identity in selected], :, new_column] = values[:, :, old_column]
                present.add(canonical_keys[new_column])
            hashes[path_string] = core.digest(path)
        offsets = [index for index, key in enumerate(canonical_keys) if key in present]
        keys = [canonical_keys[index] for index in offsets]
        projected = calls[:, :, offsets]
        write_gt(output / f"{resource}.vcf.gz", keys, ids, projected, phased=resource in ("hgdp1kg", "gse77508"))
        blocks[resource] = (ids, calls, present)
        all_ids += ids
        reports[resource] = {
            "people": len(ids), "variants": len(keys), "called": {
                identity: int((projected[index] >= 0).all(axis=0).sum()) for index, identity in enumerate(ids)},
            "source_sha256": hashes, "projection_sha256": core.digest(output / f"{resource}.vcf.gz"),
        }
        print(f"chr{chrom}/{resource}: archived {len(ids)} people at {len(keys)} SNPs without GT filling", flush=True)
    common = set.intersection(*(blocks[resource][2] for resource in RESOURCE_IDS[:4]))
    common_offsets = [index for index, key in enumerate(canonical_keys) if key in common]
    common_keys = [canonical_keys[index] for index in common_offsets]
    common_calls = np.concatenate([blocks[resource][1][:, :, common_offsets] for resource in RESOURCE_IDS])
    blocks.clear()
    write_gt(output / "identity.vcf.gz", common_keys, all_ids, common_calls)
    del common_calls
    plink = ["plink2", "--threads", "1", "--memory", "1536"]
    core.run(plink + ["--vcf", output / "identity.vcf.gz", "--double-id", "--set-all-var-ids", "@:#:$r:$a",
                      "--maf", ".05", "--make-pgen", "--out", output / "all"], output / "plink-import")
    core.run(plink + ["--pfile", output / "all", "--indep-pairwise", "200kb", "1", ".2",
                      "--out", output / "ld"], output / "plink-prune")
    core.run(plink + ["--pfile", output / "all", "--extract", output / "ld.prune.in",
                      "--make-pgen", "--out", output / "pruned"], output / "plink-export")
    for suffix in (".pgen", ".pvar", ".psam"):
        (output / "all").with_suffix(suffix).unlink()
    for mask in ("lonni", "frederik"):
        shutil.copyfile(Path(spec["observed_dataset"]) / f"genotypes/all/chr{chrom}/master/{mask}.sites.tsv",
                        output / f"{mask}.sites.tsv")
    core.write(output / "projection.json", {
        "chromosome": chrom, "resources": reports, "people": all_ids,
        "identity_variants_before_ld": len(common_keys),
        "identity_variants_after_ld": len((output / "ld.prune.in").read_text().splitlines()),
        "ambiguous_canonical_positions_rejected": sum(n > 1 for n in counts.values()),
        "identity_only_intersection_not_model_backbone": True,
    })
    print(f"chr{chrom}: retained all {len(all_ids)} resource identities; {len(common_keys)} common identity SNPs; native model panels preserved", flush=True)


def screen_relatedness(bundles, inventory_dir, output):
    data = core.read(inventory_dir / "inventory.json")
    spec = data["specification"]
    ordered = sorted(bundles, key=lambda path: core.read(path / "projection.json")["chromosome"])
    if [core.read(path / "projection.json")["chromosome"] for path in ordered] != list(range(1, 23)):
        raise ValueError("All autosomes required for full-resource relatedness")
    axes = [core.read(path / "projection.json")["people"] for path in ordered]
    if any(axis != axes[0] for axis in axes):
        raise ValueError("Resource identities vary between autosomes")
    output.mkdir()
    (output / "merge.list").write_text("".join(str(path.resolve() / "pruned") + "\n" for path in ordered))
    plink = ["plink2", "--threads", "4", "--memory", "28672"]
    core.run(plink + ["--pmerge-list", output / "merge.list", "--make-pgen", "--out", output / "genome"], output / "merge")
    core.run(plink + ["--pfile", output / "genome", "--make-king-table", "counts",
                      "--out", output / "pairs"], output / "king")
    related, near_duplicates, minimum, checked, seen_ids = [], [], None, 0, set()
    group_min = {}
    with (output / "pairs.kin0").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t", skipinitialspace=True)
        for raw in reader:
            row = {key.strip(): value.strip() for key, value in raw.items()}
            left, right, n = row["IID1"], row["IID2"], int(row["NSNP"])
            if left == right or left not in data["nodes"] or right not in data["nodes"]:
                raise ValueError("Invalid KING identity pair")
            checked += 1
            seen_ids.update((left, right))
            minimum = n if minimum is None else min(minimum, n)
            group = "/".join(sorted((data["nodes"][left]["resource"], data["nodes"][right]["resource"])))
            group_min[group] = min(group_min.get(group, n), n)
            value = float(row["KINSHIP"])
            if not np.isfinite(value):
                raise ValueError(f"Nonfinite KING estimate: {left}/{right}")
            if value >= spec["kinship_threshold"]:
                related.append({"first": left, "second": right, "called_markers": n, "kinship": value})
            if value >= spec["near_duplicate_threshold"]:
                near_duplicates.append({"first": left, "second": right, "called_markers": n, "kinship": value})
    expected = len(axes[0]) * (len(axes[0]) - 1) // 2
    if checked != expected or seen_ids != set(axes[0]) or minimum < spec["minimum_pair_markers"]:
        raise ValueError(f"Insufficient all-pair identity evidence: {checked}/{expected} pairs; minimum called SNPs {minimum}; per-resource {group_min}")
    with (output / "pairs.kin0").open("rb") as source, gzip.open(output / "pairs.kin0.gz", "wb", compresslevel=1) as target:
        shutil.copyfileobj(source, target)
    (output / "pairs.kin0").unlink()
    core.write(output / "screen.json", {
        "people": axes[0], "pairs_checked": checked, "minimum_pair_called_markers": minimum,
        "minimum_pair_markers_by_resource": group_min, "related_pairs": related,
        "near_duplicate_or_mz_pairs": near_duplicates, "kinship_threshold": spec["kinship_threshold"],
        "near_duplicate_threshold": spec["near_duplicate_threshold"],
        "simulation_edges_used": False, "pair_table_sha256": core.digest(output / "pairs.kin0.gz"),
        "scope": "biological pedigrees plus all-resource genotype screening; not proof of absolute unrelatedness",
    })
    print(f"Full-resource KING: {checked:,} pairs, minimum {minimum:,} SNPs per pair, {len(related)} flagged edges", flush=True)


def import_preparation(prepared, inventory_dir, output):
    current = core.read(inventory_dir / "inventory.json")
    original = core.read(prepared / "inventory/inventory.json")
    if current["nodes"] != original["nodes"] or current["resources"] != original["resources"]:
        raise ValueError("Prepared source people, labels or genotype paths differ")
    for key in ("reference_root", "project", "observed_dataset", "template_dataset", "public_qc", "ancestry_axis"):
        if current["specification"][key] != original["specification"][key]:
            raise ValueError(f"Prepared source protocol differs: {key}")
    for filename, expected in original["input_sha256"].items():
        if not Path(filename).is_absolute() or core.digest(filename) != expected:
            raise ValueError(f"Prepared scientific input changed: {filename}")
    screen = core.read(prepared / "kinship/screen.json")
    n = len(current["nodes"])
    if (screen["pairs_checked"] != n * (n - 1) // 2
            or set(screen["people"]) != set(current["nodes"])
            or screen["minimum_pair_called_markers"] < current["specification"]["minimum_pair_markers"]
            or screen["kinship_threshold"] != current["specification"]["kinship_threshold"]
            or screen["near_duplicate_threshold"] != current["specification"]["near_duplicate_threshold"]
            or core.digest(prepared / "kinship/pairs.kin0.gz") != screen["pair_table_sha256"]):
        raise ValueError("Prepared all-pair genetic evidence is incomplete, changed or incompatible")
    original_sources = {}
    fingerprints = {}
    for chrom in range(1, 23):
        folder = prepared / f"archive/chr{chrom}/projection"
        report = core.read(folder / "projection.json")
        audit_path = prepared / f"verification/resources/chr{chrom}.json"
        audit = core.read(audit_path)
        if not audit["verified"] or audit["chromosome"] != chrom or report["chromosome"] != chrom:
            raise ValueError("Prepared all-autosome original-source GT audits are incomplete")
        for resource, item in report["resources"].items():
            if core.digest(folder / f"{resource}.vcf.gz") != item["projection_sha256"]:
                raise ValueError("An independently audited source archive changed")
            original_sources.update(item["source_sha256"])
        fingerprints[str(audit_path)] = core.digest(audit_path)
    for filename, expected in original_sources.items():
        if core.digest(filename) != expected:
            raise ValueError(f"Original resource GT changed since its independent audit: {filename}")
    producer = prepared.parents[1]
    subprocess.run(["sha256sum", "--check", str(producer / "source.sha256")],
                   cwd=producer / "source", check=True, stdout=subprocess.DEVNULL)
    output.mkdir()
    shutil.copytree(prepared / "archive", output / "archive")
    shutil.copytree(prepared / "kinship", output / "kinship")
    (output / "verification").mkdir()
    shutil.copytree(prepared / "verification/resources", output / "verification/resources")
    core.write(output / "verification/preparation_import.json", {
        "verified": True, "prepared_source": str(prepared), "source_people": n,
        "original_source_sha256": original_sources, "independent_source_audit_sha256": fingerprints,
        "prepared_producer_manifest_sha256": core.digest(producer / "source.sha256"),
        "original_workflow_completion_not_assumed": True,
        "scope": "Only audited native resource genotypes and complete all-pair screening are reused; roles and descendants are rebuilt",
    })
    print(f"VERIFIED PREPARATION IMPORT: all {n} resource identities, 22 original-GT audits and {screen['pairs_checked']:,} screened pairs", flush=True)


def assign_roles(nodes, groups, design, seed):
    ownership, selected = {}, {role: [] for role in PRIMARY_ROLES}
    index_by_person = {identity: index for index, group in enumerate(groups) for identity in group}
    for role in PRIMARY_ROLES:
        for source, needed in design["founders"][role].items():
            choices = []
            for index, group in enumerate(groups):
                candidates = [identity for identity in group if nodes[identity]["resource"] == "hgdp1kg"
                              and nodes[identity]["label_role"] == "source" and nodes[identity]["source"] == source]
                if index not in ownership and candidates:
                    choices.append((index, candidates))
            choices.sort(key=lambda item: (-len(groups[item[0]]), rank(seed, ",".join(groups[item[0]]))))
            if len(choices) < needed:
                raise ValueError(f"Insufficient biological founder families for {role}/{source}")
            for index, candidates in choices[:needed]:
                ownership[index] = role
                selected[role].append(min(candidates, key=lambda identity: rank(seed, identity)))
    # Non-HGDP cohorts must contribute real held-out people, not only reference rows.
    cycle = ("validation", "test", "calibration", "train", "reference")
    strata = defaultdict(list)
    for index, group in enumerate(groups):
        if index in ownership:
            continue
        source_people = [identity for identity in group if nodes[identity]["label_role"] == "source"]
        if source_people:
            representative = min(source_people, key=lambda identity: (
                nodes[identity]["resource"] == "hgdp1kg", rank(seed, identity)))
            node = nodes[representative]
            if node["resource"] != "hgdp1kg":
                strata[(node["resource"], node["source"])].append(index)
    for stratum, indices in sorted(strata.items()):
        indices.sort(key=lambda index: rank(seed, ",".join(groups[index])))
        for offset, index in enumerate(indices):
            ownership[index] = cycle[offset % len(cycle)]
    for index, group in enumerate(groups):
        ownership.setdefault(index, "reference" if any(nodes[p]["label_role"] == "source" for p in group) else "affinity_only")
    for group in groups:
        index = index_by_person[group[0]]
        role = ownership[index]
        source_people = [person for person in group if nodes[person]["label_role"] == "source"]
        representative = min(source_people, key=lambda person: (
            RESOURCE_IDS.index(nodes[person]["resource"]), rank(seed, person))) if source_people else None
        for person in group:
            node = nodes[person]
            node["component_role"] = role
            if node["label_role"] != "source":
                node["role"], node["reason"] = "affinity_only", "No measured additive ancestry truth; population-affinity use only"
            elif role == "reference" and person != representative:
                node["role"], node["reason"] = "related_reference_stress", "Biological relative/alias of a retained fitting reference; excluded from primary accuracy, retained for stress tests"
            else:
                node["role"], node["reason"] = role, "Whole biological family assigned without model predictions"
            node["kinship_component"] = "bio_" + rank(seed, ",".join(group))[:16]
            node["source_kinship_component"] = node["kinship_component"]
    return selected, ownership


def make_cohorts(inventory_dir, screen_dir, destination, output):
    data = core.read(inventory_dir / "inventory.json")
    nodes = copy.deepcopy(data["nodes"])
    screen = core.read(screen_dir / "screen.json")
    links = [[identity] for identity in nodes] + data["public_family_groups"]
    links += [[row["first"], row["second"]] for row in screen["related_pairs"]]
    grouped = components(links)
    selected, _ = assign_roles(nodes, grouped, data["founder_design"], data["specification"]["seed"])
    reference_counts = Counter(node["source"] for node in nodes.values() if node["role"] == "reference"
                               and node["resource"] == "hgdp1kg")
    if any(reference_counts[source] < data["founder_design"]["minimum_reference_per_source"]
           for source in data["specification"]["ancestry_axis"]):
        raise ValueError(f"Insufficient independent HGDP/1KG reference families: {dict(reference_counts)}")
    # Very high KING is not a new independent person; retain alternate records outside primary scoring.
    duplicates = components([[row["first"], row["second"]] for row in screen["near_duplicate_or_mz_pairs"]])
    seed_ids = {identity for people in selected.values() for identity in people}
    for group in duplicates:
        primary = sorted(group, key=lambda person: (person not in seed_ids, RESOURCE_IDS.index(nodes[person]["resource"]),
                                                    rank(data["specification"]["seed"], person)))
        for person in primary[1:]:
            if person in seed_ids:
                raise ValueError("Near-duplicate/MZ seed founders must not be counted twice")
            if nodes[person]["role"] in PRIMARY_ROLES:
                nodes[person]["role"] = "identity_replicate_stress"
                nodes[person]["reason"] = "Near-duplicate or monozygotic-twin genotype record; same family role, not another independent primary person"
    nodes, trees = build_trees(selected, nodes, data["specification"]["ancestry_axis"], data["specification"]["seed"])
    for identity, node in nodes.items():
        if node["kind"] == "simulated_person":
            parent_components = {nodes[donor]["kinship_component"] for donor in node["founder_donors"]}
            node["biological_founder_components"] = sorted(parent_components)
    # New simulations connect donor families only within their already frozen role.
    new_groups = components(grouped + [[identity, node["father_id"], node["mother_id"]] for identity, node in nodes.items()
                                       if node["kind"] == "simulated_person"])
    for group in new_groups:
        family = "pedbio_" + rank(data["specification"]["seed"], ",".join(group))[:16]
        for person in group:
            nodes[person]["evaluation_family"] = family
    output.mkdir()
    core.write(output / "individuals.json", {
        "schema_version": 1, "sources": data["specification"]["ancestry_axis"], "individuals": nodes,
        "biological_components_ignore_old_simulations": True,
        "fraction_units": "0_to_1; null means unknown, not zero",
    })
    fields = ["sample_id", "raw_sample_id", "resource", "label_role", "source", "fine", "affinity_group",
              "role", "component_role", "kinship_component", "evaluation_family", "father_id", "mother_id",
              "generation", "in_fixed_production_imputation_reference", "reason"]
    write_rows(output / "people.tsv", [{key: node.get(key) for key in fields} for node in nodes.values()
                                       if node["kind"] == "real_person"], fields)
    for panel in data["specification"]["panels"]:
        folder = output / panel
        folder.mkdir()
        included_resources = {"hgdp1kg"} | ({panel} if panel != "observed" else set())
        active = {identity for identity, node in nodes.items() if node["kind"] == "simulated_person"
                  or (node["resource"] in included_resources and node["role"] in core.ROLES)}
        panel_nodes = copy.deepcopy(nodes)
        for identity, node in panel_nodes.items():
            node["included_in_pipeline"] = identity in active
            node["genotype_available"] = identity in active
            if node["kind"] == "simulated_person":
                node["component_role"] = node["role"]
            if identity in active:
                node["kinship_component"] = node["evaluation_family"]
        samples = [{
            "sample_id": identity, "role": node["role"], "family_id": node["evaluation_family"],
            "donors": ",".join(node["founder_donors"]), "source": node["source"], "fine": node["fine"] or "",
            "track": "pure" if node["generation"] == 0 else f"mixed_g{node['generation']}",
        } for identity, node in panel_nodes.items() if identity in active]
        samples.sort(key=lambda row: (core.ROLES.index(row["role"]), row["sample_id"]))
        core.validate_roles(samples, data["specification"]["ancestry_axis"])
        masks = ["lonni", "frederik"] if panel == "observed" else [panel]
        aliases = write_simulation_plan(folder, panel_nodes, trees)
        write_rows(folder / "samples.tsv", samples, SAMPLE_FIELDS)
        core.write(folder / "plan.json", {
            "specification": data["founder_design"], "sources": data["specification"]["ancestry_axis"],
            "spans": data["template"]["spans"], "nodes": panel_nodes, "aliases": aliases, "masks": masks,
            "selected_founders": selected, "reference": [r["sample_id"] for r in samples if r["role"] == "reference"],
            "parent": str(destination / f"panels/{panel}/input"), "template": data["template"],
            "terminal_trees": trees, "independent_test": False, "panel": panel,
            "role_generation_counts": {role: dict(Counter(str(panel_nodes[r["sample_id"]]["generation"])
                                                            for r in samples if r["role"] == role))
                                       for role in core.ROLES},
        })
    core.write(output / "split.json", {
        "role_counts": dict(Counter(node["role"] for node in nodes.values() if node["kind"] == "real_person")),
        "selected_founder_counts": {role: len(people) for role, people in selected.items()},
        "biological_components": len(grouped), "new_evaluation_families": len(new_groups),
        "near_duplicate_or_mz_groups": duplicates, "old_simulation_edges_used_for_biological_split": False,
        "selection_uses_model_predictions": False,
    })
    print("All-resource biological-only split:", core.read(output / "split.json")["role_counts"], flush=True)


def base_panel(projection, cohort, inventory_dir, chrom, output):
    data = core.read(inventory_dir / "inventory.json")
    plan = core.read(cohort / "plan.json")
    panel = plan["panel"]
    samples = core.rows(cohort / "samples.tsv")
    real = [row for row in samples if plan["nodes"][row["sample_id"]]["kind"] == "real_person"]
    resources = ["hgdp1kg"] + ([] if panel == "observed" else [panel])
    blocks = {}
    for resource in resources:
        people = [row["sample_id"] for row in real if plan["nodes"][row["sample_id"]]["resource"] == resource]
        if not people:
            raise ValueError(f"No active primary people from {resource}/{panel}")
        blocks[resource] = (*read_gt(projection / f"{resource}.vcf.gz", chrom, people), people)
    common = set.intersection(*(set(keys) for keys, _, _ in blocks.values()))
    if panel == "observed":
        common &= {(row["chromosome"], int(row["position"]), row["ref"], row["alt"])
                   for mask in ("lonni", "frederik") for row in core.rows(projection / f"{mask}.sites.tsv")}
    start, end = plan["spans"][str(chrom)]
    keys = sorted(key for key in common if start <= key[1] < end)
    initially_shared = len(keys)
    complete = np.ones(len(keys), dtype=bool)
    aligned = {}
    for resource, (original, values, people) in blocks.items():
        index = {key: offset for offset, key in enumerate(original)}
        aligned[resource] = values[:, :, [index[key] for key in keys]]
        complete &= (aligned[resource] >= 0).all(axis=(0, 1))
    keys = [key for key, keep in zip(keys, complete) if keep]
    if len(keys) < data["specification"]["minimum_pipeline_markers"]:
        raise ValueError(f"{panel}/chr{chrom}: only {len(keys)} jointly callable model markers; source archive retained unchanged")
    for resource in aligned:
        aligned[resource] = aligned[resource][:, :, complete]
    people = [row["sample_id"] for row in real]
    calls = np.empty((len(people), 2, len(keys)), dtype=np.int8)
    positions = {identity: offset for offset, identity in enumerate(people)}
    output.mkdir()
    scratch = output.parent / "phase_only"
    scratch.mkdir()
    for resource, (_, _, identities_in_resource) in blocks.items():
        for offset, identity in enumerate(identities_in_resource):
            calls[positions[identity]] = aligned[resource][offset]
    phase_checks = []
    if panel != "observed":
        reference_ids = [row["sample_id"] for row in real if row["role"] == "reference"
                         and plan["nodes"][row["sample_id"]]["resource"] == "hgdp1kg"]
        if not reference_ids:
            raise ValueError("Phase-only adapter requires the disjoint HGDP fitting reference")
        core.write_vcf(scratch / "phase_reference.vcf.gz", keys, reference_ids,
                       calls[[positions[identity] for identity in reference_ids]])
        for role in core.ROLES:
            ids = [row["sample_id"] for row in real if row["role"] == role
                   and plan["nodes"][row["sample_id"]]["resource"] != "hgdp1kg"]
            if not ids:
                continue
            original = calls[[positions[identity] for identity in ids]].copy()
            target = scratch / f"{role}.vcf.gz"
            core.write_vcf(target, keys, ids, original, phased=False)
            prefix = scratch / f"{role}.phase"
            core.run(["java", core.java_heap(), "-jar", plan["template"]["tools"]["beagle_jar"],
                      f"gt={target}", f"ref={scratch / 'phase_reference.vcf.gz'}",
                      f"map={plan['template']['maps'][str(chrom)]}", f"out={prefix}",
                      "impute=false", "nthreads=1", f"seed={data['specification']['seed']}"],
                     scratch / f"{role}.invocation")
            actual_keys, phased = core.vcf(str(prefix) + ".vcf.gz", chrom, ids, True)
            if actual_keys != keys or not np.array_equal(phased.sum(axis=1), original.sum(axis=1)):
                raise ValueError("Phase-only processing changed a typed GT or added a locus")
            calls[[positions[identity] for identity in ids]] = phased
            phase_checks.append({"role": role, "people": len(ids), "changed_diploid_calls": 0,
                                 "reference_contains_calibration_validation_test": False})
    core.write_vcf(output / "all.phased.vcf.gz", keys, people, calls)
    retained_keys = set(keys)
    for mask in plan["masks"]:
        selected = (core.rows(projection / f"{mask}.sites.tsv") if panel == "observed"
                    else [dict(zip(SITE_FIELDS, key)) for key in keys])
        retained = [row for row in selected if
                    (str(row["chromosome"]), int(row["position"]), row["ref"], row["alt"]) in retained_keys]
        if len(retained) < data["specification"]["minimum_pipeline_markers"]:
            raise ValueError(f"{panel}/{mask}/chr{chrom}: insufficient retained route markers")
        write_rows(output / f"{mask}.sites.tsv", retained, SITE_FIELDS)
    core.write(output / "base.json", {
        "panel": panel, "chromosome": chrom, "people": people, "variants": len(keys),
        "initially_shared_variants": initially_shared, "derived_complete_panel_sites_removed": int((~complete).sum()),
        "archive_genotypes_or_missingness_changed": False, "genotypes_filled": 0,
        "phase_only_by_role": phase_checks, "source_projection_sha256": {
            resource: core.digest(projection / f"{resource}.vcf.gz") for resource in resources},
        "output_sha256": core.digest(output / "all.phased.vcf.gz"),
        "policy": "Complete-case model view only; unchanged source GT/missingness remains in the resource archive",
    })
    print(f"{panel}/chr{chrom}: {len(people)} real people including pedigree founders; {len(keys)} model markers; no GT filling", flush=True)


def simulate_panel(cohort, base, chrom, output):
    simulate(cohort, chrom, output, source_master=base)


def export_panel(cohort, master, chrom, output):
    plan = core.read(cohort / "plan.json")
    release(master, cohort, chrom, output, masks=plan["masks"])


def unscored_cases(projection, cohorts, inventory_dir, chrom, output):
    inventory_data = core.read(inventory_dir / "inventory.json")
    nodes = core.read(cohorts / "individuals.json")["individuals"]
    output.mkdir()
    records = {}
    for resource in RESOURCE_IDS:
        ids = inventory_data["resources"][resource]["people"]
        for track, accepted in (("affinity", {"affinity_only"}),
                                ("relatedness_stress", {"related_reference_stress", "identity_replicate_stress"})):
            selected = [identity for identity in ids if nodes[identity]["role"] in accepted]
            name = f"{resource}.{track}"
            records[name] = {"people": len(selected), "primary_accuracy_eligible": False}
            if selected:
                from prepare_kinship import subset_vcf
                subset_vcf(projection / f"{resource}.vcf.gz", output / f"{name}.vcf.gz", selected, output / name)
    core.write(output / "cases.json", {"chromosome": chrom, "tracks": records,
                                      "unknown_ancestry_is_not_zero": True})


def finalize_cohort(cohort, bundles, destination, output):
    plan = core.read(cohort / "plan.json")
    ordered = sorted(bundles, key=lambda path: core.read(path / "audit.json")["chromosome"])
    if [core.read(path / "audit.json")["chromosome"] for path in ordered] != list(range(1, 23)):
        raise ValueError("Every platform cohort requires all 22 autosomes")
    output.mkdir()
    metadata, configs = output / "metadata", output / "configs"
    metadata.mkdir()
    configs.mkdir()
    samples = core.rows(cohort / "samples.tsv")
    shutil.copyfile(cohort / "samples.tsv", metadata / "samples.tsv")
    segments = [row for bundle in ordered for row in core.rows(bundle / "segments.tsv")]
    write_rows(metadata / "segments.tsv", segments, SEGMENT_FIELDS)
    values = ancestry_totals(segments, plan["sources"], plan["spans"])
    nodes = copy.deepcopy(plan["nodes"])
    for identity, node in nodes.items():
        node["realized_ancestry"] = values.get(identity)
    for node in nodes.values():
        node["parent_ancestry"] = {
            key: ({"sample_id": node[key], "expected_ancestry": nodes[node[key]]["expected_ancestry"],
                   "realized_ancestry": nodes[node[key]]["realized_ancestry"]}
                  if node.get(key) in nodes else None) for key in ("father_id", "mother_id")}
    core.write(metadata / "individuals.json", {"schema_version": 1, "panel": plan["panel"], "sources": plan["sources"], "individuals": nodes,
                                              "fraction_units": "0_to_1; null means unknown"})
    config = copy.deepcopy(plan["template"])
    config.update(seed=plan["specification"]["seed"], independent_test=False,
                  samples=str(destination / "metadata/samples.tsv"), segments=str(destination / "metadata/segments.tsv"))
    for kind in ("reference", "training"):
        config[kind] = {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/{kind}.vcf.gz")
                        for chrom in range(1, 23)}
    config["routes"] = [
        {"id": f"harmonized_{mask}", "stage": "harmonized",
         "vcfs": {str(chrom): str(destination / f"genotypes/pipeline/chr{chrom}/bundle/query_{mask}.vcf.gz")
                  for chrom in range(1, 23)}} for mask in plan["masks"]]
    core.validate_config(config)
    core.write(configs / "sweep.json", config)
    compact = copy.deepcopy(config)
    compact["models"] = [model for model in config["models"] if model["id"] in (
        "flare_full_em_g010", "flare_cap16_em_g010", "gnomix_fast_w050_c2_1", "frequency_control", "pca_ridge_control")]
    compact["interpretations"] = [head for head in config["interpretations"] if head["kind"] == "native"]
    if len(compact["models"]) != 5 or len(compact["interpretations"]) != 1:
        raise ValueError("Compact comparison template differs")
    core.validate_config(compact)
    core.write(configs / "compact.json", compact)
    smoke = copy.deepcopy(compact)
    smoke.update(mode="smoke", chromosomes=[22])
    for key in ("reference", "training", "maps", "spans"):
        smoke[key] = {"22": smoke[key]["22"]}
    smoke["routes"] = [{**route, "vcfs": {"22": route["vcfs"]["22"]}} for route in smoke["routes"]]
    core.validate_config(smoke)
    core.write(configs / "smoke.json", smoke)
    core.write(configs / "dataset.json", {
        "schema_version": 1, "panel": plan["panel"], "sources": plan["sources"], "spans": plan["spans"],
        "samples": "metadata/samples.tsv", "segments": "metadata/segments.tsv",
        "individuals": "metadata/individuals.json", "chromosomes": list(range(1, 23)),
        "role_counts": dict(Counter(row["role"] for row in samples)),
        "role_generation_counts": plan["role_generation_counts"], "independent_test": False,
        "imputation_performed": False, "models_fitted": False,
        "model_routes": plan["masks"],
        "scope": "Assay-specific model view of a full-resource biological-family split, not a universal marker intersection",
    })


def finalize_master(inventory_dir, cohorts, projections, releases, config_releases, destination, output):
    inventory_data = core.read(inventory_dir / "inventory.json")
    spec = inventory_data["specification"]
    nodes = core.read(cohorts / "individuals.json")["individuals"]
    real = {identity: node for identity, node in nodes.items() if node["kind"] == "real_person"}
    if set(real) != set(inventory_data["nodes"]):
        raise ValueError("A full-resource identity was dropped")
    cohort_data = {core.read(path / "dataset.json")["panel"]: path for path in config_releases}
    cohort_nodes = {core.read(path / "individuals.json")["panel"]: core.read(path / "individuals.json")["individuals"]
                    for path in releases}
    if set(cohort_data) != set(spec["panels"]) or set(cohort_nodes) != set(spec["panels"]):
        raise ValueError("Missing or duplicate platform release")
    panels = {panel: core.read(path / "dataset.json") for panel, path in cohort_data.items()}
    for identity, node in nodes.items():
        if node["kind"] == "simulated_person":
            vectors = [cohort_nodes[panel][identity]["realized_ancestry"] for panel in spec["panels"]]
            if any(not np.allclose(list(vector.values()), list(vectors[0].values()), rtol=0, atol=1e-12)
                   for vector in vectors[1:]):
                raise ValueError("The same simulated person has different ancestry across assay panels")
            node["realized_ancestry"] = vectors[0]
        else:
            node["realized_ancestry"] = node["expected_ancestry"]
    for node in nodes.values():
        node["parent_ancestry"] = {
            key: ({"sample_id": node[key], "expected_ancestry": nodes[node[key]]["expected_ancestry"],
                   "realized_ancestry": nodes[node[key]]["realized_ancestry"]}
                  if node.get(key) in nodes else None) for key in ("father_id", "mother_id")}
    output.mkdir()
    core.write(output / "individuals.json", {
        "schema_version": 1, "sources": spec["ancestry_axis"], "individuals": nodes,
        "fraction_units": "0_to_1; null means unknown", "unknown_affinity_ancestry": "null, never one-hot or zero",
        "pure_truth": "Operational population labels; not measured biological admixture",
    })
    rows = []
    for node in nodes.values():
        row = {key: node.get(key) for key in ("sample_id", "resource", "kind", "role", "component_role", "source",
                                              "fine", "generation", "father_id", "mother_id", "evaluation_family",
                                              "kinship_component", "in_fixed_production_imputation_reference", "reason")}
        row["founder_donors"] = ",".join(node["founder_donors"])
        for kind in ("expected", "realized"):
            for source in spec["ancestry_axis"]:
                vector = node.get(f"{kind}_ancestry")
                row[f"{kind}_{source}"] = vector[source] if vector is not None else None
        rows.append(row)
    write_rows(output / "individuals.tsv", rows, list(rows[0]))
    source_hashes = {}
    for bundle in projections:
        report = core.read(bundle / "projection.json")
        for item in report["resources"].values():
            source_hashes.update(item["source_sha256"])
    core.write(output / "dataset.json", {
        "schema_version": 1, "dataset_id": destination.parent.parent.name, "genome_build": "GRCh38",
        "sample_identity_records": len(real), "resources": {resource: {
            "people": value["count"], "original_full_vcfs": value["input_files"],
            "target_panel_archive": f"archive/chr{{chromosome}}/projection/{resource}.vcf.gz",
            "roles": dict(Counter(real[identity]["role"] for identity in value["people"])),
        } for resource, value in inventory_data["resources"].items()},
        "model_cohorts": {panel: {"descriptor": f"panels/{panel}/configs/dataset.json",
                                 "compact": f"panels/{panel}/configs/compact.json",
                                 "sweep": f"panels/{panel}/configs/sweep.json",
                                 "role_counts": value["role_counts"]} for panel, value in panels.items()},
        "individuals": "master/individuals.json", "individuals_table": "master/individuals.tsv",
        "biological_only_split": True, "old_simulation_links_used_for_biological_split": False,
        "reference_relatives_kept_as_unscored_stress": True, "affinity_people_have_known_ancestry_truth": False,
        "source_sha256": source_hashes, "models_fitted": False,
        "limitations": [
            "Reused donors are not an untouched independent external test.",
            "KING screening and high-kinship/MZ flags do not prove absolute unrelatedness or exact identity.",
            "Related people within held-out families are correlated; use family-aware uncertainty.",
            "Cross-platform marker sets differ; do not interpret them as a matched-marker causal experiment.",
            "Complete-case model views remove missing sites only in that derived cohort; resource archives keep missingness unchanged.",
            "GSE77508 supplied phased calls lack independently certified original per-call missingness provenance.",
            "Affinity populations have unknown additive fractions and cannot supply percentage-accuracy truth.",
            "Production imputation-reference membership is retained; these observed-stage results cannot certify independent imputed-stage accuracy.",
            "External full-source VCFs remain at recorded immutable locations; target-panel archives and model cohorts are copied into this release.",
        ],
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("inventory", "project", "kinship", "import-preparation", "split", "base", "simulate",
                                       "export", "cases", "finalize-cohort", "finalize-master"))
    for name in ("specification", "inventory", "cohorts", "cohort", "screen", "projection", "base",
                 "master", "destination", "output", "prepared"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--bundles", nargs="+", type=Path)
    parser.add_argument("--releases", nargs="+", type=Path)
    parser.add_argument("--config-releases", nargs="+", type=Path)
    args = parser.parse_args()
    if args.step == "inventory":
        inventory(args.specification, args.output)
    elif args.step == "project":
        project_chromosome(args.inventory, args.chromosome, args.output)
    elif args.step == "kinship":
        screen_relatedness(args.bundles, args.inventory, args.output)
    elif args.step == "import-preparation":
        import_preparation(args.prepared, args.inventory, args.output)
    elif args.step == "split":
        make_cohorts(args.inventory, args.screen, args.destination, args.output)
    elif args.step == "base":
        base_panel(args.projection, args.cohort, args.inventory, args.chromosome, args.output)
    elif args.step == "simulate":
        simulate_panel(args.cohort, args.base, args.chromosome, args.output)
    elif args.step == "export":
        export_panel(args.cohort, args.master, args.chromosome, args.output)
    elif args.step == "cases":
        unscored_cases(args.projection, args.cohorts, args.inventory, args.chromosome, args.output)
    elif args.step == "finalize-cohort":
        finalize_cohort(args.cohort, args.bundles, args.destination, args.output)
    else:
        finalize_master(args.inventory, args.cohorts, args.bundles, args.releases, args.config_releases, args.destination, args.output)


if __name__ == "__main__":
    main()
