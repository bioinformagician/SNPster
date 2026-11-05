"""Independent full-resource genotype, biological-family and published-interface checks."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import subprocess
import tempfile

import numpy as np

import core
from audit_pedigree import chromosome as audit_pedigree_chromosome
from prepare_kinship import ancestry_totals


def text_calls(path, regions=None, people=None):
    command = ["bcftools", "query"]
    if regions:
        command += ["-R", str(regions)]
    if people:
        command += ["-s", ",".join(people)]
    command += ["-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)]
    result = {}
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[2], fields[3])
            if key in result:
                raise ValueError("Duplicate exact audit variant")
            result[key] = np.asarray(["./." if value == "." else value.replace("|", "/") for value in fields[4:]], dtype="S3")
        if process.wait():
            raise RuntimeError(f"Independent source audit query failed: {path}")
    return result


def project(projection, inventory, output):
    data = core.read(inventory / "inventory.json")
    report = core.read(projection / "projection.json")
    chrom = report["chromosome"]
    checks = {}
    with tempfile.TemporaryDirectory(prefix="snpster-resource-audit-") as temporary:
        regions = Path(temporary) / "sites.tsv"
        for resource, info in data["resources"].items():
            path = projection / f"{resource}.vcf.gz"
            if core.digest(path) != report["resources"][resource]["projection_sha256"]:
                raise ValueError("Published archive differs from producer fingerprint")
            actual = text_calls(path)
            regions.write_text("".join(f"{chrom}\t{key[1]}\t{key[1]}\n" for key in actual))
            observed_ids = subprocess.check_output(["bcftools", "query", "-l", str(path)], text=True).splitlines()
            if observed_ids != info["people"]:
                raise ValueError("A resource identity was lost, reordered or duplicated")
            offsets = {identity: index for index, identity in enumerate(observed_ids)}
            checked, absent_missing = 0, 0
            expected_hashes = data.get("public_source_expected_sha256", {}).get(str(chrom), {})
            for number, filename in enumerate(info["input_files"][str(chrom)]):
                if core.digest(filename) != report["resources"][resource]["source_sha256"][filename]:
                    raise ValueError("Original source GT changed during projection")
                if filename in expected_hashes and core.digest(filename) != expected_hashes[filename]:
                    raise ValueError("Original SGDP calls differ from independently normalized source")
                mapping = info["raw_to_canonical"][number]
                raw = [identity for identity in mapping if mapping[identity] in offsets
                       and not (resource == "sgdp" and number == 1 and identity == "S_Chane-1.DG")]
                expected = text_calls(filename, regions, raw)
                positions = [offsets[mapping[identity]] for identity in raw]
                for key, values in actual.items():
                    cells = values[positions]
                    if key in expected:
                        if not np.array_equal(cells, expected[key]):
                            raise ValueError(f"Archive altered observed GT/missingness: {resource}/{key}")
                        checked += len(cells)
                    else:
                        if not np.all(cells == b"./."):
                            raise ValueError("An unobserved source marker was assigned nonmissing GT")
                        absent_missing += len(cells)
            checks[resource] = {"people": len(observed_ids), "observed_diploid_cells_checked": checked,
                                "absent_record_cells_preserved_as_missing": absent_missing}
    core.write(output, {"verified": True, "chromosome": chrom, "resources": checks,
                        "missing_genotypes_filled": 0, "original_source_genotypes_changed": 0})


def base(base_dir, projection, cohort, chrom, output):
    plan = core.read(cohort / "plan.json")
    rows = core.rows(cohort / "samples.tsv")
    people = [row["sample_id"] for row in rows if plan["nodes"][row["sample_id"]]["kind"] == "real_person"]
    keys, phased = core.vcf(base_dir / "all.phased.vcf.gz", chrom, people, True)
    by_resource = defaultdict(list)
    for index, identity in enumerate(people):
        by_resource[plan["nodes"][identity]["resource"]].append((index, identity))
    regions = output.parent / f"base-audit-{plan['panel']}-{chrom}.regions"
    regions.write_text("".join(f"{chrom}\t{key[1]}\t{key[1]}\n" for key in keys))
    try:
        for resource, members in by_resource.items():
            original = text_calls(projection / f"{resource}.vcf.gz", regions, [identity for _, identity in members])
            for column, key in enumerate(keys):
                for offset, (index, identity) in enumerate(members):
                    text = original[key][offset].decode("ascii")
                    if "." in text:
                        raise ValueError("A source missing GT entered the complete model panel")
                    alleles = [int(value) for value in text.split("/")]
                    if int(phased[index, :, column].sum()) != sum(alleles):
                        raise ValueError(f"Phase-only panel changed a source dosage: {identity}/{key}")
                    if resource == "hgdp1kg" and list(phased[index, :, column]) != alleles:
                        raise ValueError("Provided HGDP founder homologs changed")
    finally:
        regions.unlink()
    core.write(output, {"verified": True, "panel": plan["panel"], "chromosome": chrom,
                        "source_diploid_cells_verified": len(people) * len(keys),
                        "source_missing_gt_filled": 0, "phase_only_changed_diploid_calls": 0})


def verify(dataset, output):
    descriptor = core.read(dataset / "master/dataset.json")
    inventory = core.read(dataset / "inventory/inventory.json")
    nodes = core.read(dataset / descriptor["individuals"])["individuals"]
    real = {identity: node for identity, node in nodes.items() if node["kind"] == "real_person"}
    if set(real) != set(inventory["nodes"]) or len(real) != descriptor["sample_identity_records"]:
        raise ValueError("The master inventory does not cover every local resource identity")
    previous = core.read(output) if output.exists() else None
    if previous:
        for relative, expected in previous["files_sha256"].items():
            if core.digest(dataset / relative) != expected:
                raise ValueError(f"A verified durable file changed: {relative}")
    for filename, expected in {**inventory["input_sha256"], **descriptor["source_sha256"]}.items():
        if core.digest(filename) != expected:
            raise ValueError(f"A scientific source changed: {filename}")
    owners = defaultdict(set)
    for identity, node in real.items():
        owners[node["kinship_component"]].add(node["component_role"])
        if node["role"] == "affinity_only" and (node["expected_ancestry"] is not None or node["realized_ancestry"] is not None):
            raise ValueError("Unknown affinity ancestry was converted into invented truth")
    if any(len(roles) != 1 for roles in owners.values()):
        raise ValueError("A biological family crosses primary split ownership")
    screen = core.read(dataset / "kinship/screen.json")
    n = len(real)
    if screen["pairs_checked"] != n * (n - 1) // 2 or screen["minimum_pair_called_markers"] < inventory["specification"]["minimum_pair_markers"]:
        raise ValueError("All-pair identity evidence is incomplete or underpowered")
    if core.digest(dataset / "kinship/pairs.kin0.gz") != screen["pair_table_sha256"]:
        raise ValueError("Durable all-pair KING evidence changed")
    for edge in screen["related_pairs"]:
        if real[edge["first"]]["component_role"] != real[edge["second"]]["component_role"]:
            raise ValueError("An actual flagged biological relationship crosses primary roles")
    for group in inventory["public_family_groups"]:
        if len({real[person]["component_role"] for person in group}) != 1:
            raise ValueError("Published biological family split across roles")
    for chrom in range(1, 23):
        archive = dataset / f"archive/chr{chrom}/projection"
        report = core.read(archive / "projection.json")
        audit = core.read(dataset / f"verification/resources/chr{chrom}.json")
        if not audit["verified"] or audit["chromosome"] != chrom:
            raise ValueError("Missing source-genotype audit")
        for resource, record in report["resources"].items():
            if core.digest(archive / f"{resource}.vcf.gz") != record["projection_sha256"]:
                raise ValueError("Durable resource archive changed")
    interfaces, alleles = [], 0
    for panel in descriptor["model_cohorts"]:
        root = dataset / f"panels/{panel}"
        plan = core.read(dataset / f"cohorts/{panel}/plan.json")
        rows = core.rows(root / "metadata/samples.tsv")
        core.validate_roles(rows, inventory["specification"]["ancestry_axis"])
        truth = ancestry_totals(core.rows(root / "metadata/segments.tsv"), plan["sources"], plan["spans"])
        for row in rows:
            identity = row["sample_id"]
            node = nodes[identity]
            if row["role"] != node["role"] or node["role"] not in core.ROLES:
                raise ValueError("A stress/affinity record entered primary model scoring")
            np.testing.assert_allclose(list(truth[identity].values()),
                                       [node["realized_ancestry"][source] for source in truth[identity]], rtol=0, atol=1e-12)
            if node["kind"] == "simulated_person":
                parents = [nodes[node[key]] for key in ("father_id", "mother_id")]
                if any(parent["role"] != node["role"] or parent["generation"] != node["generation"] - 1 for parent in parents):
                    raise ValueError("Invalid simulated pedigree role/generation")
                if set(parents[0]["founder_donors"]) & set(parents[1]["founder_donors"]):
                    raise ValueError("Founder reuse/inbreeding in a simulated tree")
        for chrom in range(1, 23):
            master = root / f"genotypes/all/chr{chrom}/master"
            base_dir = root / f"input/chr{chrom}/base"
            base_report = core.read(base_dir / "base.json")
            if core.digest(base_dir / "all.phased.vcf.gz") != base_report["output_sha256"]:
                raise ValueError("Durable phase-only source panel changed")
            source_audit = core.read(root / f"verification/chromosomes/chr{chrom}.source.json")
            if not source_audit["verified"] or source_audit["source_missing_gt_filled"]:
                raise ValueError("Missing independent phase-only/unchanged-source-GT audit")
            producer = core.read(master / "audit.json")
            for filename, expected in producer["output_sha256"].items():
                if core.digest(master / filename) != expected:
                    raise ValueError("A durable pedigree output changed")
            audit = core.read(root / f"verification/chromosomes/chr{chrom}.json")
            if not audit["verified"] or not audit["mendelian_alleles_verified"]:
                raise ValueError("Missing exact independent founder-haplotype audit")
            alleles += audit["founder_haplotype_alleles_independently_checked"]
        with tempfile.TemporaryDirectory(prefix="snpster-full-resource-interface-") as temporary:
            root_tmp = Path(temporary)
            core.freeze(root / "configs/smoke.json", root_tmp / "frozen")
            for mask in plan["masks"]:
                core.prepare(root_tmp / "frozen", f"harmonized_{mask}", 22, root_tmp / mask)
                interfaces.append(f"{panel}/harmonized_{mask}/chr22")
    excluded = {str(output.relative_to(dataset))}
    files = {str(path.relative_to(dataset)): core.digest(path) for path in dataset.rglob("*")
             if path.is_file() and str(path.relative_to(dataset)) not in excluded}
    if previous and set(files) != set(previous["files_sha256"]):
        raise ValueError("A durable file was added or removed after verification")
    output.parent.mkdir(exist_ok=True, parents=True)
    core.write(output, {"verified": True, "sample_identity_records": n, "models_fitted": False,
                        "independent_test": False, "all_pair_relationships_checked": screen["pairs_checked"],
                        "minimum_pair_called_markers": screen["minimum_pair_called_markers"],
                        "founder_haplotype_alleles_independently_checked": alleles,
                        "actual_pipeline_interfaces_passed": interfaces,
                        "resources": descriptor["resources"], "files_sha256": files})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("project", "panel", "verify"))
    for name in ("projection", "inventory", "base", "master", "bundle", "cohort", "dataset", "output"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--chromosome", type=int)
    args = parser.parse_args()
    if args.step == "project":
        project(args.projection, args.inventory, args.output)
    elif args.step == "panel":
        source_audit = args.output.with_suffix(".source.json")
        base(args.base, args.projection, args.cohort, args.chromosome, source_audit)
        audit_pedigree_chromosome(args.master, args.bundle, args.cohort, args.chromosome,
                                 args.output, source_master=args.base)
    else:
        verify(args.dataset, args.output)


if __name__ == "__main__":
    main()
