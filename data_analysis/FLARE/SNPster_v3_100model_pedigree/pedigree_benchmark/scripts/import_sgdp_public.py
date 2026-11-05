#!/usr/bin/env python3
"""Import only public, metadata-disjoint SGDP sequencing likelihood controls."""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.request

from audit_production_ancestry import digest, write_json


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def aliases(row):
    values = set()
    for field in ("Illumina_ID", "Sample_ID", "Sample_ID(Aliases)", "SGDP_ID"):
        row[field].encode("ascii")
        values.update(value for value in re.split(r"[,;\s]+", row[field]) if value)
    return values


def select_public(records, forbidden, protocol):
    selected = []
    for row in records:
        if (row["Embargo"] != protocol["allowed_access"]["Embargo"]
                or row[list(row)[-1]] != protocol["allowed_access"]["category"]
                or row["Population_ID"] not in protocol["candidate_population_sources"]):
            continue
        names = aliases(row)
        if names & forbidden:
            continue
        selected.append(dict(
            sample_id=row["SGDP_ID"] + ".DG", sgdp_id=row["SGDP_ID"],
            illumina_id=row["Illumina_ID"], population=row["Population_ID"],
            source=protocol["candidate_population_sources"][row["Population_ID"]],
            aliases=sorted(names), access="FullyPublic/X",
        ))
    if len({row["sample_id"] for row in selected}) != len(selected):
        raise ValueError("Duplicate public sample identities")
    seen = set()
    for row in selected:
        names = set(row["aliases"])
        if seen & names:
            raise ValueError("Aliases connect two selected SGDP individuals")
        seen.update(names)
    return sorted(selected, key=lambda row: row["sample_id"])


def metadata(project, protocol):
    preflight_path = project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
    preflight = json.loads(preflight_path.read_text())
    forbidden = set(preflight["references"][0]["sample_ids"])
    for cohort in preflight["cohorts"].values():
        forbidden.update(row["sample_id"] for row in (
            cohort["reference_overlap"] + cohort["direct_id_disjoint_people"]
        ))
    raw = fetch(protocol["metadata_url"])
    Path("metadata.original.txt").write_bytes(raw)
    records = list(csv.DictReader(
        io.StringIO(raw.decode(protocol["metadata_encoding"])), delimiter="\t"
    ))
    selected = select_public(records, forbidden, protocol)
    counts = dict(Counter(row["source"] for row in selected))
    if counts != protocol["expected_metadata_disjoint_counts"]:
        raise ValueError(f"Public donor counts changed: {counts}")
    base = protocol["likelihood_base_url"]
    Path("likelihood_readme.txt").write_bytes(fetch(base + "Readme"))
    Path("publisher_md5sums.txt").write_bytes(fetch(base + "md5sums"))
    write_json("selected.json", dict(
        people=selected, counts=counts, genotype_relatedness_verified=False,
        source_metadata_sha256=digest("metadata.original.txt"),
        production_preflight_sha256=digest(preflight_path),
        production_unchanged=True,
    ))


def download(chromosome, protocol, checksums):
    name = f"chr.sgdp.pub.{chromosome}.bcf"
    md5s = {}
    for line in checksums.read_text().splitlines():
        fields = line.split()
        if len(fields) != 2:
            raise ValueError("Invalid publisher checksum row")
        md5s[fields[1].lstrip("*").removeprefix("./")] = fields[0]
    expected = md5s[name]
    url = protocol["likelihood_base_url"] + name
    checksum = hashlib.md5()
    partial = Path(name + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as output:
        expected_bytes = int(response.headers["Content-Length"])
        headers = dict(response.headers)
        while block := response.read(1024 * 1024):
            checksum.update(block)
            output.write(block)
    if partial.stat().st_size != expected_bytes or checksum.hexdigest() != expected:
        raise ValueError("SGDP download byte length or publisher MD5 mismatch")
    partial.rename(name)
    write_json("download.json", dict(
        url=url, bytes=expected_bytes, publisher_md5=expected,
        sha256=digest(name), response_headers=headers,
    ))
    print(f"Verified public chromosome {chromosome}: {expected_bytes} bytes", flush=True)


def raw_control(protocol):
    base = protocol["raw_control_base_url"]
    listing = fetch(base + "vcfs.variants.public_samples.279samples.file.list")
    Path("public_tar_members.txt").write_bytes(listing)
    desired = protocol["raw_control_illumina_id"] + ".annotated.nh2.variants.vcf.gz"
    members, offset = {}, 0
    for line in listing.decode("ascii").splitlines():
        fields = line.split()
        if len(fields) != 6 or fields[5] in members:
            raise ValueError("Invalid or duplicate public tar member")
        size = int(fields[2])
        members[fields[5]] = (offset, size)
        offset += 512 + ((size + 511) // 512) * 512
    archive_url = base + "vcfs.variants.public_samples.279samples.tar"

    def extract(name, output):
        start, size = members[name]
        request = urllib.request.Request(archive_url, headers={
            "Range": f"bytes={start}-{start + 512 + size - 1}",
            "Accept-Encoding": "identity",
        })
        with urllib.request.urlopen(request, timeout=120) as response:
            content_range = response.headers.get("Content-Range", "")
            if response.status != 206 or not content_range.startswith(
                f"bytes {start}-{start + 512 + size - 1}/"
            ):
                raise ValueError("Public tar server did not honor the exact byte range")
            header = tarfile.TarInfo.frombuf(response.read(512), "utf-8", "strict")
            if header.name != name or header.size != size:
                raise ValueError("Public tar index/header mismatch")
            remaining = size
            with Path(output).open("wb") as handle:
                while remaining:
                    block = response.read(min(1024 * 1024, remaining))
                    if not block:
                        raise ValueError("Truncated public tar member")
                    handle.write(block)
                    remaining -= len(block)
        return dict(url=archive_url, member=name, offset=start, bytes=size)

    provenance = extract(desired, "raw_control.vcf.gz")
    extract(desired + ".md5sum", "raw_control.md5sum")
    expected = Path("raw_control.md5sum").read_text().split()[0]
    with Path("raw_control.vcf.gz").open("rb") as handle:
        checksum = hashlib.md5()
        while block := handle.read(1024 * 1024):
            checksum.update(block)
    if checksum.hexdigest() != expected:
        raise ValueError("Original SGDP control publisher MD5 mismatch")
    write_json("raw_control.json", dict(
        **provenance, publisher_md5=expected, sha256=digest("raw_control.vcf.gz")
    ))


def call_from_likelihoods(pl, ad, protocol):
    if "." in pl or "." in ad:
        return "./.", ".", "."
    likelihoods = [int(value) for value in pl.split(",")]
    depths = [int(value) for value in ad.split(",")]
    if len(likelihoods) != 3 or len(depths) != 2 or min(likelihoods + depths) < 0:
        raise ValueError("Invalid biallelic diploid PL/AD shape or value")
    depth = sum(depths)
    ranked = sorted(enumerate(likelihoods), key=lambda item: item[1])
    gap = ranked[1][1] - ranked[0][1]
    if (not protocol["min_depth"] <= depth <= protocol["max_depth"]
            or gap < protocol["min_likelihood_gap_phred"]):
        return "./.", str(gap), str(depth)
    return ("0/0", "0/1", "1/1")[ranked[0][0]], str(gap), str(depth)


def basecalls(source, selected_path, protocol, chromosome):
    people = json.loads(selected_path.read_text())["people"]
    ids = [row["sample_id"] for row in people]
    actual = subprocess.check_output(["bcftools", "query", "-l", str(source)], text=True).splitlines()
    if not set(ids) <= set(actual) or len(actual) != len(set(actual)):
        raise ValueError("Public BCF lacks selected identities or has duplicate names")
    command = ["bcftools", "query", "-s", ",".join(ids),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%AD]\n", str(source)]
    write_json("likelihood_query.command.json", command)
    called, sites = Counter(), 0
    with Path("basecalls.vcf").open("w") as output:
        output.write("##fileformat=VCFv4.2\n")
        output.write("##source=SGDP_mpileup_PL_AD_only_no_imputed_GT_DS_GP\n")
        output.write(f"##contig=<ID={chromosome}>\n")
        for name, number, kind, description in (
            ("GT", "1", "String", "Unphased maximum-likelihood sequencing genotype"),
            ("GQ", "1", "Integer", "Gap between smallest two raw PL values"),
            ("DP", "1", "Integer", "Sum of raw biallelic AD"),
            ("PL", "G", "Integer", "Original sequencing genotype likelihoods"),
            ("AD", "R", "Integer", "Original sequencing allelic depths"),
        ):
            output.write(f'##FORMAT=<ID={name},Number={number},Type={kind},Description="{description}">\n')
        output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
        try:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                chrom, pos, ref, alt = fields[:4]
                if len(fields) != 4 + 2 * len(ids) or chrom != str(chromosome):
                    raise ValueError("Invalid likelihood chromosome/sample axis")
                if ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
                    continue
                values = []
                for i, sample in enumerate(ids):
                    pl, ad = fields[4 + 2 * i:6 + 2 * i]
                    gt, gq, depth = call_from_likelihoods(pl, ad, protocol)
                    called[sample] += gt != "./."
                    values.append(":".join((gt, gq, depth, pl, ad)))
                output.write("\t".join(
                    [chrom, pos, ".", ref, alt, ".", ".", ".", "GT:GQ:DP:PL:AD"] + values
                ) + "\n")
                sites += 1
            if process.wait():
                raise ValueError("Raw PL/AD query failed")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    if not sites:
        raise ValueError("No biallelic sequencing SNPs")
    fractions = {sample: called[sample] / sites for sample in ids}
    failing = {
        sample: value for sample, value in fractions.items()
        if value < protocol["minimum_called_fraction"]
    }
    if failing:
        write_json("failed_callability.json", dict(sites=sites, failing=failing))
        raise ValueError(f"Sequencing callability below predeclared threshold: {failing}")
    subprocess.run(["bgzip", "basecalls.vcf"], check=True)
    subprocess.run(["bcftools", "index", "basecalls.vcf.gz"], check=True)
    write_json("basecalls.json", dict(
        chromosome=chromosome, sites=sites, called=dict(called),
        called_fraction=fractions,
        source_sha256=digest(source), output_sha256=digest("basecalls.vcf.gz"),
        genotypes_from="mpileup_PL_and_AD_only", input_imputed_GT_DS_GP_used=False,
        coordinate_harmonization_complete=False, relatedness_verified=False,
    ))


def verify(calls, control, selected_path, protocol, chromosome):
    subprocess.run(["bcftools", "index", str(control)], check=True)
    sample = protocol["raw_control_sgdp_id"] + ".DG"
    raw = subprocess.check_output([
        "bcftools", "query", "-r", str(chromosome),
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%FL]\n", str(control),
    ], text=True)
    truth = {}
    for line in raw.splitlines():
        chrom, pos, ref, alt, gt, gq, level = line.split("\t")
        if (len(ref) != 1 or len(alt) != 1 or not level.isdigit() or gq == "."
                or int(level) < protocol["raw_control_min_filter_level"]
                or int(gq) < protocol["raw_control_min_genotype_quality"]
                or gt not in ("0/0", "0/1", "1/0", "1/1")):
            continue
        key = (chrom, pos, ref, alt)
        if key in truth:
            raise ValueError("Duplicate original-call allele key")
        truth[key] = sorted(gt.split("/"))
    text = subprocess.check_output([
        "bcftools", "query", "-s", sample,
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(calls),
    ], text=True)
    matches, total = 0, 0
    for line in text.splitlines():
        *key, gt = line.split("\t")
        if tuple(key) in truth and gt != "./.":
            total += 1
            matches += sorted(gt.split("/")) == truth[tuple(key)]
    if total < protocol["raw_control_min_shared_calls"]:
        raise ValueError(f"Insufficient raw-call validation anchors: {total}")
    concordance = matches / total
    if concordance < protocol["raw_control_min_concordance"]:
        raise ValueError(f"Independent sequencing calls fail concordance: {concordance}")
    result = dict(
        status="raw_control_concordance_passed_not_ancestry_execution_ready",
        chromosome=chromosome, n_compared=total, concordance=concordance,
        raw_control_sha256=digest(control), basecalls_sha256=digest(calls),
        selected_sha256=digest(selected_path),
        required_ancestry_resolutions=protocol["required_ancestry_resolutions"],
        scientific_ancestry_launch_authorized=False,
        remaining_gates=protocol["remaining_gates"],
    )
    write_json("verification.json", result)
    print(json.dumps(result, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("metadata", "download", "raw_control", "basecalls", "verify"))
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--chromosome", type=int, default=22)
    parser.add_argument("--checksums", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--selected", type=Path)
    parser.add_argument("--control", type=Path)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    if args.stage == "metadata":
        metadata(args.project, protocol)
    elif args.stage == "download":
        download(args.chromosome, protocol, args.checksums)
    elif args.stage == "raw_control":
        raw_control(protocol)
    elif args.stage == "basecalls":
        basecalls(args.source, args.selected, protocol, args.chromosome)
    else:
        verify(args.source, args.control, args.selected, protocol, args.chromosome)


if __name__ == "__main__":
    main()
