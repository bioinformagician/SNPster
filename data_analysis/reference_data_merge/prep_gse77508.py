#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import re
import subprocess
from pathlib import Path


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def comp(a):
    return a.translate(COMPLEMENT)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_snp_field(value):
    """
    Parse Illumina manifest SNP field like:
        [A/G]
        [C/T]
    """
    if not value:
        return None

    m = re.match(r"\[([ACGT])/([ACGT])\]", value.strip())

    if not m:
        return None

    return m.group(1), m.group(2)


def find_manifest_header(path):
    """
    Illumina CSV manifests often contain metadata before [Assay].
    Returns CSV header + rows starting at the actual assay section.
    """
    fh = open(path, "rt", newline="")

    for line in fh:
        if line.strip() == "[Assay]":
            break

    reader = csv.DictReader(fh)

    return fh, reader


def read_manifest(path):
    """
    Returns:
        rsid -> metadata

    Expected useful columns are usually:
        Name
        SNP
        Chr
        MapInfo

    Different OmniExpress releases vary slightly.
    """

    fh, reader = find_manifest_header(path)

    manifest = {}

    for row in reader:

        name = (
            row.get("Name")
            or row.get("IlmnID")
            or row.get("SNP_Name")
        )

        snp_field = row.get("SNP")

        chrom = (
            row.get("Chr")
            or row.get("Chromosome")
        )

        pos = (
            row.get("MapInfo")
            or row.get("Position")
        )

        if not name or not chrom or not pos:
            continue

        # GSE matrix already uses rsIDs as Name.
        rsid = name.strip()

        if not rsid.startswith("rs"):
            continue

        alleles = parse_snp_field(snp_field)

        if alleles is None:
            continue

        a1, a2 = alleles

        # Initially discard strand-ambiguous SNPs.
        if {a1, a2} in ({"A", "T"}, {"C", "G"}):
            continue

        chrom = chrom.replace("chr", "")

        if chrom not in {str(i) for i in range(1, 23)}:
            continue

        try:
            pos = int(pos)
        except ValueError:
            continue

        manifest[rsid] = {
            "chrom": chrom,
            "pos": pos,
            "A": a1,
            "B": a2,
        }

    fh.close()

    print(f"Manifest SNPs retained: {len(manifest):,}")

    return manifest


def read_matrix_header(matrix_path):
    with open_text(matrix_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")

    genotype_columns = []

    for i, col in enumerate(header):
        if col.endswith(".GType"):
            sample = col.removesuffix(".GType")
            genotype_columns.append((i, sample))

    print(f"Samples found: {len(genotype_columns)}")

    return genotype_columns


def query_reference_variant(reference_vcf, chrom, pos, rsid=None):
    """
    Query HGDP+1kGP at one position.

    Returns:
        (ref, alt)
    or None.

    We use the existing HGDP+1kGP reference as the final authority
    for GRCh38 REF/ALT.
    """

    region = f"{chrom}:{pos}-{pos}"

    cmd = [
        "bcftools",
        "query",
        "-r",
        region,
        "-f",
        "%ID\t%REF\t%ALT\n",
        reference_vcf,
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )

    for line in result.stdout.splitlines():
        vid, ref, alt = line.split("\t")

        if "," in alt:
            continue

        if rsid and vid not in {rsid, "."}:
            # Prefer matching rsID but don't necessarily require it.
            continue

        if len(ref) != 1 or len(alt) != 1:
            continue

        return ref, alt

    # Second pass without requiring ID match
    for line in result.stdout.splitlines():
        vid, ref, alt = line.split("\t")

        if "," in alt:
            continue

        if len(ref) == 1 and len(alt) == 1:
            return ref, alt

    return None


def orient_illumina_alleles(a, b, ref, alt):
    """
    Determine whether Illumina A/B correspond directly or on the
    complementary strand to HGDP+1kGP REF/ALT.

    Returns:
        allele_A_genomic,
        allele_B_genomic
    """

    if {a, b} == {ref, alt}:
        return a, b

    ca = comp(a)
    cb = comp(b)

    if {ca, cb} == {ref, alt}:
        return ca, cb

    return None


def gtype_to_gt(gtype, allele_a, allele_b, ref, alt):
    """
    Convert AA / AB / BB to VCF GT relative to REF / ALT.
    """

    gtype = gtype.strip()

    if gtype in {"", "NC", "--", "00", "NoCall"}:
        return "./."

    if gtype == "AA":
        alleles = [allele_a, allele_a]

    elif gtype in {"AB", "BA"}:
        alleles = [allele_a, allele_b]

    elif gtype == "BB":
        alleles = [allele_b, allele_b]

    else:
        return "./."

    gt = []

    for allele in alleles:

        if allele == ref:
            gt.append("0")

        elif allele == alt:
            gt.append("1")

        else:
            return "./."

    return "/".join(gt)


def write_vcf_header(handle, samples, chrom):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##source=GSE77508_OmniExpress\n")
    handle.write(f"##contig=<ID={chrom}>\n")
    handle.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )

    cols = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        *samples,
    ]

    handle.write("\t".join(cols) + "\n")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--matrix",
        required=True,
        help="GSE77508_ProcessedMatrix.txt.gz",
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="Illumina OmniExpress CSV manifest",
    )

    parser.add_argument(
        "--reference-pattern",
        required=True,
        help="HGDP+1kGP reference VCF pattern with {chrom}",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
    )

    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading manifest...")

    manifest = read_manifest(args.manifest)

    genotype_columns = read_matrix_header(args.matrix)

    samples = [x[1] for x in genotype_columns]

    #
    # Open one VCF per chromosome
    #

    handles = {}

    for chrom in map(str, range(1, 23)):

        path = out_dir / f"GSE77508.chr{chrom}.vcf"

        fh = open(path, "w")

        write_vcf_header(
            fh,
            samples,
            chrom,
        )

        handles[chrom] = fh

    stats = {
        "matrix_rows": 0,
        "not_manifest": 0,
        "not_reference": 0,
        "allele_mismatch": 0,
        "written": 0,
    }

    print("Processing genotype matrix...")

    with open_text(args.matrix) as fh:

        reader = csv.reader(fh, delimiter="\t")

        header = next(reader)

        for row in reader:

            stats["matrix_rows"] += 1

            rsid = row[0]

            meta = manifest.get(rsid)

            if meta is None:
                stats["not_manifest"] += 1
                continue

            chrom = meta["chrom"]
            pos = meta["pos"]

            ref_vcf = args.reference_pattern.format(
                chrom=chrom
            )

            reference = query_reference_variant(
                ref_vcf,
                chrom,
                pos,
                rsid,
            )

            if reference is None:
                stats["not_reference"] += 1
                continue

            ref, alt = reference

            oriented = orient_illumina_alleles(
                meta["A"],
                meta["B"],
                ref,
                alt,
            )

            if oriented is None:
                stats["allele_mismatch"] += 1
                continue

            allele_a, allele_b = oriented

            genotypes = []

            for idx, sample in genotype_columns:

                gt = gtype_to_gt(
                    row[idx],
                    allele_a,
                    allele_b,
                    ref,
                    alt,
                )

                genotypes.append(gt)

            record = [
                chrom,
                str(pos),
                rsid,
                ref,
                alt,
                ".",
                "PASS",
                ".",
                "GT",
                *genotypes,
            ]

            handles[chrom].write(
                "\t".join(record) + "\n"
            )

            stats["written"] += 1

            if stats["matrix_rows"] % 10_000 == 0:
                print(
                    f"Processed {stats['matrix_rows']:,} rows, "
                    f"written {stats['written']:,}"
                )

    for fh in handles.values():
        fh.close()

    print("\nFinished")
    print("--------")

    for key, value in stats.items():
        print(f"{key:20s}: {value:,}")


if __name__ == "__main__":
    main()