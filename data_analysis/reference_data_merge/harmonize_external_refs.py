#!/usr/bin/env python3

import argparse
import gzip
import subprocess
from collections import defaultdict
from pathlib import Path

from pyliftover import LiftOver


BASES = {"A", "C", "G", "T"}
COMP = str.maketrans("ACGT", "TGCA")


def complement(x):
    return x.translate(COMP)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalize_chrom(chrom):
    chrom = chrom.removeprefix("chr")

    if chrom in {str(i) for i in range(1, 23)}:
        return chrom

    return None


def read_samples(vcf):
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]

    raise RuntimeError(f"No VCF header found: {vcf}")


def collect_source_positions(inputs):
    """
    Collect unique GRCh37 autosomal positions across all input datasets.
    """

    positions = set()

    for name, vcf in inputs.items():

        n_records = 0

        with open_text(vcf) as fh:

            for line in fh:

                if line.startswith("#"):
                    continue

                fields = line.split("\t", 3)

                chrom = normalize_chrom(fields[0])

                if chrom is None:
                    continue

                pos = int(fields[1])

                positions.add((chrom, pos))
                n_records += 1

        print(
            f"{name}: {n_records:,} autosomal records"
        )

    print(
        f"Unique GRCh37 positions across datasets: "
        f"{len(positions):,}"
    )

    return positions


def liftover_positions(source_positions, chain):
    """
    Lift GRCh37 1-based positions to GRCh38.

    pyliftover uses 0-based coordinates internally.
    """

    print("\nLoading liftover chain...")

    lo = LiftOver(str(chain))

    mapping = {}

    stats = {
        "source_positions": len(source_positions),
        "lifted": 0,
        "unmapped": 0,
        "multimapped": 0,
        "non_autosomal_target": 0,
    }

    for i, (chrom, pos) in enumerate(source_positions, 1):

        hits = lo.convert_coordinate(
            f"chr{chrom}",
            pos - 1
        )

        valid = []

        for hit in hits:

            target_chrom = normalize_chrom(hit[0])

            if target_chrom is None:
                continue

            target_pos = int(hit[1]) + 1
            strand = hit[2]

            valid.append(
                (target_chrom, target_pos, strand)
            )

        # Deduplicate identical mappings.
        valid = list(dict.fromkeys(valid))

        if len(valid) == 0:

            if hits:
                stats["non_autosomal_target"] += 1
            else:
                stats["unmapped"] += 1

            mapping[(chrom, pos)] = None

        elif len(valid) > 1:

            stats["multimapped"] += 1
            mapping[(chrom, pos)] = None

        else:

            mapping[(chrom, pos)] = valid[0]
            stats["lifted"] += 1

        if i % 100_000 == 0:
            print(
                f"Lifted {i:,}/{len(source_positions):,}"
            )

    print("\nLiftover stats:")

    for k, v in stats.items():
        print(f"{k:24s}: {v:,}")

    return mapping


def collect_target_positions(lift_map):
    target = defaultdict(set)

    for result in lift_map.values():

        if result is None:
            continue

        chrom, pos, strand = result
        target[chrom].add(pos)

    return target


def build_hgdp_lookup(
    target_positions,
    reference_pattern,
    work_dir,
):
    """
    Query HGDP only at lifted GRCh38 positions.

    Returns:
      (chr, pos) -> [(REF, ALT), ...]
    """

    print("\nBuilding canonical HGDP+1kGP allele lookup...")

    lookup = defaultdict(list)

    region_dir = work_dir / "regions"
    region_dir.mkdir(parents=True, exist_ok=True)

    for chrom in range(1, 23):

        chrom = str(chrom)

        positions = sorted(
            target_positions.get(chrom, [])
        )

        print(
            f"  chr{chrom}: querying "
            f"{len(positions):,} lifted positions"
        )

        if not positions:
            continue

        regions = region_dir / f"chr{chrom}.regions.tsv"

        with open(regions, "w") as fh:
            for pos in positions:
                fh.write(
                    f"{chrom}\t{pos}\t{pos}\n"
                )

        ref_vcf = reference_pattern.format(
            chrom=chrom
        )

        cmd = [
            "bcftools",
            "query",
            "-R",
            str(regions),
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\n",
            ref_vcf,
        ]

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        for line in proc.stdout:

            rchrom, pos, ref, alt = (
                line.rstrip().split("\t")
            )

            if "," in alt:
                continue

            ref = ref.upper()
            alt = alt.upper()

            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in BASES
                or alt not in BASES
            ):
                continue

            key = (
                normalize_chrom(rchrom),
                int(pos),
            )

            pair = (ref, alt)

            if pair not in lookup[key]:
                lookup[key].append(pair)

        stderr = proc.stderr.read()
        rc = proc.wait()

        if rc != 0:
            raise RuntimeError(
                f"bcftools failed for chr{chrom}:\n"
                f"{stderr}"
            )

    print(
        f"\nHGDP positions matched: "
        f"{len(lookup):,}"
    )

    return lookup


def determine_orientation(
    src_ref,
    src_alt,
    dst_ref,
    dst_alt,
):
    """
    Return whether GT allele indices need swapping.

    False:
       0 remains 0, 1 remains 1

    True:
       0 becomes 1, 1 becomes 0

    None:
       alleles incompatible
    """

    if (
        src_ref == dst_ref
        and src_alt == dst_alt
    ):
        return False

    if (
        src_ref == dst_alt
        and src_alt == dst_ref
    ):
        return True

    cref = complement(src_ref)
    calt = complement(src_alt)

    if (
        cref == dst_ref
        and calt == dst_alt
    ):
        return False

    if (
        cref == dst_alt
        and calt == dst_ref
    ):
        return True

    return None


def remap_gt(gt, swap):

    if gt in {".", "./.", ".|."}:
        return "./."

    separator = "|" if "|" in gt else "/"

    alleles = gt.replace("|", "/").split("/")

    result = []

    for allele in alleles:

        if allele == ".":
            result.append(".")
            continue

        if allele not in {"0", "1"}:
            return "./."

        if swap:
            result.append(
                "1" if allele == "0" else "0"
            )
        else:
            result.append(allele)

    return separator.join(result)


def write_header(handle, dataset, chrom, samples):

    handle.write("##fileformat=VCFv4.2\n")
    handle.write(
        f"##source={dataset}_lifted_and_harmonized\n"
    )
    handle.write("##reference=GRCh38\n")
    handle.write(
        f"##contig=<ID={chrom}>\n"
    )

    handle.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype">\n'
    )

    columns = [
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

    handle.write("\t".join(columns) + "\n")


def harmonize_dataset(
    name,
    vcf,
    lift_map,
    hgdp_lookup,
    output_root,
):

    print(f"\n===== {name} =====")

    samples = read_samples(vcf)

    print(f"Samples: {len(samples):,}")

    dataset_dir = output_root / name
    raw_dir = dataset_dir / "raw"
    final_dir = dataset_dir / "vcf"

    raw_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    final_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    handles = {}

    for chrom in range(1, 23):

        chrom = str(chrom)

        fh = open(
            raw_dir / f"{name}.chr{chrom}.vcf",
            "w"
        )

        write_header(
            fh,
            name,
            chrom,
            samples,
        )

        handles[chrom] = fh

    stats = {
        "input_variants": 0,
        "non_autosomal": 0,
        "liftover_failed": 0,
        "not_in_hgdp": 0,
        "non_biallelic_snp": 0,
        "ambiguous_dropped": 0,
        "allele_mismatch": 0,
        "ambiguous_hgdp_match": 0,
        "duplicate_target": 0,
        "written": 0,
    }

    seen = set()

    with open_text(vcf) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            stats["input_variants"] += 1

            fields = line.rstrip().split("\t")

            chrom37 = normalize_chrom(
                fields[0]
            )

            if chrom37 is None:
                stats["non_autosomal"] += 1
                continue

            pos37 = int(fields[1])

            lifted = lift_map.get(
                (chrom37, pos37)
            )

            if lifted is None:
                stats["liftover_failed"] += 1
                continue

            chrom38, pos38, lift_strand = lifted

            source_ref = fields[3].upper()
            source_alt = fields[4].upper()

            if (
                "," in source_alt
                or len(source_ref) != 1
                or len(source_alt) != 1
                or source_ref not in BASES
                or source_alt not in BASES
            ):
                stats["non_biallelic_snp"] += 1
                continue

            # Drop strand-ambiguous markers initially.
            if {source_ref, source_alt} in (
                {"A", "T"},
                {"C", "G"},
            ):
                stats["ambiguous_dropped"] += 1
                continue

            candidates = hgdp_lookup.get(
                (chrom38, pos38)
            )

            if not candidates:
                stats["not_in_hgdp"] += 1
                continue

            compatible = []

            for dst_ref, dst_alt in candidates:

                swap = determine_orientation(
                    source_ref,
                    source_alt,
                    dst_ref,
                    dst_alt,
                )

                if swap is not None:

                    compatible.append(
                        (
                            dst_ref,
                            dst_alt,
                            swap,
                        )
                    )

            if len(compatible) == 0:

                stats["allele_mismatch"] += 1
                continue

            if len(compatible) > 1:

                stats["ambiguous_hgdp_match"] += 1
                continue

            dst_ref, dst_alt, swap = compatible[0]

            target_key = (
                chrom38,
                pos38,
                dst_ref,
                dst_alt,
            )

            if target_key in seen:
                stats["duplicate_target"] += 1
                continue

            format_fields = fields[8].split(":")

            if "GT" not in format_fields:
                continue

            gt_idx = format_fields.index("GT")

            genotypes = []

            for sample in fields[9:]:

                parts = sample.split(":")

                if gt_idx >= len(parts):
                    gt = "./."
                else:
                    gt = parts[gt_idx]

                genotypes.append(
                    remap_gt(gt, swap)
                )

            source_id = fields[2]

            record = [
                chrom38,
                str(pos38),
                source_id,
                dst_ref,
                dst_alt,
                ".",
                "PASS",
                ".",
                "GT",
                *genotypes,
            ]

            handles[chrom38].write(
                "\t".join(record) + "\n"
            )

            seen.add(target_key)

            stats["written"] += 1

            if stats["input_variants"] % 100_000 == 0:

                print(
                    f"Processed "
                    f"{stats['input_variants']:,}; "
                    f"written "
                    f"{stats['written']:,}"
                )

    for fh in handles.values():
        fh.close()

    print("\nSorting/compressing/indexing...")

    for chrom in range(1, 23):

        raw = (
            raw_dir
            / f"{name}.chr{chrom}.vcf"
        )

        final = (
            final_dir
            / f"{name}.chr{chrom}.grch38.vcf.gz"
        )

        subprocess.run(
            [
                "bcftools",
                "sort",
                str(raw),
                "-Oz",
                "-o",
                str(final),
            ],
            check=True,
        )

        subprocess.run(
            [
                "bcftools",
                "index",
                "-f",
                str(final),
            ],
            check=True,
        )

        raw.unlink()

    stats_file = (
        dataset_dir
        / "harmonization_stats.tsv"
    )

    with open(stats_file, "w") as fh:

        fh.write("metric\tvalue\n")

        for key, value in stats.items():
            fh.write(
                f"{key}\t{value}\n"
            )

    print("\n--------")

    for key, value in stats.items():
        print(
            f"{key:24s}: {value:,}"
        )


def parse_inputs(values):

    result = {}

    for value in values:

        name, path = value.split("=", 1)

        result[name] = Path(path)

    return result


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        action="append",
        required=True,
    )

    parser.add_argument(
        "--chain",
        required=True,
    )

    parser.add_argument(
        "--reference-pattern",
        required=True,
    )

    parser.add_argument(
        "--out-root",
        required=True,
    )

    args = parser.parse_args()

    inputs = parse_inputs(args.input)

    out_root = Path(args.out_root)

    out_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    source_positions = collect_source_positions(
        inputs
    )

    lift_map = liftover_positions(
        source_positions,
        Path(args.chain),
    )

    target_positions = collect_target_positions(
        lift_map
    )

    hgdp_lookup = build_hgdp_lookup(
        target_positions,
        args.reference_pattern,
        out_root,
    )

    for name, vcf in inputs.items():

        harmonize_dataset(
            name,
            vcf,
            lift_map,
            hgdp_lookup,
            out_root,
        )


if __name__ == "__main__":
    main()