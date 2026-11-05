#!/usr/bin/env python3
"""Normalize public sequencing calls without changing diploid genotype evidence."""
import argparse
from bisect import bisect_right
from collections import Counter
from functools import lru_cache
import gzip
import importlib.metadata
import json
import mmap
from pathlib import Path
import subprocess

from pyliftover import LiftOver
import zstandard

from audit_production_ancestry import digest, write_json
from verify_sgdp_full_manifest import require_hashes

COMPLEMENT = str.maketrans("ACGT", "TGCA")


class Fasta:
    def __init__(self, path):
        self.entries = {}
        for line in Path(str(path) + ".fai").read_text().splitlines():
            name, length, offset, bases, width = line.split("\t")[:5]
            self.entries[name] = tuple(map(int, (length, offset, bases, width)))
        self.handle = Path(path).open("rb")
        self.data = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)

    def base(self, chromosome, position):
        length, offset, bases, width = self.entries[chromosome]
        if not 1 <= position <= length:
            raise ValueError(f"FASTA coordinate out of bounds: {chromosome}:{position}")
        zero = position - 1
        return chr(self.data[offset + zero // bases * width + zero % bases]).upper()

    def close(self):
        self.data.close()
        self.handle.close()


def index_fasta(path):
    entries, name, length, offset, bases, width, short = [], None, 0, 0, 0, 0, False
    with path.open("rb") as handle:
        while line := handle.readline():
            if line.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, bases, width))
                name = line[1:].split()[0].decode("ascii")
                if name in {entry[0] for entry in entries}:
                    raise ValueError("Duplicate FASTA contig")
                length, offset, bases, width, short = 0, handle.tell(), 0, 0, False
            else:
                value = line.rstrip(b"\r\n")
                if name is None or not value or short:
                    raise ValueError("Invalid FASTA line layout")
                if not bases:
                    bases, width = len(value), len(line)
                elif len(value) != bases or len(line) != width:
                    short = True
                length += len(value)
        if name is not None:
            entries.append((name, length, offset, bases, width))
    Path(str(path) + ".fai").write_text("".join(
        "\t".join(map(str, entry)) + "\n" for entry in entries))


class IndependentChain:
    """Direct UCSC block arithmetic, separate from the production lift-over library."""
    def __init__(self, path, chromosome):
        blocks = []
        with gzip.open(path, "rt") as handle:
            selected = False
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "chain":
                    if len(fields) != 13 or fields[4] != "+":
                        raise ValueError("Unsupported chain header")
                    selected = fields[2] == "chr" + chromosome
                    source, target = int(fields[5]), int(fields[10])
                    target_name, target_size, strand = fields[7], int(fields[8]), fields[9]
                else:
                    size = int(fields[0])
                    if selected:
                        blocks.append((source, source + size, target_name, target, target_size, strand))
                    if len(fields) == 3:
                        source += size + int(fields[1])
                        target += size + int(fields[2])
                    elif len(fields) != 1:
                        raise ValueError("Invalid chain block")
        self.blocks = sorted(blocks)
        self.starts = [block[0] for block in self.blocks]
        self.maximum = []
        for block in self.blocks:
            self.maximum.append(max(block[1], self.maximum[-1] if self.maximum else 0))

    def convert_coordinate(self, chromosome, zero):
        results = []
        index = bisect_right(self.starts, zero) - 1
        while index >= 0 and self.maximum[index] > zero:
            first, last, name, target, size, strand = self.blocks[index]
            if first <= zero < last:
                mapped = target + zero - first
                results.append((name, size - 1 - mapped if strand == "-" else mapped, strand))
            index -= 1
        return results


def mapping(chromosome, position, ref, alt, lift, source, target):
    if source.base(chromosome, position) != ref:
        return None, "source_reference_mismatch"
    hits = lift.convert_coordinate("chr" + chromosome, position - 1) or []
    if not hits:
        return None, "unmapped"
    if len(hits) != 1:
        return None, "ambiguous_mapping"
    name, zero, strand = hits[0][:3]
    if name != "chr" + chromosome:
        return None, "different_contig"
    if strand not in ("+", "-"):
        raise ValueError("Unsupported mapping strand")
    new_ref, new_alt = (ref, alt) if strand == "+" else (
        ref.translate(COMPLEMENT), alt.translate(COMPLEMENT))
    reference = target.base(chromosome, zero + 1)
    if reference == new_ref:
        swapped = False
    elif reference == new_alt:
        new_ref, new_alt, swapped = new_alt, new_ref, True
    else:
        return None, "target_reference_mismatch"
    return (zero + 1, new_ref, new_alt, strand, swapped), None


@lru_cache(maxsize=100000)
def swap_call(value):
    gt, gq, dp, pl, ad = value.split(":")
    if gt != "./.":
        gt = "/".join(sorted(str(1 - int(allele)) for allele in gt.split("/")))
    pl = ",".join(reversed(pl.split(",")))
    ad = ",".join(reversed(ad.split(",")))
    return ":".join((gt, gq, dp, pl, ad))


def vcf_records(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t")


def prepare(project, source_zst, target, chain):
    protocol = json.loads((project / "config/sgdp_normalization.json").read_text())
    manifest_path = project / protocol["input_manifest"]
    manifest = json.loads(manifest_path.read_text())
    if (manifest["status"] != "all_autosome_public_genotypes_independently_verified_not_ancestry_ready"
            or manifest["chromosomes"] != list(range(1, 23))
            or manifest["scientific_ancestry_launch_authorized"] is not False):
        raise ValueError("Verified all-autosome public genotype manifest required")
    require_hashes(manifest["inputs_sha256"])
    if digest(chain) != protocol["chain_sha256"]:
        raise ValueError("Production chain identity changed")
    fasta = Path("hs37d5.fa")
    with source_zst.open("rb") as compressed, fasta.open("wb") as output:
        zstandard.ZstdDecompressor().copy_stream(compressed, output)
    index_fasta(fasta)
    write_json("preflight.json", dict(
        status="sequencing_inputs_and_normalization_resources_bound",
        input_manifest=str(manifest_path), input_manifest_sha256=digest(manifest_path),
        source_compressed_sha256=digest(source_zst), source_fasta_sha256=digest(fasta),
        source_index_sha256=digest(str(fasta) + ".fai"),
        target_fasta=str(target), target_fasta_sha256=digest(target),
        target_index_sha256=digest(str(target) + ".fai"), chain=str(chain), chain_sha256=digest(chain),
        sample_ids=manifest["sample_ids"], pyliftover_version=importlib.metadata.version("pyliftover"),
        zstandard_version=importlib.metadata.version("zstandard"),
        production_preprocessing_changed=False, scientific_ancestry_launch_authorized=False,
    ))
    print("Verified sequencing inputs; source FASTA prepared and both assemblies/chain bound.", flush=True)


def normalize(calls, chromosome, source_path, target_path, chain, protocol, output):
    output.mkdir()
    source, target = Fasta(source_path), Fasta(target_path)
    lift = LiftOver(str(chain))
    counts, sites = Counter(), 0
    unsorted = output / "mapped.vcf"
    with gzip.open(calls, "rt") as original, unsorted.open("w") as destination, \
            gzip.open(output / "rejected.tsv.gz", "wt") as rejected:
        rejected.write("original_row\tchromosome\tposition\tref\talt\treason\n")
        ids = None
        for line in original:
            if line.startswith("##contig="):
                destination.write(f"##contig=<ID={chromosome},length={target.entries[chromosome][0]}>\n")
            elif line.startswith("#CHROM"):
                ids = line.rstrip("\n").split("\t")[9:]
                for name, kind in (("ORIGROW", "Integer"), ("SRCPOS", "Integer"),
                                   ("SRCREF", "String"), ("SRCALT", "String"), ("MAPSTRAND", "String"),
                                   ("REFSWAP", "Integer")):
                    destination.write(f'##INFO=<ID={name},Number=1,Type={kind},Description="Sequencing coordinate provenance">\n')
                destination.write("##reference=GRCh38\n" + line)
            elif line.startswith("#"):
                destination.write(line)
            else:
                fields = line.rstrip("\n").split("\t")
                chrom, pos, _, ref, alt = fields[:5]
                sites += 1
                if chrom != chromosome or len(fields) != 9 + len(ids) or fields[8] != "GT:GQ:DP:PL:AD":
                    raise ValueError("Unexpected source genotype locus, sample or format axis")
                mapped, reason = mapping(chrom, int(pos), ref, alt, lift, source, target)
                if reason is not None:
                    counts[reason] += 1
                    rejected.write(f"{sites}\t{chrom}\t{pos}\t{ref}\t{alt}\t{reason}\n")
                    continue
                position, new_ref, new_alt, strand, swapped = mapped
                counts["mapped"] += 1
                counts["reverse_strand"] += strand == "-"
                counts["ref_alt_swaps"] += swapped
                fields[1], fields[3], fields[4] = str(position), new_ref, new_alt
                fields[7] = f"ORIGROW={sites};SRCPOS={pos};SRCREF={ref};SRCALT={alt};MAPSTRAND={strand};REFSWAP={int(swapped)}"
                if swapped:
                    fields[9:] = [swap_call(value) for value in fields[9:]]
                destination.write("\t".join(fields) + "\n")
                if sites % 1000000 == 0:
                    print(f"chr{chromosome}: normalized {sites:,} sequencing loci", flush=True)
    source.close()
    target.close()
    if ids is None or not sites:
        raise ValueError("Empty sequencing inputs")
    sorted_path = output / "sorted.vcf.gz"
    subprocess.run(["bcftools", "sort", "-m", "256M", "-Oz", "-o", str(sorted_path),
                    str(unsorted)], check=True)
    unsorted.unlink()
    collision_rows, retained, called, key_groups = [], 0, Counter(), {}
    final = output / "normalized.vcf"

    def emit_group(destination):
        nonlocal retained
        for records in key_groups.values():
            if len(records) > 1:
                collision_rows.extend(int(dict(item.split("=", 1) for item in fields[7].split(";"))["ORIGROW"])
                                      for fields in records)
                continue
            fields = records[0]
            destination.write("\t".join(fields) + "\n")
            for sample, value in zip(ids, fields[9:]):
                called[sample] += not value.startswith("./.")
            retained += 1

    with gzip.open(sorted_path, "rt") as stream, final.open("w") as destination:
        previous = None
        for line in stream:
            if line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[1] != previous:
                emit_group(destination)
                key_groups.clear()
                previous = fields[1]
            key_groups.setdefault((fields[3], fields[4]), []).append(fields)
        emit_group(destination)
    sorted_path.unlink()
    if retained / sites < protocol["minimum_retained_fraction"]:
        raise ValueError(f"Insufficient mapped sequencing coverage: {retained}/{sites}")
    if any(called[sample] / retained < protocol["minimum_called_fraction"] for sample in ids):
        raise ValueError("Normalized sequencing calls fail predeclared sample callability")
    write_json(output / "collision_rows.json", sorted(collision_rows))
    subprocess.run(["bgzip", str(final)], check=True)
    calls_out = output / "normalized.vcf.gz"
    subprocess.run(["bcftools", "index", str(calls_out)], check=True)
    write_json(output / "normalization.json", dict(
        status="producer_normalized_independent_verification_required", chromosome=int(chromosome),
        source_vcf=str(calls), source_vcf_sha256=digest(calls),
        source_sites=sites, retained_sites=retained, counts=dict(counts),
        collision_excluded_rows=len(collision_rows), sample_ids=ids, called=dict(called),
        output_sha256=digest(calls_out), scientific_ancestry_launch_authorized=False,
    ))
    print(f"chr{chromosome}: {retained:,}/{sites:,} records retained; independent verification required.", flush=True)


def verify_normalized(original, folder, chromosome, source_path, target_path, chain, output):
    summary = json.loads((folder / "normalization.json").read_text())
    calls = folder / "normalized.vcf.gz"
    if digest(original) != summary["source_vcf_sha256"] or digest(calls) != summary["output_sha256"]:
        raise ValueError("Source or normalized genotype identity changed")
    source, target = Fasta(source_path), Fasta(target_path)
    independent = IndependentChain(chain, chromosome)
    collisions = set(json.loads((folder / "collision_rows.json").read_text()))
    restored = Path(f"chr{chromosome}.source_order.tsv")
    previous_position, keys = 0, set()
    with restored.open("w") as ordered:
        process = subprocess.Popen(["sort", "-n", "-k1,1", "-S", "128M", "-T", "."],
                                   stdin=subprocess.PIPE, stdout=ordered, text=True, env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"})
        try:
            for fields in vcf_records(calls):
                position = int(fields[1])
                if fields[0] != chromosome or position < previous_position:
                    raise ValueError("Normalized coordinates are not chromosome-sorted")
                if position != previous_position:
                    keys.clear()
                key = tuple(fields[3:5])
                if key in keys or target.base(chromosome, position) != fields[3]:
                    raise ValueError("Duplicate normalized allele key or invalid GRCh38 REF")
                keys.add(key)
                previous_position = position
                info = dict(item.split("=", 1) for item in fields[7].split(";"))
                process.stdin.write(info["ORIGROW"] + "\t" + "\t".join(fields) + "\n")
            process.stdin.close()
            if process.wait():
                raise ValueError("Independent source-order restoration failed")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    checked, source_sites, exclusions, called = 0, 0, Counter(), Counter()
    with restored.open() as ordered:
        pending = ordered.readline()
        for row, fields in enumerate(vcf_records(original), 1):
            source_sites = row
            expected, reason = mapping(chromosome, int(fields[1]), fields[3], fields[4],
                                       independent, source, target)
            selected_row = int(pending.split("\t", 1)[0]) if pending else None
            if reason is not None or row in collisions:
                if selected_row == row:
                    raise ValueError("Invalid or colliding source mapping was retained")
                exclusions[reason if reason else "target_collision"] += 1
                continue
            if selected_row != row:
                raise ValueError(f"Valid source row missing or duplicated: chr{chromosome}, row{row}")
            actual = pending.rstrip("\n").split("\t")[1:]
            info = dict(item.split("=", 1) for item in actual[7].split(";"))
            position, ref, alt, strand, swapped = expected
            if (actual[0] != chromosome or int(actual[1]) != position or actual[3:5] != [ref, alt]
                    or actual[8] != fields[8] or len(actual) != len(fields)
                    or info != dict(ORIGROW=str(row), SRCPOS=fields[1], SRCREF=fields[3], SRCALT=fields[4],
                                    MAPSTRAND=strand, REFSWAP=str(int(swapped)))):
                raise ValueError("Independent chain-block mapping or provenance disagrees")
            for sample, old, new in zip(summary["sample_ids"], fields[9:], actual[9:]):
                old_gt, old_gq, old_dp, old_pl, old_ad = old.split(":")
                gt, gq, dp, pl, ad = new.split(":")
                if (gq, dp) != (old_gq, old_dp):
                    raise ValueError("Sequencing quality/depth changed")
                old_alleles = [fields[3], fields[4]]
                if strand == "-":
                    old_alleles = [allele.translate(COMPLEMENT) for allele in old_alleles]
                if old_gt == "./.":
                    if gt != "./.":
                        raise ValueError("Missing sequencing genotype was filled")
                elif (gt not in ("0/0", "0/1", "1/1")
                      or sorted([ref, alt][int(i)] for i in gt.split("/")) !=
                      sorted(old_alleles[int(i)] for i in old_gt.split("/"))):
                    raise ValueError("Diploid nucleotide genotype changed during normalization")
                if (pl.split(",") != (old_pl.split(",")[::-1] if swapped else old_pl.split(","))
                        or ad.split(",") != (old_ad.split(",")[::-1] if swapped else old_ad.split(","))):
                    raise ValueError("Raw sequencing likelihood/depth axis changed")
                called[sample] += gt != "./."
            checked += 1
            pending = ordered.readline()
        if pending:
            raise ValueError("Normalized records remain beyond original source axis")
    restored.unlink()
    source.close()
    target.close()
    ids = subprocess.check_output(["bcftools", "query", "-l", str(calls)], text=True).splitlines()
    if (ids != summary["sample_ids"] or source_sites != summary["source_sites"]
            or checked != summary["retained_sites"] or dict(called) != summary["called"]
            or int(subprocess.check_output(["bcftools", "index", "-n", str(calls)], text=True)) != checked
            or source_sites != checked + sum(exclusions.values())):
        raise ValueError("Independent normalized site/sample/index accounting differs")
    write_json(output, dict(
        status="coordinate_allele_and_genotype_preservation_verified", chromosome=int(chromosome),
        source_sites=source_sites, sites=checked, sample_ids=ids,
        genotype_cells_verified=checked * len(ids), called=dict(called), exclusions=dict(exclusions),
        normalized_sha256=digest(calls), normalized_index_sha256=digest(str(calls) + ".csi"),
        normalization_summary_sha256=digest(folder / "normalization.json"),
        source_vcf_sha256=digest(original), scientific_ancestry_launch_authorized=False,
    ))
    print(f"chr{chromosome}: independently verified {checked:,} loci and all diploid calls.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "normalize", "verify"))
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--calls", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--chromosome")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.project, args.source, args.target, args.chain)
    elif args.stage == "normalize":
        normalize(args.calls, args.chromosome, args.source, args.target, args.chain,
                  json.loads((args.project / "config/sgdp_normalization.json").read_text()), args.output)
    else:
        verify_normalized(args.calls, args.folder, args.chromosome, args.source, args.target,
                          args.chain, args.output)
