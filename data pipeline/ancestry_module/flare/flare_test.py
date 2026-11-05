
"""
Run FLARE local ancestry inference on phased Beagle-imputed VCFs.

Expected setup
--------------
Target:
    One Beagle-imputed VCF per chromosome.
    Genotypes must be phased and contain no missing alleles.

Reference:
    One phased HGDP+1kGP reference VCF per chromosome.

Reference panel:
    Two whitespace-separated columns:

        SAMPLE_ID    POPULATION

    Example:
        HGDP00001    EUR
        HGDP00002    EUR
        HGDP00003    SAS
        HGDP00004    EAS

Genetic map:
    PLINK-format genetic map in cM.
    Can either be a single multi-chromosome map or a pattern containing
    {chrom}, e.g.:

        /data/maps/plink.chr{chrom}.GRCh38.map

Usage example
-------------
python run_flare.py \
    --gt \
        IMPID307.chr1.imputed.vcf.gz \
        IMPID307.chr2.imputed.vcf.gz \
        IMPID307.chr3.imputed.vcf.gz \
    --reference-pattern \
        "/data/references/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz" \
    --ref-panel \
        /data/references/hgdp_1kg_panel_flare.txt \
    --map \
        /data/references/plink.GRCh38.map \
    --flare-jar \
        /data/flare.jar \
    --out-dir flare_output \
    --threads 8 \
    --memory-gb 32
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path


AUTOSOMES = {str(x) for x in range(1, 23)}


def open_text(path: Path):
    """Open normal or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def infer_chromosome(vcf_path: Path) -> str:
    """
    Read the first variant in a VCF and return its chromosome.

    chr1 -> 1
    1    -> 1
    """
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            chrom = line.split("\t", 1)[0]
            chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)

            if chrom not in AUTOSOMES:
                raise ValueError(
                    f"{vcf_path}: expected autosome 1-22, found chromosome '{chrom}'"
                )

            return chrom

    raise ValueError(f"No variants found in {vcf_path}")


def chromosome_sort_key(chrom: str) -> int:
    return int(chrom)


def validate_target_vcf(
    vcf_path: Path,
    max_records: int = 10_000,
) -> None:
    """
    Quick validation of the Beagle VCF.

    Checks a subset of records to ensure:
      * GT exists
      * genotypes are phased
      * genotypes are not missing

    FLARE itself will perform complete validation.
    """

    checked = 0

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                raise ValueError(
                    f"{vcf_path}: VCF contains no genotype/sample columns."
                )

            format_fields = fields[8].split(":")

            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                raise ValueError(
                    f"{vcf_path}: FORMAT field does not contain GT."
                )

            for sample in fields[9:]:
                sample_fields = sample.split(":")

                if gt_index >= len(sample_fields):
                    raise ValueError(
                        f"{vcf_path}: malformed genotype FORMAT field."
                    )

                gt = sample_fields[gt_index]

                if "." in gt:
                    raise ValueError(
                        f"{vcf_path}: missing genotype encountered: {gt}"
                    )

                if "|" not in gt:
                    raise ValueError(
                        f"{vcf_path}: unphased genotype encountered: {gt}"
                    )

            checked += 1

            if checked >= max_records:
                break

    if checked == 0:
        raise ValueError(f"{vcf_path}: contains no variant records.")

    print(
        f"[OK] {vcf_path.name}: "
        f"checked {checked:,} phased, non-missing records."
    )


def resolve_pattern(pattern: str, chrom: str) -> Path:
    """
    Resolve a path containing {chrom}.

    If {chrom} isn't present, return path unchanged.
    """
    return Path(pattern.format(chrom=chrom))


def run_command(command: list[str]) -> None:
    print("\n[COMMAND]")
    print(" ".join(command))
    print()

    subprocess.run(
        command,
        check=True,
    )


def run_flare(
    *,
    flare_jar: Path,
    target_vcf: Path,
    reference_vcf: Path,
    ref_panel: Path,
    genetic_map: Path,
    out_prefix: Path,
    threads: int,
    memory_gb: int,
    probs: bool,
    seed: int,
    model: Path | None = None,
    em: bool = True,
) -> None:

    command = [
        "java",
        f"-Xmx{memory_gb}g",
        "-jar",
        str(flare_jar),
        f"gt={target_vcf}",
        f"ref={reference_vcf}",
        f"ref-panel={ref_panel}",
        f"map={genetic_map}",
        f"out={out_prefix}",
        f"nthreads={threads}",
        f"seed={seed}",
        f"probs={'true' if probs else 'false'}",

        # We're feeding the Beagle-imputed dataset rather than
        # the original microarray.
        "array=false",

        f"em={'true' if em else 'false'}",
    ]

    if model is not None:
        command.append(f"model={model}")

    run_command(command)


def parse_global_ancestry(
    path: Path,
) -> tuple[list[str], list[list[str]]]:

    with gzip.open(path, "rt") as handle:
        lines = [
            line.strip().split()
            for line in handle
            if line.strip()
        ]

    if not lines:
        raise ValueError(f"Empty FLARE global ancestry file: {path}")

    return lines[0], lines[1:]


def write_chromosome_summary(
    results: list[dict],
    output_path: Path,
) -> None:
    """
    Combine per-chromosome FLARE global ancestry estimates
    into one long-format TSV.

    This does NOT average chromosomes, because eventually you
    should derive genome-wide ancestry from the local ancestry
    calls/tract lengths.
    """

    records = []

    for result in results:
        global_file = Path(result["global_ancestry"])

        header, rows = parse_global_ancestry(global_file)

        ancestries = header[1:]

        for row in rows:
            sample = row[0]
            values = row[1:]

            for ancestry, value in zip(ancestries, values):
                records.append(
                    (
                        result["chromosome"],
                        sample,
                        ancestry,
                        value,
                    )
                )

    with open(output_path, "w") as handle:
        handle.write(
            "chromosome\tsample\tancestry\tproportion\n"
        )

        for chrom, sample, ancestry, value in records:
            handle.write(
                f"{chrom}\t{sample}\t{ancestry}\t{value}\n"
            )


def main():

    parser = argparse.ArgumentParser(
        description="Run FLARE on Beagle-imputed chromosome VCFs."
    )

    parser.add_argument(
        "--gt",
        nargs="+",
        required=True,
        help="One or more Beagle-imputed chromosome VCFs.",
    )

    parser.add_argument(
        "--reference-pattern",
        required=True,
        help=(
            "Reference VCF path. May contain {chrom}, e.g. "
            "'hgdp1kgp_chr{chrom}.vcf.gz'"
        ),
    )

    parser.add_argument(
        "--ref-panel",
        required=True,
        help="Two-column FLARE reference panel file.",
    )

    parser.add_argument(
        "--map",
        required=True,
        help=(
            "PLINK-format genetic map. May contain {chrom}, "
            "or may be one multi-chromosome map."
        ),
    )

    parser.add_argument(
        "--flare-jar",
        required=True,
    )

    parser.add_argument(
        "--out-dir",
        required=True,
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--memory-gb",
        type=int,
        default=16,
    )

    parser.add_argument(
        "--train-chrom",
        default="1",
        help=(
            "Chromosome used to estimate FLARE model parameters. "
            "Default: 1."
        ),
    )

    parser.add_argument(
        "--model",
        default=None,
        help=(
            "Optional pretrained FLARE model. "
            "If supplied, EM is disabled for every chromosome."
        ),
    )

    parser.add_argument(
        "--probs",
        action="store_true",
        help="Write ancestry posterior probabilities.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260925,
    )

    parser.add_argument(
        "--skip-validation",
        action="store_true",
    )

    args = parser.parse_args()

    flare_jar = Path(args.flare_jar)
    ref_panel = Path(args.ref_panel)
    output_dir = Path(args.out_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in [flare_jar, ref_panel]:
        if not required.exists():
            raise FileNotFoundError(required)

    #
    # Discover chromosomes
    #

    targets: dict[str, Path] = {}

    for filename in args.gt:

        target = Path(filename)

        if not target.exists():
            raise FileNotFoundError(target)

        chrom = infer_chromosome(target)

        if chrom in targets:
            raise ValueError(
                f"Multiple target VCFs supplied for chromosome {chrom}"
            )

        targets[chrom] = target

    targets = dict(
        sorted(
            targets.items(),
            key=lambda x: chromosome_sort_key(x[0]),
        )
    )

    print("\nDetected chromosomes:")
    for chrom, vcf in targets.items():
        print(f"  chr{chrom}: {vcf}")

    #
    # Quick Beagle validation
    #

    if not args.skip_validation:
        print("\nValidating phased Beagle genotypes...")

        for target in targets.values():
            validate_target_vcf(target)

    #
    # Resolve reference/map files
    #

    inputs = {}

    for chrom, target in targets.items():

        reference = resolve_pattern(
            args.reference_pattern,
            chrom,
        )

        genetic_map = resolve_pattern(
            args.map,
            chrom,
        )

        if not reference.exists():
            raise FileNotFoundError(
                f"Reference VCF missing for chr{chrom}: {reference}"
            )

        if not genetic_map.exists():
            raise FileNotFoundError(
                f"Genetic map missing for chr{chrom}: {genetic_map}"
            )

        inputs[chrom] = {
            "target": target,
            "reference": reference,
            "map": genetic_map,
        }

    #
    # Determine model strategy
    #

    external_model = (
        Path(args.model)
        if args.model is not None
        else None
    )

    if external_model is not None:

        if not external_model.exists():
            raise FileNotFoundError(external_model)

        print(
            f"\nUsing pretrained FLARE model: "
            f"{external_model}"
        )

        model = external_model

    else:

        train_chrom = str(args.train_chrom)

        if train_chrom not in inputs:
            train_chrom = next(iter(inputs))

            print(
                f"\nRequested training chromosome was not supplied. "
                f"Using chr{train_chrom}."
            )

        train_input = inputs[train_chrom]

        train_dir = output_dir / f"chr{train_chrom}"
        train_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        train_prefix = (
            train_dir /
            f"flare.chr{train_chrom}"
        )

        print(
            f"\nEstimating FLARE model parameters using chr{train_chrom}..."
        )

        run_flare(
            flare_jar=flare_jar,
            target_vcf=train_input["target"],
            reference_vcf=train_input["reference"],
            ref_panel=ref_panel,
            genetic_map=train_input["map"],
            out_prefix=train_prefix,
            threads=args.threads,
            memory_gb=args.memory_gb,
            probs=args.probs,
            seed=args.seed,
            model=None,
            em=True,
        )

        model = Path(
            f"{train_prefix}.model"
        )

        if not model.exists():
            raise RuntimeError(
                f"FLARE did not create expected model: {model}"
            )

    #
    # Run remaining chromosomes
    #

    results = []

    for chrom, input_files in inputs.items():

        chromosome_dir = output_dir / f"chr{chrom}"
        chromosome_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        prefix = (
            chromosome_dir /
            f"flare.chr{chrom}"
        )

        expected_vcf = Path(
            f"{prefix}.anc.vcf.gz"
        )

        #
        # If this chromosome generated the model above,
        # FLARE has already run.
        #
        already_done = (
            external_model is None
            and chrom == str(args.train_chrom)
            and expected_vcf.exists()
        )

        if not already_done:

            print(
                f"\nRunning FLARE chr{chrom}..."
            )

            run_flare(
                flare_jar=flare_jar,
                target_vcf=input_files["target"],
                reference_vcf=input_files["reference"],
                ref_panel=ref_panel,
                genetic_map=input_files["map"],
                out_prefix=prefix,
                threads=args.threads,
                memory_gb=args.memory_gb,
                probs=args.probs,
                seed=args.seed,
                model=model,
                em=False,
            )

        global_file = Path(
            f"{prefix}.global.anc.gz"
        )

        anc_vcf = Path(
            f"{prefix}.anc.vcf.gz"
        )

        if not global_file.exists():
            raise RuntimeError(
                f"Missing FLARE global output: {global_file}"
            )

        if not anc_vcf.exists():
            raise RuntimeError(
                f"Missing FLARE ancestry VCF: {anc_vcf}"
            )

        results.append(
            {
                "chromosome": chrom,
                "target": str(input_files["target"]),
                "reference": str(input_files["reference"]),
                "ancestry_vcf": str(anc_vcf),
                "global_ancestry": str(global_file),
            }
        )

    #
    # Combined chromosome summary
    #

    summary_file = (
        output_dir /
        "global_ancestry_by_chromosome.tsv"
    )

    write_chromosome_summary(
        results,
        summary_file,
    )

    metadata = {
        "model": str(model),
        "reference_panel": str(ref_panel),
        "chromosomes": results,
    }

    with open(
        output_dir / "flare_run.json",
        "w",
    ) as handle:
        json.dump(
            metadata,
            handle,
            indent=2,
        )

    print("\nFLARE completed successfully.")
    print(f"Model:   {model}")
    print(f"Summary: {summary_file}")
    print(f"Output:  {output_dir}")


if __name__ == "__main__":
    main()