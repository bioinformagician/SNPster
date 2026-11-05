#!/usr/bin/env python3
"""Recover redundant public-source SNP records without weakening genotype checks."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from audit_production_ancestry import digest, write_json
from verify_public_sgdp import verify
from verify_sgdp_full_manifest import require_hashes, source_hashes


def select_record(records):
    unique = set(records)
    supported = {line for line in unique
                 if any(value != "./.:.:.:.:." for value in line.rstrip("\n").split("\t")[9:])}
    if len(supported) > 1:
        raise ValueError("Conflicting duplicate sequencing records; no automatic repair allowed")
    return next(iter(supported)) if supported else records[0]


def duplicates(calls, chromosome):
    process = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(calls)],
        stdout=subprocess.PIPE, text=True,
    )
    previous, seen, repeated, sites = 0, set(), Counter(), 0
    try:
        for line in process.stdout:
            chrom, pos, ref, alt = line.rstrip("\n").split("\t")
            position = int(pos)
            if chrom != str(chromosome) or position < previous or position < 1:
                raise ValueError("Invalid or unsorted public genotype coordinates")
            if position != previous:
                seen.clear()
            key = (chrom, pos, ref, alt)
            if key in seen:
                repeated[key] += 1
            seen.add(key)
            previous = position
            sites += 1
        if process.wait():
            raise ValueError("Public duplicate-key scan failed")
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait()
    return sites, repeated


def bound_hashes(record, root, smoke):
    hashes = {}
    for base, manifest in ((root, record["sha256"]), (smoke, record["shared_sha256"])):
        for relative, expected in manifest.items():
            filename = (base / relative).resolve()
            if not filename.is_relative_to(base.resolve()):
                raise ValueError("Verification evidence escapes published input directory")
            hashes[str(filename)] = expected
    return hashes


def repair(root, smoke, chromosome, repeated, output):
    regions = output / "duplicate_regions.tsv"
    regions.write_text("".join(
        f"{chrom}\t{pos}\t{pos}\n" for chrom, pos in
        sorted({key[:2] for key in repeated}, key=lambda key: int(key[1]))
    ))
    text = subprocess.check_output(
        ["bcftools", "view", "-H", "-R", str(regions), str(root / "basecalls/basecalls.vcf.gz")],
        text=True,
    )
    records = {key: [] for key in repeated}
    for line in text.splitlines(keepends=True):
        fields = line.rstrip("\n").split("\t")
        key = (fields[0], fields[1], fields[3], fields[4])
        if key in records:
            records[key].append(line)
    if any(len(records[key]) != count + 1 for key, count in repeated.items()):
        raise ValueError("Indexed duplicate extraction disagrees with full scan")
    chosen = {key: select_record(lines) for key, lines in records.items()}
    evidence = [dict(
        key=list(key), input_records=lines, retained_record=chosen[key],
        removed_records=len(lines) - 1,
    ) for key, lines in records.items()]
    write_json(output / "duplicate_evidence.json", evidence)
    normalized = output / "normalized"
    for folder in ("basecalls", "source", "verification"):
        (normalized / folder).mkdir(parents=True)
    for relative in ("source/download.json", "basecalls/likelihood_query.command.json"):
        shutil.copyfile(root / relative, normalized / relative)
    selected = json.loads((smoke / "metadata/selected.json").read_text())
    ids = [person["sample_id"] for person in selected["people"]]
    original_summary = json.loads((root / "basecalls/basecalls.json").read_text())
    original_calls = root / "basecalls/basecalls.vcf.gz"
    if digest(original_calls) != original_summary["output_sha256"]:
        raise ValueError("Original genotype evidence hash mismatch")
    counts, retained, removed, emitted = Counter(), 0, 0, set()
    retained_body = hashlib.sha256()
    target = normalized / "basecalls/basecalls.vcf"
    with gzip.open(original_calls, "rt") as source, target.open("w") as destination:
        for line in source:
            if line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(ids) or fields[8] != "GT:GQ:DP:PL:AD":
                raise ValueError("Unexpected original genotype format/sample axis")
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in chosen:
                if line not in records[key]:
                    raise ValueError("Duplicate record changed since evidence extraction")
                if line != chosen[key] or key in emitted:
                    removed += 1
                    continue
                emitted.add(key)
            destination.write(line)
            retained_body.update(line.encode("ascii"))
            for sample, value in zip(ids, fields[9:]):
                counts[sample] += value.split(":", 1)[0] != "./."
            retained += 1
    if (removed != sum(repeated.values()) or emitted != set(chosen)
            or retained + removed != original_summary["sites"]):
        raise ValueError("Duplicate-repair site accounting differs")
    fractions = {sample: counts[sample] / retained for sample in ids}
    protocol = Path(__file__).resolve().parents[2] / "config/sgdp_public_import.json"
    if min(fractions.values()) < json.loads(protocol.read_text())["minimum_called_fraction"]:
        raise ValueError("Corrected calls fail original callability threshold")
    subprocess.run(["bgzip", str(target)], check=True)
    calls = target.with_suffix(".vcf.gz")
    subprocess.run(["bcftools", "index", str(calls)], check=True)
    check_body = hashlib.sha256()
    with gzip.open(calls, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                check_body.update(line.encode("ascii"))
    if check_body.digest() != retained_body.digest():
        raise ValueError("Retained genotype records changed during compression")
    summary = dict(original_summary, sites=retained, called=dict(counts),
                   called_fraction=fractions, output_sha256=digest(calls))
    write_json(normalized / "basecalls/basecalls.json", summary)
    raw = output / "raw_control.vcf.gz"
    shutil.copyfile(smoke / "raw_control/raw_control.vcf.gz", raw)
    subprocess.run([
        "python3", str(Path(__file__).with_name("import_sgdp_public.py")), "verify",
        "--protocol", str(protocol), "--chromosome", str(chromosome),
        "--source", str(calls.resolve()), "--control", str(raw.resolve()),
        "--selected", str((smoke / "metadata/selected.json").resolve()),
    ], cwd=normalized / "verification", check=True)
    raw.unlink()
    raw.with_suffix(raw.suffix + ".csi").unlink()
    return normalized, evidence


def recover(root, smoke, project, chromosome, previous_report, output):
    output.mkdir()
    if previous_report is not None:
        record = json.loads(previous_report.read_text())
        if (record["status"] != "published_genotypes_verified_not_ancestry_execution_ready"
                or record["chromosome"] != chromosome
                or record["scientific_ancestry_launch_authorized"] is not False):
            raise ValueError("Invalid previously verified chromosome report")
        require_hashes(bound_hashes(record, root, smoke))
        shutil.copyfile(previous_report, output / "independent.json")
        details = dict(
            chromosome=chromosome, status="reused_independently_verified_chromosome",
            original_root=str(root), original_sites=record["sites"], removed_records=0,
            previous_report=str(previous_report), previous_report_sha256=digest(previous_report),
        )
    else:
        sites, repeated = duplicates(root / "basecalls/basecalls.vcf.gz", chromosome)
        data_root, evidence = (repair(root, smoke, chromosome, repeated, output)
                               if repeated else (root, []))
        verify(data_root, project, chromosome, smoke, output / "independent.json")
        details = dict(
            chromosome=chromosome, status="corrected_redundant_duplicates" if repeated else "no_duplicates",
            original_root=str(root), original_sites=sites, removed_records=sum(repeated.values()),
            original_vcf_sha256=digest(root / "basecalls/basecalls.vcf.gz"),
            duplicate_keys=[item["key"] for item in evidence],
            correction_policy="one_unchanged_supported_record_or_one_identical_missing_record; conflicts_fail",
        )
    write_json(output / "recovery.json", details)


def collect(import_run, smoke, previous, recovery_run, folders, output):
    reports = [json.loads((folder / "independent.json").read_text()) for folder in folders]
    if len(reports) != 22 or sorted(row["chromosome"] for row in reports) != list(range(1, 23)):
        raise ValueError("Exactly one independent report per autosome is required")
    ids = [row["sample_id"] for row in json.loads((smoke / "metadata/selected.json").read_text())["people"]]
    qc = {row["chromosome"]: row for row in
          json.loads((import_run / "summary/full_import.json").read_text())["raw_call_qc"]}
    inputs, data_roots, reused, removed = {}, {}, [], 0
    for run in (import_run, previous, recovery_run):
        inputs.update(source_hashes(run / "source_sha256.txt"))
        inputs[str(run / "source_sha256.txt")] = digest(run / "source_sha256.txt")
    require_hashes(inputs)
    for folder, record in zip(folders, reports):
        chromosome = record["chromosome"]
        details = json.loads((folder / "recovery.json").read_text())
        original_root = smoke if chromosome == 22 else import_run / f"chromosomes/chr{chromosome}"
        if (details["chromosome"] != chromosome or details["original_root"] != str(original_root)
                or details["original_sites"] - details["removed_records"] != record["sites"]
                or record["sites"] < 1 or record["sample_ids"] != ids
                or record["status"] != "published_genotypes_verified_not_ancestry_execution_ready"
                or record["scientific_ancestry_launch_authorized"] is not False
                or record["genotype_cells_independently_rederived"] != record["sites"] * len(ids)
                or record["original_raw_control_concordance"] != qc[chromosome]["concordance"]
                or record["original_raw_control_compared"] != qc[chromosome]["n_compared"]):
            raise ValueError(f"Recovery or genotype verification differs: chr{chromosome}")
        published = recovery_run / f"chromosomes/chr{chromosome}"
        corrected = details["status"] == "corrected_redundant_duplicates"
        root = folder / "normalized" if corrected else original_root
        require_hashes(bound_hashes(record, root, smoke))
        data_roots[str(chromosome)] = str(published / "normalized" if corrected else original_root)
        for filename, expected in bound_hashes(record, root, smoke).items():
            path = Path(filename)
            destination = published / path.relative_to(folder.resolve()) if path.is_relative_to(folder.resolve()) else path
            inputs[str(destination)] = expected
        if "original_vcf_sha256" in details:
            original = str(original_root / "basecalls/basecalls.vcf.gz")
            require_hashes({original: details["original_vcf_sha256"]})
            inputs[original] = details["original_vcf_sha256"]
        if "previous_report" in details:
            require_hashes({details["previous_report"]: details["previous_report_sha256"]})
            inputs[details["previous_report"]] = details["previous_report_sha256"]
            reused.append(chromosome)
        for path in folder.rglob("*"):
            if path.is_file():
                inputs[str(published / path.relative_to(folder))] = digest(path)
        removed += details["removed_records"]
    result = dict(
        status="all_autosome_public_genotypes_independently_verified_not_ancestry_ready",
        import_run=str(import_run), smoke_run=str(smoke), chromosomes=list(range(1, 23)),
        sample_ids=ids, data_roots=data_roots, reused_verified_chromosomes=sorted(reused),
        removed_redundant_records=removed, sites=sum(row["sites"] for row in reports),
        genotype_cells_independently_rederived=sum(row["genotype_cells_independently_rederived"] for row in reports),
        original_raw_control_compared=sum(row["original_raw_control_compared"] for row in reports),
        minimum_chromosome_raw_concordance=min(row["original_raw_control_concordance"] for row in reports),
        chromosome_reports=sorted(reports, key=lambda row: row["chromosome"]),
        coordinate_harmonization_complete=False, genotype_relatedness_verified=False,
        scientific_ancestry_launch_authorized=False, inputs_sha256=inputs,
    )
    write_json(output, result)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("inputs_sha256", "chromosome_reports")}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("recover", "collect"))
    parser.add_argument("--root", type=Path)
    parser.add_argument("--smoke", type=Path, required=True)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--chromosome", type=int)
    parser.add_argument("--previous-report", type=Path)
    parser.add_argument("--import-run", type=Path)
    parser.add_argument("--previous-run", type=Path)
    parser.add_argument("--recovery-run", type=Path)
    parser.add_argument("--folders", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "recover":
        recover(args.root, args.smoke, args.project, args.chromosome, args.previous_report, args.output)
    else:
        collect(args.import_run, args.smoke, args.previous_run, args.recovery_run, args.folders, args.output)
