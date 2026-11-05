#!/usr/bin/env python3

import csv
import re
import subprocess
from pathlib import Path
from collections import Counter

# ============================================================
# Paths
# ============================================================

HGDP_META = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/"
    "gnomad_meta_updated.tsv"
)

HGDP_VCF = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/"
    "hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5."
    "SNV_biallelic.numericCHR.vcf.gz"
)

GSE77508_META = Path(
    "/srv/dependencies/ancestry/gse77508/"
    "GSE77508_metadata.tsv"
)

GSE77508_VCF = Path(
    "/srv/dependencies/ancestry/gse77508/"
    "GSE77508.chr1.phased.vcf.gz"
)

GSE73996_META = Path(
    "/srv/dependencies/ancestry/gse73996/"
    "GSE73996_metadata.tsv"
)

GSE73996_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/gse73996/"
    "selected/phased/gse73996.chr1.phased.vcf.gz"
)

KAZAKH_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/kazakh/"
    "vcf/phased/kazakh.chr1.phased.vcf.gz"
)

# Original -> sanitized mapping created during previous merge.
ORIGINAL_NAMES = Path(
    "/srv/dependencies/ancestry/flare_reference_v1/output/"
    "original_sample_names.txt"
)

SANITIZED_NAMES = Path(
    "/srv/dependencies/ancestry/flare_reference_v1/output/"
    "sanitized/sample_names.txt"
)

# Current dense merged reference: already sanitized.
MERGED_VCF = Path(
    "/home/frederik/github_projects/SNPster/data_analysis/"
    "reference_data_merge/imputed/merged_reference/output/"
    "flare_reference.chr1.vcf.gz"
)

OUTDIR = Path(
    "/home/frederik/github_projects/SNPster/data_analysis/"
    "reference_data_merge/imputed/merged_reference/output"
)

FINE_PANEL = OUTDIR / "flare_reference.fine.ref.panel"

FINE_METADATA = OUTDIR / "fine_population_metadata.tsv"

SAMPLE_METADATA = OUTDIR / "fine_sample_metadata.tsv"


# ============================================================
# Helpers
# ============================================================

def vcf_samples(path):
    out = subprocess.check_output(
        ["bcftools", "query", "-l", str(path)],
        text=True
    )
    return [
        x.rstrip("\n")
        for x in out.splitlines()
        if x.strip()
    ]


def clean(value):
    if value is None:
        return ""

    value = value.strip()

    if value.upper() in {
        "",
        "NA",
        "N/A",
        "NONE",
        "NULL",
        ".",
    }:
        return ""

    return value


def safe_label(value):
    value = clean(value)

    value = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        value
    )

    return value.strip("_")


def find_column(columns, possibilities):
    lookup = {
        x.lower(): x
        for x in columns
    }

    for candidate in possibilities:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]

    return None


# ============================================================
# Original -> sanitized sample IDs
# ============================================================

original_names = [
    x.rstrip("\n")
    for x in ORIGINAL_NAMES.read_text().splitlines()
    if x.strip()
]

sanitized_names = [
    x.rstrip("\n")
    for x in SANITIZED_NAMES.read_text().splitlines()
    if x.strip()
]

if len(original_names) != len(sanitized_names):
    raise RuntimeError(
        "original/sanitized sample-name counts differ"
    )

rename = dict(
    zip(
        original_names,
        sanitized_names
    )
)

if len(rename) != len(original_names):
    raise RuntimeError(
        "Duplicate original sample IDs in rename mapping"
    )


# ============================================================
# Assign fine populations using ORIGINAL IDs
# ============================================================

assignments = {}


def add(
    sample,
    population,
    source,
    region="",
    project=""
):

    population = clean(population)
    region = clean(region)
    project = clean(project)

    if not population:
        raise RuntimeError(
            f"No population label for {sample} ({source})"
        )

    if sample in assignments:
        raise RuntimeError(
            f"Duplicate assignment: {sample}"
        )

    # Keep the dataset/source namespace.
    #
    # FLARE2 can later learn that, for example,
    # GSE73996_Kazakh and Kazakh2025_Kazakh
    # are closely related instead of us forcing them
    # into one panel beforehand.
    label = (
        f"{safe_label(source)}_"
        f"{safe_label(population)}"
    )

    assignments[sample] = {
        "fine_label": label,
        "source": source,
        "population": population,
        "region": region,
        "project": project,
    }


# ============================================================
# HGDP + 1kGP
#
# Use:
#   sample ID:
#       s
#
#   fine population:
#       hgdp_tgp_meta.Population
#
#   coarse context:
#       hgdp_tgp_meta.Genetic.region
#
#   project:
#       hgdp_tgp_meta.Project
# ============================================================

hgdp_vcf_samples = set(
    vcf_samples(HGDP_VCF)
)

with open(HGDP_META) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    required = {
        "s",
        "hgdp_tgp_meta.Population",
        "hgdp_tgp_meta.Genetic.region",
        "hgdp_tgp_meta.Project",
    }

    missing_columns = (
        required - set(reader.fieldnames)
    )

    if missing_columns:
        raise RuntimeError(
            "gnomAD metadata missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    for row in reader:

        sample = clean(row["s"])

        if sample not in hgdp_vcf_samples:
            continue

        population = clean(
            row["hgdp_tgp_meta.Population"]
        )

        region = clean(
            row["hgdp_tgp_meta.Genetic.region"]
        )

        project = clean(
            row["hgdp_tgp_meta.Project"]
        )

        # Fallbacks for any unusual row.
        if not population:

            population = clean(
                row.get("population", "")
            )

        if not population:

            population = clean(
                row.get(
                    "project_meta.v2_known_subpop",
                    ""
                )
            )

        if not population:

            population = clean(
                row.get(
                    "project_meta.v2_known_pop",
                    ""
                )
            )

        if not region:

            region = clean(
                row.get("region", "")
            )

        add(
            sample=sample,
            population=population,
            source=(
                "HGDP"
                if project.upper() == "HGDP"
                else "1KG"
            ),
            region=region,
            project=project,
        )


# Check HGDP+1kGP coverage immediately.

missing_hgdp = (
    hgdp_vcf_samples - set(assignments)
)

if missing_hgdp:
    raise RuntimeError(
        "HGDP+1kGP samples missing metadata:\n"
        + "\n".join(
            sorted(missing_hgdp)[:100]
        )
    )


# ============================================================
# GSE77508
# ============================================================

gse77508_samples = set(
    vcf_samples(GSE77508_VCF)
)

with open(GSE77508_META) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    for row in reader:

        sample = clean(
            row["sample"]
        )

        if sample not in gse77508_samples:
            continue

        add(
            sample=sample,
            population=row["population"],
            source="GSE77508",
            region=row["region"],
            project="GSE77508",
        )


missing_77508 = (
    gse77508_samples - set(assignments)
)

if missing_77508:
    raise RuntimeError(
        "GSE77508 samples missing metadata:\n"
        + "\n".join(
            sorted(missing_77508)[:100]
        )
    )


# ============================================================
# GSE73996
# ============================================================

gse73996_samples = set(
    vcf_samples(GSE73996_VCF)
)

with open(GSE73996_META) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    columns = reader.fieldnames

    sample_col = find_column(
        columns,
        [
            "sample",
            "sample_id",
            "sampleid",
            "id",
            "title",
        ]
    )

    pop_col = find_column(
        columns,
        [
            "population",
            "pop",
        ]
    )

    region_col = find_column(
        columns,
        [
            "region",
            "geographic_region",
        ]
    )

    if sample_col is None:
        raise RuntimeError(
            "Cannot identify GSE73996 sample column"
        )

    if pop_col is None:
        raise RuntimeError(
            "Cannot identify GSE73996 population column"
        )

    for row in reader:

        sample = clean(
            row[sample_col]
        )

        if sample not in gse73996_samples:
            continue

        region = (
            row[region_col]
            if region_col
            else ""
        )

        add(
            sample=sample,
            population=row[pop_col],
            source="GSE73996",
            region=region,
            project="GSE73996",
        )


missing_73996 = (
    gse73996_samples - set(assignments)
)

if missing_73996:
    raise RuntimeError(
        "GSE73996 samples missing metadata:\n"
        + "\n".join(
            sorted(missing_73996)[:100]
        )
    )


# ============================================================
# Kazakh 2025
# ============================================================

kazakh_samples = set(
    vcf_samples(KAZAKH_VCF)
)

for sample in kazakh_samples:

    add(
        sample=sample,
        population="Kazakh",
        source="Kazakh2025",
        region="Central_Asia",
        project="PRJEB89820",
    )


# ============================================================
# Map original source IDs -> sanitized merged IDs
# ============================================================

sanitized_assignments = {}

unmapped_original = []

for original_sample, meta in assignments.items():

    if original_sample not in rename:

        unmapped_original.append(
            original_sample
        )

        continue

    sanitized = rename[
        original_sample
    ]

    if sanitized in sanitized_assignments:

        raise RuntimeError(
            f"Duplicate sanitized ID: {sanitized}"
        )

    sanitized_assignments[
        sanitized
    ] = meta


if unmapped_original:
    raise RuntimeError(
        "Original source IDs missing from "
        "rename mapping:\n"
        + "\n".join(
            sorted(unmapped_original)[:100]
        )
    )


# ============================================================
# Validate against dense merged reference
# ============================================================

merged_samples = vcf_samples(
    MERGED_VCF
)

merged_set = set(
    merged_samples
)

assigned_set = set(
    sanitized_assignments
)

missing = (
    merged_set - assigned_set
)

extra = (
    assigned_set - merged_set
)

if missing:

    raise RuntimeError(
        "Merged-reference samples without "
        "fine population:\n"
        + "\n".join(
            sorted(missing)[:100]
        )
    )

if extra:

    raise RuntimeError(
        "Fine-panel samples not in merged "
        "reference:\n"
        + "\n".join(
            sorted(extra)[:100]
        )
    )

if len(merged_samples) != 4886:

    raise RuntimeError(
        f"Expected 4886 merged samples, "
        f"found {len(merged_samples)}"
    )


# ============================================================
# Write FLARE2 fine reference panel
#
# Exactly:
#
# sanitized_sample_id   fine_population
# ============================================================

with open(FINE_PANEL, "w") as handle:

    for sample in merged_samples:

        label = (
            sanitized_assignments[
                sample
            ]["fine_label"]
        )

        handle.write(
            f"{sample}\t{label}\n"
        )


# ============================================================
# Sample-level metadata
# ============================================================

with open(
    SAMPLE_METADATA,
    "w"
) as handle:

    handle.write(
        "sample\tfine_population\t"
        "source\tpopulation\tregion\tproject\n"
    )

    for sample in merged_samples:

        meta = sanitized_assignments[
            sample
        ]

        handle.write(
            "\t".join([
                sample,
                meta["fine_label"],
                meta["source"],
                meta["population"],
                meta["region"],
                meta["project"],
            ])
            + "\n"
        )


# ============================================================
# Unique fine-population metadata
# ============================================================

population_info = {}

for meta in sanitized_assignments.values():

    label = meta["fine_label"]

    info = (
        meta["source"],
        meta["population"],
        meta["region"],
        meta["project"],
    )

    if (
        label in population_info
        and population_info[label] != info
    ):
        raise RuntimeError(
            f"Inconsistent metadata for "
            f"{label}"
        )

    population_info[
        label
    ] = info


counts = Counter(
    x["fine_label"]
    for x in sanitized_assignments.values()
)


with open(
    FINE_METADATA,
    "w"
) as handle:

    handle.write(
        "fine_population\tn_samples\t"
        "source\tpopulation\tregion\tproject\n"
    )

    for label in sorted(
        population_info
    ):

        source, population, region, project = (
            population_info[label]
        )

        handle.write(
            "\t".join([
                label,
                str(counts[label]),
                source,
                population,
                region,
                project,
            ])
            + "\n"
        )


# ============================================================
# Summary
# ============================================================

print()
print("Fine reference panel successfully created")
print("------------------------------------------")
print(f"Samples          : {len(merged_samples):,}")
print(f"Fine populations : {len(population_info):,}")
print()
print(f"Panel            : {FINE_PANEL}")
print(f"Population meta  : {FINE_METADATA}")
print(f"Sample meta      : {SAMPLE_METADATA}")
print()

print("Largest fine populations:")

for label, n in counts.most_common(30):
    print(
        f"  {label:35s} {n:5d}"
    )
