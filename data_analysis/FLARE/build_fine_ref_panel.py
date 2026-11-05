#!/usr/bin/env python3

from pathlib import Path
import csv
import subprocess
import re

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

HGDP_FINE_PANEL = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/hgdp_1kg_panel.txt"
)

GSE73996_META = Path(
    "/srv/dependencies/ancestry/gse73996/GSE73996_metadata.tsv"
)

GSE77508_META = Path(
    "/srv/dependencies/ancestry/gse77508/GSE77508_metadata.tsv"
)

HGDP_VCF = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/"
    "hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5."
    "SNV_biallelic.numericCHR.vcf.gz"
)

GSE77508_VCF = Path(
    "/srv/dependencies/ancestry/gse77508/"
    "GSE77508.chr1.phased.vcf.gz"
)

GSE73996_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/gse73996/"
    "selected/phased/gse73996.chr1.phased.vcf.gz"
)

KAZAKH_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/kazakh/"
    "vcf/phased/kazakh.chr1.phased.vcf.gz"
)

MERGED_VCF = Path(
    "/home/frederik/github_projects/SNPster/data_analysis/"
    "reference_data_merge/imputed/merged_reference/output/"
    "flare_reference.chr1.vcf.gz"
)

SANITIZED_NAMES = Path(
    "/srv/dependencies/ancestry/flare_reference_v1/output/"
    "sanitized/sample_names.txt"
)

OUT = Path(
    "/home/frederik/github_projects/SNPster/data_analysis/"
    "reference_data_merge/imputed/merged_reference/output/"
)

OUT.mkdir(parents=True, exist_ok=True)

FINE_PANEL = OUT / "flare_reference.fine.ref.panel"
BROAD_MAP = OUT / "fine_to_broad.tsv"


def samples(vcf):
    x = subprocess.check_output(
        ["bcftools", "query", "-l", str(vcf)],
        text=True,
    )
    return [s.strip() for s in x.splitlines() if s.strip()]


def safe_label(s):
    s = s.strip()
    s = re.sub(r"[^A-Za-z0-9]+", "_", s)
    return s.strip("_")


# ------------------------------------------------------------
# Original -> sanitized sample mapping
# ------------------------------------------------------------

merged_original = samples(MERGED_VCF)

sanitized = [
    x.strip()
    for x in SANITIZED_NAMES.read_text().splitlines()
    if x.strip()
]

if len(merged_original) != len(sanitized):
    raise RuntimeError("Sample count mismatch")

rename = dict(zip(merged_original, sanitized))


# ------------------------------------------------------------
# Fine assignments
# ------------------------------------------------------------

fine = {}
broad = {}


def add(sample, population, broad_group, source):
    label = f"{source}_{safe_label(population)}"

    if sample in fine:
        raise RuntimeError(f"Duplicate sample: {sample}")

    fine[sample] = label
    broad[label] = broad_group


# ------------------------------------------------------------
# HGDP + 1KG
#
# Expected original panel format:
# sample population [possibly extra columns]
# ------------------------------------------------------------

hgdp_samples = set(samples(HGDP_VCF))

with open(HGDP_FINE_PANEL) as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        sample = fields[0]
        population = fields[1]

        if sample not in hgdp_samples:
            continue

        # Keep original fine population.
        # Broad classification is intentionally not inferred here yet.
        add(
            sample,
            population,
            "UNASSIGNED",
            "HGDP1KG"
        )


# ------------------------------------------------------------
# GSE73996
#
# Try to locate sample/population/region columns by header.
# ------------------------------------------------------------

gse73996_samples = set(samples(GSE73996_VCF))

with open(GSE73996_META) as f:
    reader = csv.DictReader(f, delimiter="\t")

    cols = reader.fieldnames

    sample_col = next(
        x for x in cols
        if x.lower() in {"sample", "sample_id", "sampleid", "id", "title"}
    )

    pop_col = next(
        x for x in cols
        if x.lower() in {"population", "pop"}
    )

    region_col = next(
        (
            x for x in cols
            if x.lower() in {"region", "geographic_region"}
        ),
        None
    )

    for row in reader:

        sample = row[sample_col].strip()

        if sample not in gse73996_samples:
            continue

        population = row[pop_col].strip()

        broad_group = (
            row[region_col].strip()
            if region_col
            else "UNASSIGNED"
        )

        add(
            sample,
            population,
            broad_group,
            "GSE73996"
        )


# ------------------------------------------------------------
# GSE77508
#
# Expected metadata columns:
# sample, population, region
# ------------------------------------------------------------

gse77508_samples = set(samples(GSE77508_VCF))

with open(GSE77508_META) as f:
    reader = csv.DictReader(f, delimiter="\t")

    cols = reader.fieldnames

    sample_col = next(
        x for x in cols
        if x.lower() in {"sample", "sample_id", "sampleid", "id", "title"}
    )

    pop_col = next(
        x for x in cols
        if x.lower() in {"population", "pop"}
    )

    for row in reader:

        sample = row[sample_col].strip()

        if sample not in gse77508_samples:
            continue

        population = row[pop_col].strip()

        add(
            sample,
            population,
            "SEA",
            "GSE77508"
        )


# ------------------------------------------------------------
# Kazakh 2025
# ------------------------------------------------------------

for sample in samples(KAZAKH_VCF):

    add(
        sample,
        "Kazakh",
        "CAS",
        "Kazakh2025"
    )


# ------------------------------------------------------------
# Validate all merged samples have assignment
# ------------------------------------------------------------

missing = [
    x for x in merged_original
    if x not in fine
]

if missing:
    raise RuntimeError(
        "Missing fine labels:\n"
        + "\n".join(missing[:100])
    )


# ------------------------------------------------------------
# Write FLARE2 panel with sanitized IDs
# ------------------------------------------------------------

with open(FINE_PANEL, "w") as f:

    for original in merged_original:

        f.write(
            f"{rename[original]}\t{fine[original]}\n"
        )


# ------------------------------------------------------------
# Write fine -> broad lookup
# ------------------------------------------------------------

with open(BROAD_MAP, "w") as f:

    f.write("fine_population\tbroad_group\n")

    for label in sorted(broad):
        f.write(
            f"{label}\t{broad[label]}\n"
        )


print(f"Fine panel: {FINE_PANEL}")
print(f"Broad map:  {BROAD_MAP}")
print(f"Samples:    {len(fine)}")
print(f"Panels:     {len(set(fine.values()))}")
