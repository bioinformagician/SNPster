#!/usr/bin/env python3

import subprocess
from pathlib import Path
from collections import Counter

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

HGDP_VCF = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/"
    "hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5."
    "SNV_biallelic.numericCHR.vcf.gz"
)

HGDP_PANEL = Path(
    "/srv/dependencies/ancestry/hgdp_1kg_references/"
    "hgdp_1kg_flare_panel.txt"
)

GSE77508_VCF = Path(
    "/srv/dependencies/ancestry/gse77508/"
    "GSE77508.chr1.phased.vcf.gz"
)

GSE73996_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/gse73996/"
    "selected/phased/gse73996.chr1.phased.vcf.gz"
)

GSE73996_PANEL = Path(
    "/srv/dependencies/ancestry/gse73996/"
    "gse73996_flare_panel.txt"
)

KAZAKH_VCF = Path(
    "/srv/dependencies/ancestry/harmonized/kazakh/"
    "vcf/phased/kazakh.chr1.phased.vcf.gz"
)

MERGED_VCF = Path(
    "/srv/dependencies/ancestry/flare_reference_v1/output/"
    "flare_reference.chr1.vcf.gz"
)

OUTPUT = Path(
    "/srv/dependencies/ancestry/flare_reference_v1/output/"
    "flare_reference.ref.panel"
)

# Desired FLARE ancestry order
ANCESTRY_ORDER = [
    "EUR",
    "EAS",
    "SAS",
    "SEA",
    "OCE",
    "SSA",
    "MID",
    "CAS",
    "NAM",
    "AMR",
]

VALID = set(ANCESTRY_ORDER)


def samples(vcf):
    out = subprocess.check_output(
        ["bcftools", "query", "-l", str(vcf)],
        text=True,
    )
    return [x.strip() for x in out.splitlines() if x.strip()]


def read_panel(path):
    mapping = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 2:
                continue

            sample, ancestry = fields[0], fields[1]

            # Also safely skips accidental headers
            if ancestry not in VALID:
                continue

            if sample in mapping:
                raise RuntimeError(
                    f"Duplicate sample in {path}: {sample}"
                )

            mapping[sample] = ancestry

    return mapping


# ------------------------------------------------------------
# Read actual samples from the VCFs
# ------------------------------------------------------------

hgdp_samples = samples(HGDP_VCF)
gse77508_samples = samples(GSE77508_VCF)
gse73996_samples = samples(GSE73996_VCF)
kazakh_samples = samples(KAZAKH_VCF)
merged_samples = samples(MERGED_VCF)

hgdp_panel = read_panel(HGDP_PANEL)
gse73996_panel = read_panel(GSE73996_PANEL)


# ------------------------------------------------------------
# Construct sample -> ancestry mapping
# ------------------------------------------------------------

mapping = {}


def add(sample, ancestry, source):
    if ancestry not in VALID:
        raise RuntimeError(
            f"Invalid ancestry '{ancestry}' for {sample}"
        )

    if sample in mapping:
        raise RuntimeError(
            f"Duplicate sample across sources: {sample} ({source})"
        )

    mapping[sample] = ancestry


# HGDP+1kGP
for sample in hgdp_samples:
    if sample not in hgdp_panel:
        raise RuntimeError(
            f"HGDP sample has no panel assignment: {sample}"
        )

    add(sample, hgdp_panel[sample], "HGDP+1kGP")


# GSE77508
# GEO VCF samples are the Southeast Asian component.
for sample in gse77508_samples:
    add(sample, "SEA", "GSE77508")


# GSE73996
for sample in gse73996_samples:
    if sample not in gse73996_panel:
        raise RuntimeError(
            f"GSE73996 sample has no panel assignment: {sample}"
        )

    add(sample, gse73996_panel[sample], "GSE73996")


# Kazakh
for sample in kazakh_samples:
    add(sample, "CAS", "Kazakh")


# ------------------------------------------------------------
# Validate against merged VCF
# ------------------------------------------------------------

merged_set = set(merged_samples)
mapped_set = set(mapping)

missing = merged_set - mapped_set
extra = mapped_set - merged_set

if missing:
    raise RuntimeError(
        "Samples in merged VCF but missing from panel:\n"
        + "\n".join(sorted(missing))
    )

if extra:
    raise RuntimeError(
        "Samples in panel but absent from merged VCF:\n"
        + "\n".join(sorted(extra))
    )

if len(mapping) != len(merged_samples):
    raise RuntimeError(
        f"Panel has {len(mapping)} samples but merged VCF "
        f"has {len(merged_samples)}"
    )


# ------------------------------------------------------------
# Write in explicit ancestry order
#
# This controls FLARE's ancestry index order.
# ------------------------------------------------------------

with open(OUTPUT, "w") as out:
    for ancestry in ANCESTRY_ORDER:
        for sample in merged_samples:
            if mapping[sample] == ancestry:
                out.write(f"{sample}\t{ancestry}\n")


# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

counts = Counter(mapping.values())

print()
print("Reference panel created:")
print(OUTPUT)
print()

print("Source sample counts:")
print(f"  HGDP+1kGP : {len(hgdp_samples)}")
print(f"  GSE77508  : {len(gse77508_samples)}")
print(f"  GSE73996  : {len(gse73996_samples)}")
print(f"  Kazakh    : {len(kazakh_samples)}")
print(f"  TOTAL     : {len(mapping)}")
print()

print("FLARE ancestry counts:")

for ancestry in ANCESTRY_ORDER:
    print(f"  {ancestry:4s} {counts[ancestry]:5d}")

print()
print("Validation: OK")
