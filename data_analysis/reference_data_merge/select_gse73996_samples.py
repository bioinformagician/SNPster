#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

INPUT = Path(
    "/srv/dependencies/ancestry/gse73996/GSE73996_metadata.tsv"
)

OUTDIR = Path(
    "/srv/dependencies/ancestry/gse73996"
)

# Mapping into the fixed SNPster / PGS-Catalog-aligned groups.
REGION_TO_SUPERPOP = {

    # European
    "West Europe": "EUR",
    "East Europe": "EUR",

    # East Asian
    "East Asia": "EAS",
    "Mongolia": "EAS",

    # South Asian
    "South Asia": "SAS",

    # Southeast Asian
    "South East Asia": "SEA",

    # Oceanian
    "Oceania": "OCE",

    # Sub-Saharan African
    "Sub Saharan Africa": "SSA",

    # Greater Middle Eastern
    "Near East": "MID",
    "North Africa": "MID",
    "Caucasus": "MID",

    # Central Asian
    "Central Asia": "CAS",
}


selected = []
excluded = []


with open(INPUT, newline="") as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t",
    )

    for row in reader:

        region = row["region"].strip()

        super_pop = REGION_TO_SUPERPOP.get(region)

        if super_pop:

            row["super_pop"] = super_pop
            selected.append(row)

        else:

            row["super_pop"] = "EXCLUDED"
            excluded.append(row)


# ------------------------------------------------------
# Full selected metadata
# ------------------------------------------------------

selected_metadata = OUTDIR / "GSE73996_flare_metadata.tsv"

fields = [
    "sample",
    "population",
    "region",
    "super_pop",
    "gsm",
    "gender",
    "title",
]

with open(selected_metadata, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(selected)


# ------------------------------------------------------
# bcftools sample inclusion file
# ------------------------------------------------------

sample_file = OUTDIR / "gse73996_flare_samples.txt"

with open(sample_file, "w") as fh:

    for row in selected:
        fh.write(row["sample"] + "\n")


# ------------------------------------------------------
# FLARE 2-column reference panel
# ------------------------------------------------------

panel_file = OUTDIR / "gse73996_flare_panel.txt"

with open(panel_file, "w") as fh:

    for row in selected:
        fh.write(
            f"{row['sample']}\t{row['super_pop']}\n"
        )


# ------------------------------------------------------
# Excluded individuals, retained for audit trail
# ------------------------------------------------------

excluded_file = OUTDIR / "GSE73996_excluded_metadata.tsv"

with open(excluded_file, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(excluded)


# ------------------------------------------------------
# Summary
# ------------------------------------------------------

counts = Counter(
    row["super_pop"]
    for row in selected
)

print()
print("Selected reference individuals")
print("------------------------------")

for pop in [
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
]:
    print(f"{pop:5s}: {counts.get(pop, 0):4d}")

print("------------------------------")
print(f"Selected : {len(selected)}")
print(f"Excluded : {len(excluded)}")
print(f"Total    : {len(selected) + len(excluded)}")

print()
print(f"Metadata : {selected_metadata}")
print(f"Samples  : {sample_file}")
print(f"Panel    : {panel_file}")
print(f"Excluded : {excluded_file}")