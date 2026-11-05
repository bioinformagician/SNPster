#!/usr/bin/env python3

import gzip
import re
import csv
from collections import Counter

src = "/srv/dependencies/ancestry/gse73996/GSE73996_family.soft.gz"
out = "/srv/dependencies/ancestry/gse73996/GSE73996_metadata.tsv"

samples = []
cur = None


def finish(x):
    if not x:
        return

    title = x.get("title", "")

    m = re.search(r"\[([^\]]+)\]", title)

    if m:
        x["sample"] = m.group(1)

    if x.get("sample"):
        samples.append(x)


with gzip.open(src, "rt") as fh:

    for raw in fh:

        line = raw.rstrip("\n")

        if line.startswith("^SAMPLE = "):

            finish(cur)

            cur = {
                "gsm": line.split("=", 1)[1].strip(),
                "sample": "",
                "title": "",
                "population": "",
                "region": "",
                "gender": "",
            }

        elif cur is not None:

            if line.startswith("!Sample_title = "):
                cur["title"] = line.split("=", 1)[1].strip()

            elif line.startswith("!Sample_characteristics"):

                value = line.split("=", 1)[1].strip()

                if ":" not in value:
                    continue

                key, value = value.split(":", 1)

                key = key.strip().lower()
                value = value.strip()

                if key == "population":
                    cur["population"] = value

                elif key == "region":
                    cur["region"] = value

                elif key == "gender":
                    cur["gender"] = value

finish(cur)


with open(out, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "sample",
            "gsm",
            "population",
            "region",
            "gender",
            "title",
        ],
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(samples)


print(f"Wrote {len(samples)} samples to {out}")

print("\nPopulation counts:")
for pop, n in sorted(Counter(x["population"] for x in samples).items()):
    print(f"{pop:30s} {n}")

print("\nRegion counts:")
for region, n in sorted(Counter(x["region"] for x in samples).items()):
    print(f"{region:30s} {n}")