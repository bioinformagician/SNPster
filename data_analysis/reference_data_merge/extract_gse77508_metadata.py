#!/usr/bin/env python3

import gzip
import re
import subprocess
from pathlib import Path

SOFT = Path("GSE77508_family.soft.gz")
VCF = Path("GSE77508.chr1.phased.vcf.gz")
OUT = Path("GSE77508_metadata.tsv")


def vcf_samples():
    out = subprocess.check_output(
        ["bcftools", "query", "-l", str(VCF)],
        text=True
    )
    return [x.strip() for x in out.splitlines() if x.strip()]


records = []
current = None

with gzip.open(SOFT, "rt", errors="replace") as f:
    for raw in f:
        line = raw.rstrip("\n")

        if line.startswith("^SAMPLE = "):

            if current:
                records.append(current)

            current = {
                "gsm": line.split("=", 1)[1].strip(),
                "title": "",
                "population": "",
                "region": "",
                "gender": "",
                "description": "",
            }

        elif current is not None:

            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()

            elif line.startswith("!Sample_description = "):
                value = line.split("=", 1)[1].strip()

                if current["description"]:
                    current["description"] += "; " + value
                else:
                    current["description"] = value

            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()

                if ":" not in value:
                    continue

                key, val = value.split(":", 1)
                key = key.strip().lower()
                val = val.strip()

                if key == "population":
                    current["population"] = val

                elif key == "region":
                    current["region"] = val

                elif key == "gender":
                    current["gender"] = val

if current:
    records.append(current)


# Map GEO titles to the sample names actually used in the VCF.
#
# The processed matrix/VCF names can differ slightly from GEO titles,
# especially because the original matrix columns had ".GType".
def norm(x):
    x = x.strip()
    x = re.sub(r"\.GType$", "", x, flags=re.I)
    x = re.sub(r"[^A-Za-z0-9]+", "", x)
    return x.lower()


by_title = {}

for r in records:
    key = norm(r["title"])

    if key in by_title:
        raise RuntimeError(
            f"Duplicate normalized GEO title: {r['title']}"
        )

    by_title[key] = r


samples = vcf_samples()

matched = []
missing = []

for sample in samples:

    key = norm(sample)

    if key not in by_title:
        missing.append(sample)
        continue

    r = by_title[key]

    matched.append({
        "sample": sample,
        "gsm": r["gsm"],
        "geo_title": r["title"],
        "population": r["population"],
        "region": r["region"],
        "gender": r["gender"],
        "description": r["description"],
    })


if missing:
    print("VCF samples not matched to GEO metadata:")
    for x in missing:
        print(" ", x)

    raise RuntimeError(
        f"{len(missing)} VCF samples could not be mapped"
    )


with open(OUT, "w") as f:

    f.write(
        "sample\tgsm\tgeo_title\tpopulation\tregion\tgender\tdescription\n"
    )

    for r in matched:
        f.write(
            "\t".join([
                r["sample"],
                r["gsm"],
                r["geo_title"],
                r["population"],
                r["region"],
                r["gender"],
                r["description"],
            ]) + "\n"
        )


print(f"GEO samples parsed : {len(records)}")
print(f"VCF samples        : {len(samples)}")
print(f"Matched            : {len(matched)}")
print(f"Output             : {OUT}")
