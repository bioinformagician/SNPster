#!/usr/bin/env bash
set -euo pipefail

BASE=/output/source_model
REF=${BASE}/source_reference.autosomes.vcf.gz
PANEL=${BASE}/source_reference.ref.panel
MAP=/output/model_sweep_allchr/plink.autosomes.GRCh38.map
FLARE_MODEL=${BASE}/models/diag_K18/train.final.model
OUTROOT=${BASE}/real_user_validation

mkdir -p "$OUTROOT"

# name|input_directory|IMPID prefix
CASES=(
  "frederik|/app/imputed_data/frederik_all_re/330/qc_output|IMPID330"
  "lonni|/app/imputed_data/Lonni/325/qc_output|IMPID325"
  "murtaza|/app/imputed_data/murtaza/322/qc_output|IMPID322"
)

for CASE in "${CASES[@]}"; do
  IFS='|' read -r NAME INDIR PREFIX <<< "$CASE"

  OUTDIR="${OUTROOT}/${NAME}"
  mkdir -p "$OUTDIR"

  echo
  echo "======================================================================"
  echo "Preparing ${NAME}"
  echo "======================================================================"

  FILES=()

  for CHR in $(seq 1 22); do

    F="${INDIR}/${PREFIX}.chr${CHR}.split_QC.vcf.gz"

    if [[ ! -s "$F" ]]; then
      echo "ERROR: missing $F" >&2
      exit 1
    fi

    FILES+=("$F")

  done

  TARGET="${OUTDIR}/${NAME}.autosomes.vcf.gz"

  bcftools concat \
    --threads 20 \
    -a \
    -Oz \
    -o "$TARGET" \
    "${FILES[@]}"

  tabix -f -p vcf "$TARGET"

  NSAMPLES=$(bcftools query -l "$TARGET" | wc -l)

  if [[ "$NSAMPLES" -ne 1 ]]; then
    echo "ERROR: expected exactly one sample in $TARGET, found $NSAMPLES" >&2
    exit 1
  fi

  echo "Target sample: $(bcftools query -l "$TARGET")"

  echo
  echo "Running frozen K18 FLARE model for ${NAME}"

  flare \
    gt="$TARGET" \
    ref="$REF" \
    ref-panel="$PANEL" \
    map="$MAP" \
    model="$FLARE_MODEL" \
    update-p=false \
    em=false \
    probs=false \
    array=false \
    nthreads=20 \
    seed=20261002 \
    out="${OUTDIR}/K18"

  if [[ ! -s "${OUTDIR}/K18.global.anc.gz" ]]; then
    echo "ERROR: FLARE did not create ${OUTDIR}/K18.global.anc.gz" >&2
    exit 1
  fi

done

python3 <<'PY'
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


BASE = Path("/output/source_model")
OUTROOT = BASE / "real_user_validation"

SOURCES = [
    "EUR",
    "EAS",
    "SAS",
    "SEA",
    "GME",
    "NAM",
    "OCE",
    "SSA",
]

USERS = [
    "frederik",
    "lonni",
    "murtaza",
]


base_model = joblib.load(
    BASE
    / "crossfit_calibration"
    / "k18_extratrees_base.joblib"
)

calibrator = joblib.load(
    BASE
    / "crossfit_calibration"
    / "k18_linear_calibrator.joblib"
)


def simplex(v):

    v = np.asarray(
        v,
        dtype=float,
    )

    u = np.sort(v)[::-1]

    cssv = (
        np.cumsum(u)
        - 1.0
    )

    ind = np.arange(
        1,
        len(v) + 1,
    )

    cond = (
        u
        - cssv / ind
    ) > 0

    if not np.any(cond):
        return (
            np.ones_like(v)
            / len(v)
        )

    rho = ind[cond][-1]

    theta = (
        cssv[rho - 1]
        / rho
    )

    w = np.maximum(
        v - theta,
        0.0,
    )

    if w.sum() <= 0:
        return (
            np.ones_like(v)
            / len(v)
        )

    return (
        w
        / w.sum()
    )


def simplex_matrix(x):

    return np.vstack([
        simplex(row)
        for row in np.asarray(x)
    ])


rows = []


for user in USERS:

    path = (
        OUTROOT
        / user
        / "K18.global.anc.gz"
    )

    df = pd.read_csv(
        path,
        sep=r"\s+",
        compression="gzip",
    )

    anc_cols = [
        c
        for c in df.columns
        if c.startswith("anc_")
    ]

    if len(df) != 1:
        raise RuntimeError(
            f"{user}: expected one ancestry row, got {len(df)}"
        )

    if len(anc_cols) != 18:
        raise RuntimeError(
            f"{user}: expected 18 latent columns, got {len(anc_cols)}"
        )

    x = df[
        anc_cols
    ].to_numpy(
        dtype=float
    )

    x = (
        x
        / x.sum(
            axis=1,
            keepdims=True,
        )
    )

    raw = simplex_matrix(
        base_model.predict(x)
    )[0]

    cal = simplex_matrix(
        calibrator.predict(
            raw.reshape(1, -1)
        )
    )[0]

    sample_col = (
        "SAMPLE"
        if "SAMPLE" in df.columns
        else df.columns[0]
    )

    row = {
        "user": user,
        "sample": str(
            df.iloc[0][sample_col]
        ),
    }

    for source, value in zip(
        SOURCES,
        raw,
    ):

        row[
            f"raw_{source}"
        ] = float(value)

    for source, value in zip(
        SOURCES,
        cal,
    ):

        row[
            f"cal_{source}"
        ] = float(value)

    rows.append(row)


out = pd.DataFrame(rows)

outfile = (
    OUTROOT
    / "three_user_frozen_ancestry.tsv"
)

out.to_csv(
    outfile,
    sep="\t",
    index=False,
)


print()
print("=" * 92)
print(
    "FROZEN SNPSTER ANCESTRY — "
    "CALIBRATED PRODUCTION OUTPUT"
)
print("=" * 92)


for _, row in out.iterrows():

    print(
        f"\n{row['user'].upper()} "
        f"({row['sample']})"
    )

    vals = [
        (
            source,
            row[f"cal_{source}"],
        )
        for source in SOURCES
    ]

    vals.sort(
        key=lambda x: x[1],
        reverse=True,
    )

    for source, value in vals:

        print(
            f"  {source:3s}  "
            f"{100 * value:7.2f}%"
        )


print()
print("=" * 92)
print(
    "RAW DECODER OUTPUT — DIAGNOSTIC ONLY"
)
print("=" * 92)


for _, row in out.iterrows():

    print(
        f"\n{row['user'].upper()} "
        f"({row['sample']})"
    )

    vals = [
        (
            source,
            row[f"raw_{source}"],
        )
        for source in SOURCES
    ]

    vals.sort(
        key=lambda x: x[1],
        reverse=True,
    )

    for source, value in vals:

        print(
            f"  {source:3s}  "
            f"{100 * value:7.2f}%"
        )


print(
    f"\nSaved: {outfile}"
)
PY
