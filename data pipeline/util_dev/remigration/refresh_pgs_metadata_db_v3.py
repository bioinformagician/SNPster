import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

# If you already set excel_filepath in an earlier cell, this keeps it.
# Otherwise it falls back to the repo copy of the metadata file.
excel_filepath = globals().get(
    "excel_filepath",
    r"/home/frederik/github_projects/SNPster/data pipeline/util_dev/remigration/pgs_all_metadata.xlsx"
)

# Optional normalized metric seed used only during ETL/database refresh.
# The report itself will read all curated metadata from PostgreSQL.
metric_long_filepath = Path(
    os.getenv(
        "PGS_METRIC_LONG_FILE",
        str(Path(excel_filepath).resolve().parent / "pgs_performance_metrics_long.csv"),
    )
)


def _first_existing(series_or_df, candidates, required=True, default=None):
    cols = list(series_or_df.columns)

    def _norm(name):
        return " ".join(str(name).strip().lower().split())

    normalized_to_original = {_norm(col): col for col in cols}

    for candidate in candidates:
        if candidate in cols:
            return series_or_df[candidate]

        candidate_norm = _norm(candidate)
        if candidate_norm in normalized_to_original:
            return series_or_df[normalized_to_original[candidate_norm]]

    if required:
        raise KeyError(f"Missing required column. Expected one of: {candidates}")
    return pd.Series([default] * len(series_or_df), index=series_or_df.index)


def _to_nullable_int(series):
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def _to_bool(series):
    mapping = {
        "true": True,
        "false": False,
        "yes": True,
        "no": False,
        "y": True,
        "n": False,
        "1": True,
        "0": False,
    }
    normalized = series.astype(str).str.strip().str.lower()
    return normalized.map(mapping).where(~series.isna(), None)


def _parse_effect_with_ci(series):
    """
    Parse values like:
      1.55 [1.52,1.58]
      1.55 (1.52-1.58)
      1.55
    Returns three numeric series: point_estimate, ci_lower, ci_upper.
    """
    s = series.astype(str)

    # First numeric token is the point estimate.
    point = pd.to_numeric(
        s.str.extract(r"([+-]?(?:\d+(?:\.\d+)?|\.\d+))", expand=False),
        errors="coerce",
    )

    # Capture lower/upper inside [] or () using either comma or dash separator.
    ci_match = s.str.extract(
        r"[\[\(]\s*([+-]?(?:\d+(?:\.\d+)?|\.\d+))\s*[,\-]\s*([+-]?(?:\d+(?:\.\d+)?|\.\d+))\s*[\]\)]",
        expand=True,
    )
    ci_lower = pd.to_numeric(ci_match[0], errors="coerce")
    ci_upper = pd.to_numeric(ci_match[1], errors="coerce")

    return point, ci_lower, ci_upper


def _parse_point_estimate(series):
    """Parse the leading numeric point estimate from strings with optional CI text."""
    point, _, _ = _parse_effect_with_ci(series)
    return point


def _split_pgs_ids(value):
    """
    Split one or many PGS IDs from metadata cells like:
      PGS000001
      PGS000001, PGS000004, PGS000007
      PGS000001;PGS000004
    """
    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text:
        return []

    parts = [p.strip() for p in text.replace(";", ",").split(",")]
    return [p for p in parts if p]


def _clean_df_for_sql(df):
    cleaned = df.copy()
    cleaned = cleaned.where(pd.notna(cleaned), None)
    return cleaned


def _to_db_scalar(value):
    """Convert pandas/numpy scalar types to psycopg2-compatible Python scalars."""
    if value is None:
        return None

    # Handle pandas missing sentinels (pd.NA, NaN, NaT) as SQL NULL.
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass

    if isinstance(value, (pd.Timestamp, np.datetime64)):
        ts = pd.Timestamp(value)
        return ts.to_pydatetime()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _get_table_columns(cursor, schema, table):
    cursor.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s AND table_name = %s
        ORDER BY ordinal_position
        """,
        (schema, table),
    )
    return [row[0] for row in cursor.fetchall()]


def _bulk_insert(cursor, schema, table, dataframe):
    if dataframe.empty:
        print(f"Skipping {schema}.{table}: no rows")
        return

    table_columns = _get_table_columns(cursor, schema, table)
    columns = [c for c in dataframe.columns if c in table_columns]
    skipped = [c for c in dataframe.columns if c not in table_columns]

    if skipped:
        print(f"Skipping non-existent columns for {schema}.{table}: {skipped}")
    if not columns:
        raise ValueError(f"No matching dataframe columns for {schema}.{table}")

    tuples = [
        tuple(_to_db_scalar(v) for v in row)
        for row in dataframe[columns].itertuples(index=False, name=None)
    ]
    query = f"INSERT INTO {schema}.{table} ({', '.join(columns)}) VALUES %s"
    execute_values(cursor, query, tuples, page_size=1000)
    print(f"Inserted {len(tuples)} rows into {schema}.{table}")




def _table_exists(cursor, schema, table):
    cursor.execute("SELECT to_regclass(%s)", (f"{schema}.{table}",))
    return cursor.fetchone()[0] is not None


def _bulk_upsert(cursor, schema, table, dataframe, conflict_columns):
    """Bulk upsert parent/reference metadata without truncating FK dependants."""
    if dataframe.empty:
        print(f"Skipping {schema}.{table}: no rows")
        return

    table_columns = _get_table_columns(cursor, schema, table)
    columns = [c for c in dataframe.columns if c in table_columns]
    skipped = [c for c in dataframe.columns if c not in table_columns]

    if skipped:
        print(f"Skipping non-existent columns for {schema}.{table}: {skipped}")

    missing_conflict = [c for c in conflict_columns if c not in columns]
    if missing_conflict:
        raise ValueError(
            f"Conflict columns missing for {schema}.{table}: {missing_conflict}"
        )

    tuples = [
        tuple(_to_db_scalar(v) for v in row)
        for row in dataframe[columns].itertuples(index=False, name=None)
    ]

    update_columns = [c for c in columns if c not in conflict_columns]
    if update_columns:
        update_sql = ", ".join(
            f"{c} = EXCLUDED.{c}" for c in update_columns
        )
        conflict_sql = (
            f"ON CONFLICT ({', '.join(conflict_columns)}) "
            f"DO UPDATE SET {update_sql}"
        )
    else:
        conflict_sql = (
            f"ON CONFLICT ({', '.join(conflict_columns)}) DO NOTHING"
        )

    query = (
        f"INSERT INTO {schema}.{table} ({', '.join(columns)}) VALUES %s "
        f"{conflict_sql}"
    )
    execute_values(cursor, query, tuples, page_size=1000)
    print(f"Upserted {len(tuples)} rows into {schema}.{table}")



_ALLOWED_EFFECT_BASIS = {
    "per_sd",
    "per_sd_assumed",
    "group_contrast",
    "absolute",
    "unknown",
    "not_applicable",
}


def _normalize_metric_source(value):
    if value is None:
        return None

    source = str(value).strip().casefold()

    return {
        "dedicated_column": "dedicated",
        "dedicated": "dedicated",
        "other_metric": "other_metric",
    }.get(source)


def _normalize_effect_basis(
    metric_type,
    source,
    raw_basis,
    display_mode,
    personalization_allowed,
):
    """
    Convert parser/free-text semantics to the DB controlled vocabulary.

    Long comparison descriptions must never be stored in effect_basis;
    they are preserved separately in contrast_text/raw_metric.
    """
    metric_type = "" if metric_type is None else str(metric_type).strip()
    source = "" if source is None else str(source).strip()
    display_mode = "" if display_mode is None else str(display_mode).strip()

    if raw_basis is None:
        raw_basis_text = ""
    else:
        try:
            raw_basis_text = (
                ""
                if pd.isna(raw_basis)
                else str(raw_basis).strip()
            )
        except TypeError:
            raw_basis_text = str(raw_basis).strip()

    if raw_basis_text in _ALLOWED_EFFECT_BASIS:
        return raw_basis_text

    # Dedicated PGS Catalog OR/HR values are treated as per-SD unless
    # the curated metadata explicitly disabled personalization.
    if source == "dedicated" and metric_type in {
        "odds_ratio",
        "hazard_ratio",
    }:
        return (
            "per_sd_assumed"
            if personalization_allowed
            else "group_contrast"
        )

    if source == "dedicated" and metric_type == "beta":
        return "per_sd"

    if (
        display_mode == "group_comparison"
        or metric_type in {
            "group_odds_ratio",
            "group_hazard_ratio",
            "relative_risk",
            "relative_risk_ratio",
            "incidence_rate_ratio",
        }
    ):
        return "group_contrast"

    return "not_applicable"


def _metric_contrast_text(
    metric_type,
    raw_basis,
    metric_label_raw,
    display_mode,
):
    """
    Preserve human-readable group comparison wording without placing it
    in effect_basis.
    """
    metric_type = "" if metric_type is None else str(metric_type).strip()
    display_mode = "" if display_mode is None else str(display_mode).strip()

    is_group_comparison = (
        display_mode == "group_comparison"
        or metric_type in {
            "group_odds_ratio",
            "group_hazard_ratio",
            "relative_risk",
            "relative_risk_ratio",
            "incidence_rate_ratio",
        }
    )

    if not is_group_comparison:
        return None

    for value in (raw_basis, metric_label_raw):
        if value is None:
            continue

        try:
            if pd.isna(value):
                continue
        except TypeError:
            pass

        value = str(value).strip()

        if value:
            return value

    return None



_BETA_UNIT_TOKEN = (
    r"(?:"
    r"kg/m(?:\^?2|²)|"
    r"mmol/mol|mmol/L|mg/dL|mg/L|g/dL|g/L|"
    r"[µμu]mol/L|nmol/L|pmol/L|"
    r"ng/mL|pg/mL|[µμu]g/L|"
    r"mL/min/1\.73m(?:\^?2|²)|"
    r"mmHg|beats/min|U/L|IU/L|"
    r"mm(?:\^?2|²)?|cm(?:\^?2|²)?|"
    r"years?|days?|%"
    r")"
)


def _canonical_beta_unit(unit):
    if unit is None:
        return None

    unit = str(unit).strip().rstrip(".,;:)")

    replacements = {
        "kg/m^2": "kg/m²",
        "kg/m2": "kg/m²",
        "mm^2": "mm²",
        "mm2": "mm²",
        "cm^2": "cm²",
        "cm2": "cm²",
        "mL/min/1.73m^2": "mL/min/1.73m²",
        "mL/min/1.73m2": "mL/min/1.73m²",
        "umol/L": "µmol/L",
        "μmol/L": "µmol/L",
        "ug/L": "µg/L",
        "μg/L": "µg/L",
        "year": "years",
        "day": "days",
    }

    return replacements.get(unit, unit)


def _extract_beta_unit_from_performance_info(value):
    """
    Extract a beta unit only when the PGS Catalog explicitly describes it.
    """
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass

    text = str(value).strip()
    if not text:
        return None

    patterns = [
        rf"beta\s+(?:is\s+)?(?:in\s+)?units?\s*(?:of|=|is)?\s*({_BETA_UNIT_TOKEN})",
        rf"beta\s+unit\s*(?:is|=)?\s*({_BETA_UNIT_TOKEN})",
        rf"unit\s+of\s+measurement\s*=\s*({_BETA_UNIT_TOKEN})",
        rf"(?:increase|decrease)\s+in\s+({_BETA_UNIT_TOKEN})\s+per\s+PGS\s+SD",
    ]

    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return _canonical_beta_unit(match.group(1))

    return None


def _extract_beta_unit_from_phenotype_definition(value):
    """
    Conservative fallback using the linked PSS phenotype definition.

    No unit is inferred from the trait name. Parenthetical short units are
    case-sensitive so '(CM)' in 'ICD-10-CM' cannot become centimetres.
    """
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass

    text = str(value).strip()
    if not text:
        return None

    patterns = [
        rf"(?:represented|measured|reported|expressed)\s+in\s+({_BETA_UNIT_TOKEN})",
        rf"units?\s+(?:of|are|is)\s+({_BETA_UNIT_TOKEN})",
    ]

    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return _canonical_beta_unit(match.group(1))

    parenthetical_units = {
        "kg/m^2", "kg/m2", "kg/m²",
        "mmol/mol", "mmol/L", "mg/dL", "mg/L", "g/dL", "g/L",
        "µmol/L", "μmol/L", "umol/L", "nmol/L", "pmol/L",
        "ng/mL", "pg/mL", "µg/L", "μg/L", "ug/L",
        "mL/min/1.73m^2", "mL/min/1.73m2", "mL/min/1.73m²",
        "mmHg", "beats/min", "U/L", "IU/L",
        "mm", "mm2", "mm^2", "mm²",
        "cm", "cm2", "cm^2", "cm²",
        "years", "year", "days", "day", "%",
    }

    for candidate in re.findall(r"\(([^()]*)\)", text):
        candidate = candidate.strip()
        if candidate in parenthetical_units:
            return _canonical_beta_unit(candidate)

    return None


def _build_beta_unit_by_ppm(performance_raw, eval_sets_raw):
    """
    PPM -> beta unit.

    Priority:
      1) explicit unit in Performance Metrics 'Other Relevant Information'
      2) explicit unit in the linked evaluation-sample phenotype definition
    """
    ppm_col = _first_existing(
        performance_raw,
        ["PGS Performance Metric (PPM) ID"],
    )
    pss_col = _first_existing(
        performance_raw,
        ["PGS Sample Set (PSS)"],
        required=False,
    )
    info_col = _first_existing(
        performance_raw,
        ["PGS Performance: Other Relevant Information"],
        required=False,
    )

    pss_id_col = _first_existing(
        eval_sets_raw,
        ["PGS Sample Set (PSS)"],
    )
    phenotype_col = _first_existing(
        eval_sets_raw,
        ["Phenotype Definitions and Methods"],
        required=False,
    )

    phenotype_by_pss = {}
    for pss_id, phenotype in zip(pss_id_col, phenotype_col):
        if pss_id is None or pd.isna(pss_id):
            continue
        key = str(pss_id).strip()
        if key and key not in phenotype_by_pss:
            phenotype_by_pss[key] = phenotype

    result = {}
    n_info = 0
    n_phenotype = 0

    for ppm_id, pss_id, info in zip(ppm_col, pss_col, info_col):
        if ppm_id is None or pd.isna(ppm_id):
            continue

        ppm_key = str(ppm_id).strip()
        if not ppm_key:
            continue

        unit = _extract_beta_unit_from_performance_info(info)
        if unit:
            result[ppm_key] = unit
            n_info += 1
            continue

        if pss_id is None or pd.isna(pss_id):
            continue

        phenotype = phenotype_by_pss.get(str(pss_id).strip())
        unit = _extract_beta_unit_from_phenotype_definition(phenotype)

        if unit:
            result[ppm_key] = unit
            n_phenotype += 1

    print(
        f"Beta units detected for {len(result)} PPMs "
        f"({n_info} from performance info, "
        f"{n_phenotype} from phenotype definitions)."
    )

    return result


def _load_normalized_metric_values(
    cursor,
    csv_path: Path,
    beta_unit_by_ppm=None,
):
    """Refresh normalized metric rows from the curation ETL export, if present."""
    if not _table_exists(
        cursor,
        "data_libraries",
        "pgs_performance_metric_values",
    ):
        print(
            "Skipping normalized metric refresh: "
            "pgs_performance_metric_values does not exist yet."
        )
        return

    # Never destroy the existing normalized metric store merely because
    # the optional ETL export is unavailable. The report depends on the DB,
    # so preserving the last valid curated metric load is safer.
    if not csv_path.exists():
        print(
            f"Normalized metric file not found: {csv_path}. "
            "Preserving the existing pgs_performance_metric_values rows."
        )
        return

    metric_df = pd.read_csv(csv_path, low_memory=False)
    required = {
        "ppm_id", "source", "metric_type", "metric_label_raw",
        "estimate", "ci_lower", "ci_upper", "effect_basis",
        "personalization_allowed", "display_mode", "raw_metric",
    }
    missing = required - set(metric_df.columns)
    if missing:
        raise ValueError(
            f"Normalized metric CSV is missing columns: {sorted(missing)}"
        )

    cursor.execute("SELECT ppm_id FROM data_libraries.pgs_performance")
    valid_ppm_ids = {str(row[0]).strip() for row in cursor.fetchall()}

    cursor.execute("SELECT metric_type FROM data_libraries.performance_metric_types")
    valid_metric_types = {str(row[0]).strip() for row in cursor.fetchall()}

    metric_df["ppm_id"] = metric_df["ppm_id"].astype(str).str.strip()
    metric_df["metric_type"] = metric_df["metric_type"].astype(str).str.strip()

    # The curated long export uses "dedicated_column"; the DB vocabulary uses
    # "dedicated". Normalize before filtering so dedicated OR/HR/Beta/AUROC/
    # C-index rows are not silently discarded.
    metric_df["source"] = metric_df["source"].apply(
        _normalize_metric_source
    )

    metric_df = metric_df[
        metric_df["ppm_id"].isin(valid_ppm_ids)
        & metric_df["metric_type"].isin(valid_metric_types)
        & metric_df["source"].isin(["dedicated", "other_metric"])
    ].copy()

    # Only replace the live normalized table after the new dataset has been
    # found, parsed and filtered successfully.
    cursor.execute(
        "TRUNCATE TABLE "
        "data_libraries.pgs_performance_metric_values "
        "RESTART IDENTITY"
    )

    metric_df["metric_index"] = (
        metric_df.groupby(
            ["ppm_id", "metric_type", "source"],
            sort=False,
        ).cumcount()
        + 1
    )

    def _flag(value):
        if value is None:
            return False
        try:
            if pd.isna(value):
                return False
        except TypeError:
            pass
        if isinstance(value, bool):
            return value
        return str(value).strip().lower() in {"true", "1", "yes", "y"}

    rows = []

    for row in metric_df.itertuples(index=False):
        personalization_allowed = _flag(
            row.personalization_allowed
        )

        effect_basis = _normalize_effect_basis(
            metric_type=row.metric_type,
            source=row.source,
            raw_basis=row.effect_basis,
            display_mode=row.display_mode,
            personalization_allowed=personalization_allowed,
        )

        contrast_text = _metric_contrast_text(
            metric_type=row.metric_type,
            raw_basis=row.effect_basis,
            metric_label_raw=row.metric_label_raw,
            display_mode=row.display_mode,
        )

        if effect_basis not in _ALLOWED_EFFECT_BASIS:
            raise ValueError(
                f"Invalid effect_basis {effect_basis!r} for "
                f"{row.ppm_id} / {row.metric_type}"
            )

        display_mode = _to_db_scalar(row.display_mode)

        if (
            display_mode is not None
            and len(str(display_mode)) > 100
        ):
            raise ValueError(
                f"display_mode exceeds VARCHAR(100) for "
                f"{row.ppm_id}: {display_mode!r}"
            )

        rows.append((
            str(row.ppm_id).strip(),
            str(row.metric_type).strip(),
            int(row.metric_index),
            str(row.source).strip(),
            _to_db_scalar(row.metric_label_raw),
            _to_db_scalar(row.estimate),
            _to_db_scalar(row.ci_lower),
            _to_db_scalar(row.ci_upper),
            effect_basis,
            personalization_allowed,
            display_mode,
            (
                (beta_unit_by_ppm or {}).get(str(row.ppm_id).strip())
                if str(row.metric_type).strip() == "beta"
                else None
            ),
            contrast_text,
            _to_db_scalar(row.raw_metric),
            "2026-09-14-v2",
        ))

    # Helpful diagnostics before touching PostgreSQL.
    effect_basis_counts = (
        pd.Series(
            [item[8] for item in rows],
            dtype="object",
        )
        .value_counts(dropna=False)
        .to_dict()
    )

    source_counts = (
        pd.Series(
            [item[3] for item in rows],
            dtype="object",
        )
        .value_counts(dropna=False)
        .to_dict()
    )

    print(
        f"Normalized metric sources: {source_counts}"
    )
    print(
        f"Normalized effect bases: {effect_basis_counts}"
    )

    if rows:
        execute_values(
            cursor,
            """
            INSERT INTO data_libraries.pgs_performance_metric_values (
                ppm_id,
                metric_type,
                metric_index,
                source,
                metric_label_raw,
                estimate,
                ci_lower,
                ci_upper,
                effect_basis,
                personalization_allowed,
                display_mode,
                unit,
                contrast_text,
                raw_metric,
                parser_version
            ) VALUES %s
            """,
            rows,
            page_size=2000,
        )

    covered_ppm_ids = set(metric_df["ppm_id"])
    print(
        f"Inserted {len(rows)} normalized metric rows covering "
        f"{len(covered_ppm_ids)} current PPMs from {csv_path}."
    )

    uncovered = len(valid_ppm_ids - covered_ppm_ids)
    if uncovered:
        print(
            f"Warning: {uncovered} current PPMs have no normalized metric row. "
            "Regenerate the curation export when the PGS Catalog metadata changes."
        )


# -----------------------------
# 1) Read source sheets
# -----------------------------
pgscatalog_raw = pd.read_excel(excel_filepath, sheet_name="Scores")
publications_raw = pd.read_excel(excel_filepath, sheet_name="Publications")
ontology_raw = pd.read_excel(excel_filepath, sheet_name="EFO Traits")
performance_raw = pd.read_excel(excel_filepath, sheet_name="Performance Metrics")
score_dev_raw = pd.read_excel(excel_filepath, sheet_name="Score Development Samples")
eval_sets_raw = pd.read_excel(excel_filepath, sheet_name="Evaluation Sample Sets")

beta_unit_by_ppm = _build_beta_unit_by_ppm(
    performance_raw,
    eval_sets_raw,
)


# -----------------------------
# 2) Map to DB schema columns
# -----------------------------
pgscatalog_data = pd.DataFrame({
    "pgs_id": _first_existing(pgscatalog_raw, ["Polygenic Score (PGS) ID"]),
    "pgs_name": _first_existing(pgscatalog_raw, ["PGS Name"], required=False),
    "reported_trait": _first_existing(pgscatalog_raw, ["Reported Trait"], required=False),
    "mapped_trait_efo_label": _first_existing(pgscatalog_raw, ["Mapped Trait(s) (EFO label)"], required=False),
    "efo_id": _first_existing(pgscatalog_raw, ["Mapped Trait(s) (EFO ID)"], required=False),
    "pgs_development_method": _first_existing(pgscatalog_raw, ["PGS Development Method"], required=False),
    "pgs_development_details": _first_existing(pgscatalog_raw, ["PGS Development Details/Relevant Parameters"], required=False),
    "original_genome_build": _first_existing(pgscatalog_raw, ["Original Genome Build"], required=False),
    "number_of_variants": _to_nullable_int(_first_existing(pgscatalog_raw, ["Number of Variants"], required=False)),
    "number_of_interaction_terms": _to_nullable_int(_first_existing(pgscatalog_raw, ["Number of Interaction Terms"], required=False)),
    "type_of_variant_weight": _first_existing(pgscatalog_raw, ["Type of Variant Weight"], required=False),
    "pgp_id": _first_existing(pgscatalog_raw, ["PGS Publication (PGP) ID"], required=False),
    "publication_pmid": _to_nullable_int(_first_existing(pgscatalog_raw, ["Publication (PMID)"], required=False)),
    "publication_doi": _first_existing(pgscatalog_raw, ["Publication (doi)"], required=False),
    "score_and_results_match_original_publication": _to_bool(_first_existing(pgscatalog_raw, ["Score and results match the original publication"], required=False)),
    "ancestry_distribution_source_of_variant_associations_gwas": _first_existing(pgscatalog_raw, ["Ancestry Distribution (%) - Source of Variant Associations (GWAS)"], required=False),
    "ancestry_distribution_score_development_training": _first_existing(pgscatalog_raw, ["Ancestry Distribution (%) - Score Development/Training"], required=False),
    "ancestry_distribution_pgs_evaluation": _first_existing(pgscatalog_raw, ["Ancestry Distribution (%) - PGS Evaluation"], required=False),
    "ftp_link": _first_existing(pgscatalog_raw, ["FTP link"], required=False),
    "release_date": pd.to_datetime(_first_existing(pgscatalog_raw, ["Release Date"], required=False), errors="coerce").dt.date,
    "license_terms_of_use": _first_existing(pgscatalog_raw, ["License/Terms of Use"], required=False),
})

pgs_publications = pd.DataFrame({
    "pgp_id": _first_existing(publications_raw, ["PGS Publication/Study (PGP) ID"]),
    "first_author": _first_existing(publications_raw, ["First Author"], required=False),
    "title": _first_existing(publications_raw, ["Title"], required=False),
    "journal_name": _first_existing(publications_raw, ["Journal Name"], required=False),
    "publication_date": pd.to_datetime(_first_existing(publications_raw, ["Publication Date"], required=False), errors="coerce").dt.date,
    "release_date": pd.to_datetime(_first_existing(publications_raw, ["Release Date"], required=False), errors="coerce").dt.date,
    "authors": _first_existing(publications_raw, ["Authors"], required=False),
    "digital_object_identifier_doi": _first_existing(publications_raw, ["digital object identifier (doi)"], required=False),
    "pubmed_id_pmid": _first_existing(publications_raw, ["PubMed ID (PMID)"], required=False),
})

ontology_mappings = pd.DataFrame({
    "ontology_id": _first_existing(ontology_raw, ["Ontology Trait ID"]),
    "ontology_label": _first_existing(ontology_raw, ["Ontology Trait Label"], required=False),
    "ontology_description": _first_existing(ontology_raw, ["Ontology Trait Description"], required=False),
    "ontology_url": _first_existing(ontology_raw, ["Ontology URL"], required=False),
})

hazard_point, hazard_ci_lower, hazard_ci_upper = _parse_effect_with_ci(
    _first_existing(performance_raw, ["Hazard Ratio (HR)"], required=False)
)
odds_point, odds_ci_lower, odds_ci_upper = _parse_effect_with_ci(
    _first_existing(performance_raw, ["Odds Ratio (OR)"], required=False)
)
beta_point, beta_ci_lower, beta_ci_upper = _parse_effect_with_ci(
    _first_existing(performance_raw, ["Beta"], required=False)
)
auroc_point, auroc_ci_lower, auroc_ci_upper = _parse_effect_with_ci(
    _first_existing(
        performance_raw,
        ["Area Under the Receiver-Operating Characteristic Curve (AUROC)"],
        required=False,
    )
)
cindex_point, cindex_ci_lower, cindex_ci_upper = _parse_effect_with_ci(
    _first_existing(
        performance_raw,
        ["Concordance Statistic (C-index)"],
        required=False,
    )
)

pgs_performance = pd.DataFrame({
    "ppm_id": _first_existing(performance_raw, ["PGS Performance Metric (PPM) ID"]),
    "pgs_id": _first_existing(performance_raw, ["Evaluated Score"]),
    "pss_id": _first_existing(performance_raw, ["PGS Sample Set (PSS)"], required=False),
    "pgp_id": _first_existing(performance_raw, ["PGS Publication (PGP) ID"], required=False),
    "reported_trait": _first_existing(performance_raw, ["Reported Trait"], required=False),
    "covariates_included_in_model": _first_existing(performance_raw, ["Covariates Included in the Model"], required=False),
    "pgs_performance_other_relevant_info": _first_existing(performance_raw, ["PGS Performance: Other Relevant Information"], required=False),
    "publication_pmid": _to_nullable_int(_first_existing(performance_raw, ["Publication (PMID)"], required=False)),
    "publication_doi": _first_existing(performance_raw, ["Publication (doi)"], required=False),
    "hazard_ratio": hazard_point,
    "hazard_ratio_ci_lower": hazard_ci_lower,
    "hazard_ratio_ci_upper": hazard_ci_upper,
    "odds_ratio": odds_point,
    "odds_ratio_ci_lower": odds_ci_lower,
    "odds_ratio_ci_upper": odds_ci_upper,
    "beta": beta_point,
    "beta_ci_lower": beta_ci_lower,
    "beta_ci_upper": beta_ci_upper,
    "auroc": auroc_point,
    "auroc_ci_lower": auroc_ci_lower,
    "auroc_ci_upper": auroc_ci_upper,
    "concordance_statistic": cindex_point,
    "concordance_statistic_ci_lower": cindex_ci_lower,
    "concordance_statistic_ci_upper": cindex_ci_upper,
    "other_metric": _first_existing(performance_raw, ["Other Metric(s)"], required=False),
})

score_development_samples = pd.DataFrame({
    "pgs_id": _first_existing(score_dev_raw, ["Polygenic Score (PGS) ID"]),
    "stage_of_pgs_development": _first_existing(score_dev_raw, ["Stage of PGS Development"], required=False),
    "individuals_development": _to_nullable_int(_first_existing(score_dev_raw, ["Number of Individuals"], required=False)),
    "cases_development": _to_nullable_int(_first_existing(score_dev_raw, ["Number of Cases"], required=False)),
    "controls_development": _to_nullable_int(_first_existing(score_dev_raw, ["Number of Controls"], required=False)),
    "percent_male_development": pd.to_numeric(_first_existing(score_dev_raw, ["Percent of Participants Who are Male"], required=False), errors="coerce"),
    "broad_ancestry_category": _first_existing(score_dev_raw, ["Broad Ancestry Category", "Broad Ancestry Category(ies)", "Broad ancestry category"], required=False),
})

evaluation_sample_sets = pd.DataFrame({
    "pss_id": _first_existing(eval_sets_raw, ["PGS Sample Set (PSS)"]),
    "pgs_id": _first_existing(eval_sets_raw, ["Polygenic Score (PGS) ID", "Evaluated Score", "PGS ID"]),
    "individuals_evaluation": _to_nullable_int(_first_existing(eval_sets_raw, ["Number of Individuals"], required=False)),
    "cases_evaluation": _to_nullable_int(_first_existing(eval_sets_raw, ["Number of Cases"], required=False)),
    "controls_evaluation": _to_nullable_int(_first_existing(eval_sets_raw, ["Number of Controls"], required=False)),
    "percent_male_evaluation": pd.to_numeric(_first_existing(eval_sets_raw, ["Percent of Participants Who are Male"], required=False), errors="coerce"),
    "broad_ancestry_category": _first_existing(eval_sets_raw, ["Broad Ancestry Category", "Broad Ancestry Category(ies)", "Broad ancestry category"], required=False),
    "country_of_recruitment": _first_existing(eval_sets_raw, ["Country of Recruitment", "Country/Countries of Recruitment", "Countries of recruitment"], required=False),
    "sample_age_raw": _first_existing(eval_sets_raw, ["Sample Age"], required=False),
    "ancestry_detail": _first_existing(eval_sets_raw, ["Ancestry (e.g. French, Chinese)"], required=False),
    "additional_ancestry_description": _first_existing(eval_sets_raw, ["Additional Ancestry Description"], required=False),
    "phenotype_definitions_methods": _first_existing(eval_sets_raw, ["Phenotype Definitions and Methods"], required=False),
    "followup_time": _first_existing(eval_sets_raw, ["Followup Time"], required=False),
    "additional_sample_cohort_information": _first_existing(eval_sets_raw, ["Additional Sample/Cohort Information"], required=False),
    "cohort": _first_existing(
        eval_sets_raw,
        [
            "Cohort",
            "Cohort(s)",
            "Cohorts",
            "Cohort(s) Included",
            "Cohort Name",
        ],
        required=False,
    ),
})

# Expand many-to-many style PSS->PGS links into one row per pgs_id.
evaluation_sample_sets["pgs_id"] = evaluation_sample_sets["pgs_id"].apply(_split_pgs_ids)
evaluation_sample_sets = evaluation_sample_sets.explode("pgs_id", ignore_index=True)

# Remove blank keys so inserts do not fail on PK-only empty strings.
for dataframe, pk_col in [
    (pgscatalog_data, "pgs_id"),
    (pgs_publications, "pgp_id"),
    (ontology_mappings, "ontology_id"),
    (pgs_performance, "ppm_id"),
    (score_development_samples, "pgs_id"),
    (evaluation_sample_sets, "pss_id"),
]:
    dataframe[pk_col] = dataframe[pk_col].astype(str).str.strip()
    dataframe.drop(dataframe[dataframe[pk_col].isin(["", "nan", "None"])].index, inplace=True)

# Clean exploded pgs_id values and drop rows that still have no usable pgs_id.
evaluation_sample_sets["pgs_id"] = evaluation_sample_sets["pgs_id"].astype(str).str.strip()
evaluation_sample_sets.drop(
    evaluation_sample_sets[evaluation_sample_sets["pgs_id"].isin(["", "nan", "None"])].index,
    inplace=True,
)

evaluation_sample_sets = evaluation_sample_sets.drop_duplicates(
    subset=[
        "pss_id", "pgs_id", "individuals_evaluation", "cases_evaluation",
        "controls_evaluation", "percent_male_evaluation",
        "broad_ancestry_category", "country_of_recruitment", "cohort",
        "sample_age_raw", "ancestry_detail", "additional_ancestry_description",
        "phenotype_definitions_methods", "followup_time",
        "additional_sample_cohort_information",
    ]
)

pgscatalog_data = _clean_df_for_sql(pgscatalog_data)
pgs_publications = _clean_df_for_sql(pgs_publications)
ontology_mappings = _clean_df_for_sql(ontology_mappings)
pgs_performance = _clean_df_for_sql(pgs_performance)
score_development_samples = _clean_df_for_sql(score_development_samples)
evaluation_sample_sets = _clean_df_for_sql(evaluation_sample_sets)


# -----------------------------
# 3) Connect and load into DB
# -----------------------------
DB_USER = os.getenv("DB_USER", os.getenv("POSTGRES_USER", "postgres"))
DB_PASSWORD = os.getenv("DB_PASSWORD", os.getenv("POSTGRES_PASSWORD", "zod50902"))
DB_NAME = os.getenv("DB_NAME", os.getenv("POSTGRES_DB", "snpster_db"))
DB_HOST = os.getenv("DB_HOST", "127.0.0.1")
DB_PORT = int(os.getenv("DB_PORT", "5433"))

conn = psycopg2.connect(
    dbname=DB_NAME,
    user=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
)

try:
    with conn.cursor() as cur:
        # Compatibility for older live schemas still using legacy column names.
        eval_cols = set(_get_table_columns(cur, "data_libraries", "evaluation_sample_sets"))
        if "pgs_id" not in eval_cols and "evaluated_score" in eval_cols and "pgs_id" in evaluation_sample_sets.columns:
            evaluation_sample_sets = evaluation_sample_sets.rename(columns={"pgs_id": "evaluated_score"})
            print("Mapped evaluation_sample_sets.pgs_id -> evaluated_score for live DB compatibility")

        perf_cols = set(_get_table_columns(cur, "data_libraries", "pgs_performance"))
        if "ppm_id" not in perf_cols and "performance_id" in perf_cols and "ppm_id" in pgs_performance.columns:
            pgs_performance = pgs_performance.rename(columns={"ppm_id": "performance_id"})
            print("Mapped pgs_performance.ppm_id -> performance_id for live DB compatibility")

        # Preserve parent/catalog rows so user PRS results and curated report
        # mappings keep their foreign-key targets across metadata refreshes.
        child_tables = [
            "data_libraries.pgs_performance",
            "data_libraries.evaluation_sample_sets",
            "data_libraries.score_development_samples",
        ]
        if _table_exists(
            cur,
            "data_libraries",
            "pgs_performance_metric_values",
        ):
            child_tables.insert(
                0,
                "data_libraries.pgs_performance_metric_values",
            )

        cur.execute(
            "TRUNCATE TABLE "
            + ", ".join(child_tables)
            + " RESTART IDENTITY"
        )

        _bulk_upsert(
            cur,
            "data_libraries",
            "pgscatalog_data",
            pgscatalog_data,
            ["pgs_id"],
        )
        _bulk_upsert(
            cur,
            "data_libraries",
            "pgs_publications",
            pgs_publications,
            ["pgp_id"],
        )
        _bulk_upsert(
            cur,
            "data_libraries",
            "ontology_mappings",
            ontology_mappings,
            ["ontology_id"],
        )

        _bulk_insert(
            cur,
            "data_libraries",
            "score_development_samples",
            score_development_samples,
        )
        _bulk_insert(
            cur,
            "data_libraries",
            "evaluation_sample_sets",
            evaluation_sample_sets,
        )
        _bulk_insert(
            cur,
            "data_libraries",
            "pgs_performance",
            pgs_performance,
        )

        # Optional: refresh normalized metric rows from the curation ETL export.
        # This is an update-time dependency only; reports read PostgreSQL directly.
        _load_normalized_metric_values(
            cur,
            metric_long_filepath,
            beta_unit_by_ppm=beta_unit_by_ppm,
        )

    conn.commit()
    print("Metadata refresh complete; curated taxonomy and user PRS rows were preserved.")
except Exception:
    conn.rollback()
    raise
finally:
    conn.close()