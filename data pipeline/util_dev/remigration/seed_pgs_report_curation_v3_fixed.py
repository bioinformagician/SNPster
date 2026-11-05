import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DEFAULT_CURATION_DIR = Path(
    os.getenv(
        "PGS_CURATION_DIR",
        "/home/frederik/github_projects/SNPster/data pipeline/util_dev/remigration",
    )
)

CURATED_SCORE_FILE = Path(
    os.getenv(
        "PGS_CURATED_SCORE_FILE",
        DEFAULT_CURATION_DIR / "pgs_reports_shop_curated.csv",
    )
)
OUTCOME_MAPPING_FILE = Path(
    os.getenv(
        "PGS_OUTCOME_MAPPING_FILE",
        DEFAULT_CURATION_DIR / "evaluation_outcome_mapping.csv",
    )
)
METRIC_LONG_FILE = Path(
    os.getenv(
        "PGS_METRIC_LONG_FILE",
        DEFAULT_CURATION_DIR / "pgs_performance_metrics_long_context_v3.csv",
    )
)
PPM_CONTEXT_FILE = Path(
    os.getenv(
        "PGS_PPM_CONTEXT_FILE",
        DEFAULT_CURATION_DIR / "pgs_performance_context_curation_v3.csv",
    )
)
METRIC_OVERVIEW_FILE = Path(
    os.getenv(
        "PGS_METRIC_OVERVIEW_FILE",
        DEFAULT_CURATION_DIR / "pgs_metric_overview.csv",
    )
)

CURATION_VERSION = os.getenv("PGS_CURATION_VERSION", "2026-09-15-v3")
PARSER_VERSION = os.getenv("PGS_METRIC_PARSER_VERSION", "2026-09-15-v3")

DB_USER = os.getenv("DB_USER", os.getenv("POSTGRES_USER", "postgres"))
DB_PASSWORD = os.getenv("DB_PASSWORD", os.getenv("POSTGRES_PASSWORD", "zod50902"))
DB_NAME = os.getenv("DB_NAME", os.getenv("POSTGRES_DB", "snpster_db"))
DB_HOST = os.getenv("DB_HOST", "127.0.0.1")
DB_PORT = int(os.getenv("DB_PORT", "5433"))


CATEGORY_CONFIG = [
    (
        "cardiovascular_cholesterol",
        "Cardiovascular & Cholesterol",
        10,
        "Heart and vascular traits including coronary disease, blood pressure, "
        "lipids, rhythm-related conditions, and other cardiovascular phenotypes.",
    ),
    (
        "cancer",
        "Cancer",
        20,
        "Polygenic predisposition across different cancer types. These scores "
        "describe inherited susceptibility and do not replace clinical screening "
        "or high-penetrance cancer-gene testing.",
    ),
    (
        "metabolism_hormones",
        "Metabolism & Hormones",
        30,
        "Metabolic and endocrine traits including diabetes, body composition, "
        "thyroid-related phenotypes, and hormone-associated conditions.",
    ),
    (
        "immune_autoimmune",
        "Immune & Autoimmune",
        40,
        "Immune-mediated, inflammatory, allergic, and autoimmune conditions in "
        "which inherited variation can influence immune-system activity.",
    ),
    (
        "brain_mental_health",
        "Brain & Mental Health",
        50,
        "Neurological, cognitive, neurodegenerative, and psychiatric traits that "
        "relate to brain function, behavior, and mental health.",
    ),
    (
        "blood_biomarkers",
        "Blood & Biomarkers",
        60,
        "Blood-cell traits and circulating biomarkers, including measures related "
        "to lipid metabolism, inflammation, nutrient status, and physiology.",
    ),
    (
        "lung_respiratory",
        "Lung & Respiratory",
        70,
        "Respiratory traits including asthma, chronic lung disease, pulmonary "
        "function, and other conditions affecting breathing.",
    ),
    (
        "eye_health",
        "Eye Health",
        80,
        "Ophthalmic traits involving vision and clinically meaningful eye disease.",
    ),
    (
        "other_health_conditions",
        "Other Health Conditions",
        90,
        "Consumer-relevant health traits that do not fit cleanly into the other "
        "major report categories.",
    ),
]

CATEGORY_BY_NAME = {
    display_name: {
        "category_key": category_key,
        "display_order": display_order,
        "description": description,
    }
    for category_key, display_name, display_order, description in CATEGORY_CONFIG
}


METRIC_FAMILY = {
    "odds_ratio": "effect",
    "hazard_ratio": "effect",
    "group_odds_ratio": "effect",
    "group_hazard_ratio": "effect",
    "relative_risk": "effect",
    "relative_risk_ratio": "effect",
    "incidence_rate_ratio": "effect",
    "beta": "effect",
    "beta_other": "effect",
    "mean_difference": "effect",
    "risk_difference": "effect",
    "r2": "variance_explained",
    "incremental_r2": "variance_explained",
    "partial_r2": "variance_explained",
    "adjusted_r2": "variance_explained",
    "liability_r2": "variance_explained",
    "pseudo_r2": "model_fit",
    "auroc": "discrimination",
    "auroc_other": "discrimination",
    "c_index": "discrimination",
    "cindex_other": "discrimination",
    "auprc": "discrimination",
    "delta_auroc": "discrimination_change",
    "delta_cindex": "discrimination_change",
    "correlation": "association",
    "partial_correlation": "association",
    "spearman_rho": "association",
    "brier_score": "calibration",
    "calibration": "calibration",
    "nri": "reclassification",
    "idi": "reclassification",
    "ppv": "classification",
    "npv": "classification",
    "sensitivity": "classification",
    "specificity": "classification",
    "youden_index": "classification",
    "f_measure": "classification",
    "p_value": "statistical_significance",
    "incidence_percent": "descriptive",
    "other": "other",
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required curation file not found: {path}")


def slugify(value: str) -> str:
    value = str(value).strip().casefold()
    value = value.replace("&", " and ")
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def db_scalar(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def as_bool(value, default=False):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().casefold()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return default


def fetch_existing_ids(cur, table, id_column):
    cur.execute(f"SELECT {id_column} FROM data_libraries.{table}")
    return {str(row[0]).strip() for row in cur.fetchall()}


ALLOWED_EFFECT_BASIS = {
    "per_sd",
    "per_sd_assumed",
    "group_contrast",
    "absolute",
    "unknown",
    "not_applicable",
}


def normalize_metric_source(value):
    """Map parser source labels to the DB controlled vocabulary."""
    if value is None:
        return None

    source = str(value).strip().casefold()

    return {
        "dedicated_column": "dedicated",
        "dedicated": "dedicated",
        "other_metric": "other_metric",
    }.get(source)


def normalize_effect_basis(
    metric_type,
    source,
    raw_basis,
    display_mode,
    personalization_allowed,
):
    """Return only values accepted by pgs_performance_metric_values.effect_basis."""
    metric_type = "" if metric_type is None else str(metric_type).strip()
    source = "" if source is None else str(source).strip()
    display_mode = "" if display_mode is None else str(display_mode).strip()

    if raw_basis is None:
        raw_basis_text = ""
    else:
        try:
            raw_basis_text = "" if pd.isna(raw_basis) else str(raw_basis).strip()
        except TypeError:
            raw_basis_text = str(raw_basis).strip()

    if raw_basis_text in ALLOWED_EFFECT_BASIS:
        return raw_basis_text

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


def group_contrast_text(
    metric_type,
    source,
    raw_basis,
    metric_label_raw,
    display_mode,
):
    """Preserve comparison wording separately from effect_basis."""
    is_group_comparison = (
        str(display_mode).strip() == "group_comparison"
        or str(metric_type).strip() in {
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


# -----------------------------------------------------------------------------
# Load files
# -----------------------------------------------------------------------------

for required_path in [
    CURATED_SCORE_FILE,
    OUTCOME_MAPPING_FILE,
    METRIC_LONG_FILE,
    PPM_CONTEXT_FILE,
    METRIC_OVERVIEW_FILE,
]:
    require_file(required_path)

curated_scores = pd.read_csv(CURATED_SCORE_FILE, dtype=str).fillna("")
outcome_mappings = pd.read_csv(OUTCOME_MAPPING_FILE, dtype=str).fillna("")
metric_long = pd.read_csv(METRIC_LONG_FILE, low_memory=False)
ppm_context = pd.read_csv(PPM_CONTEXT_FILE, low_memory=False)
metric_overview = pd.read_csv(METRIC_OVERVIEW_FILE, dtype=str).fillna("")

required_score_columns = {"pgs_id", "reported_trait", "category"}
missing = required_score_columns - set(curated_scores.columns)
if missing:
    raise ValueError(f"Curated score file missing columns: {sorted(missing)}")

required_outcome_columns = {
    "evaluation_outcome_raw",
    "evaluation_outcome_clean",
    "evaluation_outcome_key",
    "evaluation_target_trait",
    "evaluation_target_category",
}
missing = required_outcome_columns - set(outcome_mappings.columns)
if missing:
    raise ValueError(f"Outcome mapping file missing columns: {sorted(missing)}")

required_metric_columns = {
    "ppm_id",
    "source",
    "metric_type",
    "metric_label_raw",
    "estimate",
    "ci_lower",
    "ci_upper",
    "effect_basis",
    "personalization_allowed",
    "display_mode",
    "raw_metric",
    "metric_family",
    "metric_model_scope",
    "metric_model_scope_confidence",
    "metric_model_scope_reason",
    "r2_scale_basis",
    "r2_estimate_canonical_fraction",
    "r2_estimate_display_percent",
    "r2_ci_lower_canonical_fraction",
    "r2_ci_upper_canonical_fraction",
    "r2_scale_reason",
    "quality_plot_role",
    "quality_context_warning",
    "preferred_for_primary_quality",
}
missing = required_metric_columns - set(metric_long.columns)
if missing:
    raise ValueError(f"Long metric file missing columns: {sorted(missing)}")

required_ppm_context_columns = {
    "ppm_id",
    "has_covariates",
    "covariate_class",
    "clinical_predictors_present",
    "clinical_predictor_terms",
    "population_scope_hint",
    "population_scope_confidence",
    "population_scope_reason",
    "applicability_status",
    "applicability_requires_match",
    "applicability_restriction_types",
    "applicability_reason",
    "age_restriction_min_years",
    "age_restriction_max_years",
    "sex_restriction",
    "condition_restriction",
    "sample_age_raw",
    "sample_age_mean_years",
    "sample_age_min_years",
    "sample_age_max_years",
    "phenotype_definition_raw",
    "additional_sample_info_raw",
}
missing = required_ppm_context_columns - set(ppm_context.columns)
if missing:
    raise ValueError(f"PPM context file missing columns: {sorted(missing)}")

curated_scores["pgs_id"] = curated_scores["pgs_id"].str.strip()
curated_scores["reported_trait"] = curated_scores["reported_trait"].str.strip()
curated_scores["category"] = curated_scores["category"].str.strip()
curated_scores = curated_scores[
    curated_scores["pgs_id"].str.fullmatch(r"PGS\d{6}", na=False)
    & curated_scores["reported_trait"].ne("")
    & curated_scores["category"].ne("")
].drop_duplicates(subset=["pgs_id"], keep="first")

unexpected_categories = sorted(
    set(curated_scores["category"]) - set(CATEGORY_BY_NAME)
)
if unexpected_categories:
    raise ValueError(
        "Curated score file contains categories not configured in the seed script: "
        + ", ".join(unexpected_categories)
    )


# -----------------------------------------------------------------------------
# Seed database
# -----------------------------------------------------------------------------

conn = psycopg2.connect(
    dbname=DB_NAME,
    user=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
)

try:
    with conn.cursor() as cur:
        # -------------------------------------------------------------
        # Categories
        # -------------------------------------------------------------
        category_rows = [
            (
                category_key,
                display_name,
                description,
                display_order,
                True,
            )
            for category_key, display_name, display_order, description
            in CATEGORY_CONFIG
        ]

        execute_values(
            cur,
            """
            INSERT INTO data_libraries.report_categories (
                category_key,
                display_name,
                description,
                display_order,
                is_active
            ) VALUES %s
            ON CONFLICT (category_key)
            DO UPDATE SET
                display_name = EXCLUDED.display_name,
                description = EXCLUDED.description,
                display_order = EXCLUDED.display_order,
                is_active = EXCLUDED.is_active,
                updated_at = CURRENT_TIMESTAMP
            """,
            category_rows,
            page_size=100,
        )

        cur.execute(
            """
            SELECT category_id, category_key, display_name
            FROM data_libraries.report_categories
            """
        )
        category_lookup = {
            row[2]: {"category_id": row[0], "category_key": row[1]}
            for row in cur.fetchall()
        }

        # -------------------------------------------------------------
        # Consumer traits
        # -------------------------------------------------------------
        trait_pairs = (
            curated_scores[["reported_trait", "category"]]
            .drop_duplicates()
            .sort_values(["category", "reported_trait"])
        )

        trait_rows = []
        for row in trait_pairs.itertuples(index=False):
            category = category_lookup[row.category]
            trait_key = (
                f"{category['category_key']}__{slugify(row.reported_trait)}"
            )
            trait_rows.append(
                (
                    trait_key,
                    row.reported_trait,
                    category["category_id"],
                    True,
                )
            )

        execute_values(
            cur,
            """
            INSERT INTO data_libraries.report_traits (
                trait_key,
                display_name,
                category_id,
                is_active
            ) VALUES %s
            ON CONFLICT (trait_key)
            DO UPDATE SET
                display_name = EXCLUDED.display_name,
                category_id = EXCLUDED.category_id,
                is_active = EXCLUDED.is_active,
                updated_at = CURRENT_TIMESTAMP
            """,
            trait_rows,
            page_size=1000,
        )

        cur.execute(
            """
            SELECT
                rt.report_trait_id,
                rt.display_name,
                rc.display_name
            FROM data_libraries.report_traits AS rt
            JOIN data_libraries.report_categories AS rc
              ON rc.category_id = rt.category_id
            """
        )
        trait_lookup = {
            (row[1], row[2]): row[0]
            for row in cur.fetchall()
        }

        # Mark traits not present in the current curation inactive, rather than
        # deleting them and breaking historical references.
        active_trait_keys = [row[0] for row in trait_rows]
        cur.execute(
            """
            UPDATE data_libraries.report_traits
            SET is_active = FALSE,
                updated_at = CURRENT_TIMESTAMP
            WHERE NOT (trait_key = ANY(%s))
            """,
            (active_trait_keys,),
        )

        # -------------------------------------------------------------
        # PGS -> consumer trait mapping
        # -------------------------------------------------------------
        existing_pgs_ids = fetch_existing_ids(cur, "pgscatalog_data", "pgs_id")
        missing_pgs_ids = sorted(
            set(curated_scores["pgs_id"]) - existing_pgs_ids
        )
        if missing_pgs_ids:
            print(
                f"Warning: {len(missing_pgs_ids)} curated PGS IDs are not present "
                "in data_libraries.pgscatalog_data and will be skipped."
            )
            print("First missing IDs:", missing_pgs_ids[:20])

        score_rows = []
        for row in curated_scores.itertuples(index=False):
            if row.pgs_id not in existing_pgs_ids:
                continue
            report_trait_id = trait_lookup[(row.reported_trait, row.category)]
            score_rows.append(
                (
                    row.pgs_id,
                    report_trait_id,
                    True,
                    CURATION_VERSION,
                )
            )

        # The curated CSV is authoritative for consumer-card membership.
        cur.execute("DELETE FROM data_libraries.report_trait_scores")
        execute_values(
            cur,
            """
            INSERT INTO data_libraries.report_trait_scores (
                pgs_id,
                report_trait_id,
                include_in_report,
                curation_version
            ) VALUES %s
            """,
            score_rows,
            page_size=1000,
        )

        # -------------------------------------------------------------
        # Evaluation outcome mapping
        # -------------------------------------------------------------
        outcome_rows = []
        unresolved_targets = []

        for row in outcome_mappings.itertuples(index=False):
            raw_outcome = str(row.evaluation_outcome_raw).strip()
            clean_outcome = str(row.evaluation_outcome_clean).strip()
            outcome_key = str(row.evaluation_outcome_key).strip()
            target_trait = str(row.evaluation_target_trait).strip()
            target_category = str(row.evaluation_target_category).strip()

            if not raw_outcome:
                continue

            target_report_trait_id = None
            mapping_status = "unmapped"

            if target_trait and target_category:
                target_report_trait_id = trait_lookup.get(
                    (target_trait, target_category)
                )
                if target_report_trait_id is not None:
                    mapping_status = "mapped"
                else:
                    mapping_status = "needs_review"
                    unresolved_targets.append(
                        (raw_outcome, target_trait, target_category)
                    )

            outcome_rows.append(
                (
                    raw_outcome,
                    clean_outcome or raw_outcome,
                    outcome_key or slugify(clean_outcome or raw_outcome),
                    target_report_trait_id,
                    mapping_status,
                    CURATION_VERSION,
                )
            )

        cur.execute("DELETE FROM data_libraries.evaluation_outcome_mappings")
        execute_values(
            cur,
            """
            INSERT INTO data_libraries.evaluation_outcome_mappings (
                raw_outcome,
                cleaned_outcome,
                outcome_key,
                target_report_trait_id,
                mapping_status,
                curation_version
            ) VALUES %s
            """,
            outcome_rows,
            page_size=1000,
        )

        # -------------------------------------------------------------
        # PPM-level applicability and model-context curation
        # -------------------------------------------------------------
        existing_ppm_ids = fetch_existing_ids(cur, "pgs_performance", "ppm_id")
        ppm_context_current = ppm_context.copy()
        ppm_context_current["ppm_id"] = (
            ppm_context_current["ppm_id"].astype(str).str.strip()
        )
        ppm_context_current = ppm_context_current[
            ppm_context_current["ppm_id"].isin(existing_ppm_ids)
        ].copy()

        context_rows = []
        for row in ppm_context_current.itertuples(index=False):
            context_rows.append((
                str(row.ppm_id).strip(),
                as_bool(row.has_covariates),
                db_scalar(row.covariate_class),
                as_bool(row.clinical_predictors_present),
                db_scalar(row.clinical_predictor_terms),
                db_scalar(row.population_scope_hint),
                db_scalar(row.population_scope_confidence),
                db_scalar(row.population_scope_reason),
                db_scalar(row.applicability_status),
                as_bool(row.applicability_requires_match),
                db_scalar(row.applicability_restriction_types),
                db_scalar(row.applicability_reason),
                db_scalar(row.age_restriction_min_years),
                db_scalar(row.age_restriction_max_years),
                db_scalar(row.sex_restriction),
                db_scalar(row.condition_restriction),
                db_scalar(row.sample_age_raw),
                db_scalar(row.sample_age_mean_years),
                db_scalar(row.sample_age_min_years),
                db_scalar(row.sample_age_max_years),
                db_scalar(row.phenotype_definition_raw),
                db_scalar(row.additional_sample_info_raw),
            ))

        if context_rows:
            execute_values(
                cur,
                """
                UPDATE data_libraries.pgs_performance AS pp
                SET
                    has_covariates = v.has_covariates,
                    covariate_class = v.covariate_class,
                    clinical_predictors_present = v.clinical_predictors_present,
                    clinical_predictor_terms = v.clinical_predictor_terms,
                    population_scope_hint = v.population_scope_hint,
                    population_scope_confidence = v.population_scope_confidence,
                    population_scope_reason = v.population_scope_reason,
                    applicability_status = v.applicability_status,
                    applicability_requires_match = v.applicability_requires_match,
                    applicability_restriction_types = v.applicability_restriction_types,
                    applicability_reason = v.applicability_reason,
                    age_restriction_min_years = NULLIF(v.age_restriction_min_years::text, '')::numeric,
                    age_restriction_max_years = NULLIF(v.age_restriction_max_years::text, '')::numeric,
                    sex_restriction = NULLIF(v.sex_restriction, ''),
                    condition_restriction = NULLIF(v.condition_restriction, ''),
                    sample_age_raw = NULLIF(v.sample_age_raw, ''),
                    sample_age_mean_years = NULLIF(v.sample_age_mean_years::text, '')::numeric,
                    sample_age_min_years = NULLIF(v.sample_age_min_years::text, '')::numeric,
                    sample_age_max_years = NULLIF(v.sample_age_max_years::text, '')::numeric,
                    phenotype_definition_raw = NULLIF(v.phenotype_definition_raw, ''),
                    additional_sample_info_raw = NULLIF(v.additional_sample_info_raw, '')
                FROM (VALUES %s) AS v(
                    ppm_id,
                    has_covariates,
                    covariate_class,
                    clinical_predictors_present,
                    clinical_predictor_terms,
                    population_scope_hint,
                    population_scope_confidence,
                    population_scope_reason,
                    applicability_status,
                    applicability_requires_match,
                    applicability_restriction_types,
                    applicability_reason,
                    age_restriction_min_years,
                    age_restriction_max_years,
                    sex_restriction,
                    condition_restriction,
                    sample_age_raw,
                    sample_age_mean_years,
                    sample_age_min_years,
                    sample_age_max_years,
                    phenotype_definition_raw,
                    additional_sample_info_raw
                )
                WHERE pp.ppm_id = v.ppm_id
                """,
                context_rows,
                page_size=1000,
            )

        print(
            f"Applied curated applicability context to {len(context_rows)} PPM rows."
        )

        # -------------------------------------------------------------
        # Metric type dictionary
        # -------------------------------------------------------------
        metric_type_rows = []
        for idx, row in enumerate(metric_overview.itertuples(index=False), start=1):
            metric_type = str(row.metric_type).strip()
            if not metric_type:
                continue
            metric_type_rows.append(
                (
                    metric_type,
                    str(row.display_name).strip() or metric_type,
                    METRIC_FAMILY.get(metric_type, "other"),
                    None,
                    str(row.report_role).strip() or None,
                    str(row.recommendation).strip() or None,
                    idx,
                    True,
                )
            )

        execute_values(
            cur,
            """
            INSERT INTO data_libraries.performance_metric_types (
                metric_type,
                display_name,
                metric_family,
                default_display_mode,
                report_role,
                recommendation,
                display_order,
                is_active
            ) VALUES %s
            ON CONFLICT (metric_type)
            DO UPDATE SET
                display_name = EXCLUDED.display_name,
                metric_family = EXCLUDED.metric_family,
                report_role = EXCLUDED.report_role,
                recommendation = EXCLUDED.recommendation,
                display_order = EXCLUDED.display_order,
                is_active = EXCLUDED.is_active,
                updated_at = CURRENT_TIMESTAMP
            """,
            metric_type_rows,
            page_size=200,
        )

        # -------------------------------------------------------------
        # Normalized metric values
        # -------------------------------------------------------------
        existing_ppm_ids = fetch_existing_ids(cur, "pgs_performance", "ppm_id")
        metric_long["ppm_id"] = metric_long["ppm_id"].astype(str).str.strip()
        metric_long["metric_type"] = metric_long["metric_type"].astype(str).str.strip()

        metric_long["source"] = metric_long["source"].apply(
            normalize_metric_source
        )

        unknown_source_rows = metric_long[
            metric_long["source"].isna()
        ]

        if not unknown_source_rows.empty:
            print(
                f"Skipping {len(unknown_source_rows)} metric rows with "
                "unrecognized source labels."
            )

        metric_long_current = metric_long[
            metric_long["ppm_id"].isin(existing_ppm_ids)
            & metric_long["metric_type"].ne("")
            & metric_long["source"].isin(["dedicated", "other_metric"])
        ].copy()

        metric_long_current["metric_index"] = (
            metric_long_current.groupby(
                ["ppm_id", "metric_type", "source"],
                sort=False,
            ).cumcount()
            + 1
        )

        metric_rows = []

        for row in metric_long_current.itertuples(index=False):
            personalization_allowed = as_bool(
                row.personalization_allowed
            )

            effect_basis = normalize_effect_basis(
                metric_type=row.metric_type,
                source=row.source,
                raw_basis=row.effect_basis,
                display_mode=row.display_mode,
                personalization_allowed=personalization_allowed,
            )

            contrast_text = group_contrast_text(
                metric_type=row.metric_type,
                source=row.source,
                raw_basis=row.effect_basis,
                metric_label_raw=row.metric_label_raw,
                display_mode=row.display_mode,
            )

            if effect_basis not in ALLOWED_EFFECT_BASIS:
                raise ValueError(
                    f"Invalid effect_basis {effect_basis!r} for "
                    f"{row.ppm_id} / {row.metric_type}"
                )

            metric_rows.append(
                (
                    row.ppm_id,
                    row.metric_type,
                    int(row.metric_index),
                    row.source,
                    db_scalar(row.metric_label_raw),
                    db_scalar(row.estimate),
                    db_scalar(row.ci_lower),
                    db_scalar(row.ci_upper),
                    effect_basis,
                    personalization_allowed,
                    db_scalar(row.display_mode),
                    db_scalar(getattr(row, "unit", None)),
                    contrast_text,
                    db_scalar(row.raw_metric),
                    db_scalar(row.metric_family),
                    db_scalar(row.metric_model_scope),
                    db_scalar(row.metric_model_scope_confidence),
                    db_scalar(row.metric_model_scope_reason),
                    db_scalar(row.r2_scale_basis),
                    db_scalar(row.r2_estimate_canonical_fraction),
                    db_scalar(row.r2_estimate_display_percent),
                    db_scalar(row.r2_ci_lower_canonical_fraction),
                    db_scalar(row.r2_ci_upper_canonical_fraction),
                    db_scalar(row.r2_scale_reason),
                    db_scalar(row.quality_plot_role),
                    db_scalar(row.quality_context_warning),
                    as_bool(row.preferred_for_primary_quality),
                    PARSER_VERSION,
                )
            )

        print(
            "Metric sources to insert:",
            metric_long_current["source"].value_counts().to_dict(),
        )
        print(
            "Effect bases to insert:",
            pd.Series(
                [item[8] for item in metric_rows],
                dtype="object",
            ).value_counts().to_dict(),
        )

        if len(str(PARSER_VERSION)) > 50:
            raise ValueError(
                "PGS_METRIC_PARSER_VERSION exceeds VARCHAR(50): "
                f"{PARSER_VERSION!r}"
            )

        cur.execute("TRUNCATE TABLE data_libraries.pgs_performance_metric_values RESTART IDENTITY")
        execute_values(
            cur,
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
                metric_family,
                metric_model_scope,
                metric_model_scope_confidence,
                metric_model_scope_reason,
                r2_scale_basis,
                r2_estimate_canonical_fraction,
                r2_estimate_display_percent,
                r2_ci_lower_canonical_fraction,
                r2_ci_upper_canonical_fraction,
                r2_scale_reason,
                quality_plot_role,
                quality_context_warning,
                preferred_for_primary_quality,
                parser_version
            ) VALUES %s
            """,
            metric_rows,
            page_size=2000,
        )

    conn.commit()

    print("PGS report curation seed complete.")
    print(f"  Categories: {len(CATEGORY_CONFIG)}")
    print(f"  Consumer traits: {len(trait_rows)}")
    print(f"  Curated PGS mappings: {len(score_rows)}")
    print(f"  Evaluation outcome mappings: {len(outcome_rows)}")
    print(f"  Normalized performance metrics: {len(metric_rows)}")
    print(f"  PPM applicability contexts: {len(context_rows)}")

    if unresolved_targets:
        print(
            f"  Outcome targets needing review: {len(unresolved_targets)} "
            "(stored with target_report_trait_id = NULL)"
        )
        for item in unresolved_targets[:20]:
            print("   ", item)

except Exception:
    conn.rollback()
    raise
finally:
    conn.close()
