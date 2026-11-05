BEGIN;

-- -----------------------------------------------------------------------------
-- 1) Preserve evaluation-sample context that was previously discarded.
--    These fields are descriptive. Sample-age mean/range is NOT itself treated
--    as a user eligibility restriction.
-- -----------------------------------------------------------------------------
ALTER TABLE data_libraries.evaluation_sample_sets
    ADD COLUMN IF NOT EXISTS sample_age_raw TEXT,
    ADD COLUMN IF NOT EXISTS ancestry_detail TEXT,
    ADD COLUMN IF NOT EXISTS additional_ancestry_description TEXT,
    ADD COLUMN IF NOT EXISTS phenotype_definitions_methods TEXT,
    ADD COLUMN IF NOT EXISTS followup_time TEXT,
    ADD COLUMN IF NOT EXISTS additional_sample_cohort_information TEXT;

-- -----------------------------------------------------------------------------
-- 2) PPM-level interpretation/applicability context.
--    applicability_requires_match = TRUE means the PPM must not be used as
--    general personalized evidence unless the user's relevant context is known
--    and explicitly matches the restriction.
-- -----------------------------------------------------------------------------
ALTER TABLE data_libraries.pgs_performance
    ADD COLUMN IF NOT EXISTS has_covariates BOOLEAN,
    ADD COLUMN IF NOT EXISTS covariate_class VARCHAR(100),
    ADD COLUMN IF NOT EXISTS clinical_predictors_present BOOLEAN,
    ADD COLUMN IF NOT EXISTS clinical_predictor_terms TEXT,
    ADD COLUMN IF NOT EXISTS population_scope_hint VARCHAR(100),
    ADD COLUMN IF NOT EXISTS population_scope_confidence VARCHAR(20),
    ADD COLUMN IF NOT EXISTS population_scope_reason TEXT,
    ADD COLUMN IF NOT EXISTS applicability_status VARCHAR(100)
        NOT NULL DEFAULT 'general_or_unspecified',
    ADD COLUMN IF NOT EXISTS applicability_requires_match BOOLEAN
        NOT NULL DEFAULT FALSE,
    ADD COLUMN IF NOT EXISTS applicability_restriction_types TEXT,
    ADD COLUMN IF NOT EXISTS applicability_reason TEXT,
    ADD COLUMN IF NOT EXISTS age_restriction_min_years NUMERIC,
    ADD COLUMN IF NOT EXISTS age_restriction_max_years NUMERIC,
    ADD COLUMN IF NOT EXISTS sex_restriction VARCHAR(20),
    ADD COLUMN IF NOT EXISTS condition_restriction TEXT,
    ADD COLUMN IF NOT EXISTS sample_age_raw TEXT,
    ADD COLUMN IF NOT EXISTS sample_age_mean_years NUMERIC,
    ADD COLUMN IF NOT EXISTS sample_age_min_years NUMERIC,
    ADD COLUMN IF NOT EXISTS sample_age_max_years NUMERIC,
    ADD COLUMN IF NOT EXISTS phenotype_definition_raw TEXT,
    ADD COLUMN IF NOT EXISTS additional_sample_info_raw TEXT;

CREATE INDEX IF NOT EXISTS idx_pgs_performance_applicability
    ON data_libraries.pgs_performance (applicability_requires_match, ppm_id);

-- -----------------------------------------------------------------------------
-- 3) Metric-level v3 semantics.
-- -----------------------------------------------------------------------------
ALTER TABLE data_libraries.pgs_performance_metric_values
    ADD COLUMN IF NOT EXISTS metric_family VARCHAR(100),
    ADD COLUMN IF NOT EXISTS metric_model_scope VARCHAR(100),
    ADD COLUMN IF NOT EXISTS metric_model_scope_confidence VARCHAR(20),
    ADD COLUMN IF NOT EXISTS metric_model_scope_reason TEXT,
    ADD COLUMN IF NOT EXISTS r2_scale_basis VARCHAR(100),
    ADD COLUMN IF NOT EXISTS r2_estimate_canonical_fraction NUMERIC,
    ADD COLUMN IF NOT EXISTS r2_estimate_display_percent NUMERIC,
    ADD COLUMN IF NOT EXISTS r2_ci_lower_canonical_fraction NUMERIC,
    ADD COLUMN IF NOT EXISTS r2_ci_upper_canonical_fraction NUMERIC,
    ADD COLUMN IF NOT EXISTS r2_scale_reason TEXT,
    ADD COLUMN IF NOT EXISTS quality_plot_role VARCHAR(150),
    ADD COLUMN IF NOT EXISTS quality_context_warning TEXT,
    ADD COLUMN IF NOT EXISTS preferred_for_primary_quality BOOLEAN
        NOT NULL DEFAULT FALSE;

CREATE INDEX IF NOT EXISTS idx_metric_values_primary_quality
    ON data_libraries.pgs_performance_metric_values
        (ppm_id, metric_type, preferred_for_primary_quality)
    WHERE preferred_for_primary_quality;

-- -----------------------------------------------------------------------------
-- 4) Rebuild the report-facing view.
--    Existing columns stay in their original order; new v3 columns are appended
--    so downstream code that references the old columns remains compatible.
-- -----------------------------------------------------------------------------
CREATE OR REPLACE VIEW data_libraries.pgs_performance_curated_v AS
SELECT
    pp.ppm_id,
    pp.pgs_id,
    pp.pss_id,
    pp.pgp_id,
    pp.reported_trait,
    pp.covariates_included_in_model,
    pp.pgs_performance_other_relevant_info,
    pp.publication_pmid,
    pp.publication_doi,
    pp.hazard_ratio,
    pp.hazard_ratio_ci_lower AS hazard_ratio_ci_lower_dedicated,
    pp.hazard_ratio_ci_upper AS hazard_ratio_ci_upper_dedicated,
    pp.odds_ratio,
    pp.odds_ratio_ci_lower AS odds_ratio_ci_lower_dedicated,
    pp.odds_ratio_ci_upper AS odds_ratio_ci_upper_dedicated,
    pp.beta,
    pp.beta_ci_lower AS beta_ci_lower_dedicated,
    pp.beta_ci_upper AS beta_ci_upper_dedicated,
    pp.auroc,
    pp.auroc_ci_lower,
    pp.auroc_ci_upper,
    pp.concordance_statistic,
    pp.concordance_statistic_ci_lower,
    pp.concordance_statistic_ci_upper,
    pp.other_metric,

    COALESCE(eom.cleaned_outcome, pp.reported_trait) AS evaluation_outcome_clean,
    eom.outcome_key AS evaluation_outcome_key,
    eom.mapping_status AS evaluation_mapping_status,
    target_trait.report_trait_id AS evaluation_target_report_trait_id,
    target_trait.display_name AS evaluation_target_trait,
    target_category.display_name AS evaluation_target_category,

    metrics.odds_ratio_estimate,
    metrics.odds_ratio_ci_lower,
    metrics.odds_ratio_ci_upper,
    metrics.odds_ratio_basis,
    metrics.odds_ratio_personalization_allowed,

    metrics.hazard_ratio_estimate,
    metrics.hazard_ratio_ci_lower,
    metrics.hazard_ratio_ci_upper,
    metrics.hazard_ratio_basis,
    metrics.hazard_ratio_personalization_allowed,

    metrics.beta_estimate,
    metrics.beta_ci_lower,
    metrics.beta_ci_upper,
    metrics.beta_basis,

    metrics.auroc_estimate,
    metrics.c_index_estimate,

    metrics.r2_value,
    metrics.incremental_r2_value,
    metrics.partial_r2_value,
    metrics.pseudo_r2_value,
    metrics.adjusted_r2_value,
    metrics.liability_r2_value,

    metrics.group_odds_ratio_raw,
    metrics.group_hazard_ratio_raw,
    metrics.relative_risk_raw,
    metrics.relative_risk_ratio_raw,
    metrics.incidence_rate_ratio_raw,

    -- v3 PPM context
    pp.has_covariates,
    pp.covariate_class,
    pp.clinical_predictors_present,
    pp.clinical_predictor_terms,
    pp.population_scope_hint,
    pp.population_scope_confidence,
    pp.population_scope_reason,
    pp.applicability_status,
    pp.applicability_requires_match,
    pp.applicability_restriction_types,
    pp.applicability_reason,
    pp.age_restriction_min_years,
    pp.age_restriction_max_years,
    pp.sex_restriction,
    pp.condition_restriction,
    pp.sample_age_raw,
    pp.sample_age_mean_years,
    pp.sample_age_min_years,
    pp.sample_age_max_years,
    pp.phenotype_definition_raw,
    pp.additional_sample_info_raw,

    -- v3 quality values: only preferred PGS-quality metrics from PPMs that do
    -- not require an unresolved age/condition/sex match.
    metrics.primary_quality_auroc,
    metrics.primary_quality_c_index,
    metrics.beta_personalization_allowed

FROM data_libraries.pgs_performance AS pp
LEFT JOIN data_libraries.evaluation_outcome_mappings AS eom
    ON eom.raw_outcome = pp.reported_trait
LEFT JOIN data_libraries.report_traits AS target_trait
    ON target_trait.report_trait_id = eom.target_report_trait_id
LEFT JOIN data_libraries.report_categories AS target_category
    ON target_category.category_id = target_trait.category_id
LEFT JOIN LATERAL (
    SELECT
        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type = 'odds_ratio'
              AND mv.source = 'dedicated'
        ) AS odds_ratio_estimate,
        MAX(mv.ci_lower) FILTER (
            WHERE mv.metric_type = 'odds_ratio'
              AND mv.source = 'dedicated'
        ) AS odds_ratio_ci_lower,
        MAX(mv.ci_upper) FILTER (
            WHERE mv.metric_type = 'odds_ratio'
              AND mv.source = 'dedicated'
        ) AS odds_ratio_ci_upper,
        MAX(mv.effect_basis) FILTER (
            WHERE mv.metric_type = 'odds_ratio'
              AND mv.source = 'dedicated'
        ) AS odds_ratio_basis,
        COALESCE(
            BOOL_OR(mv.personalization_allowed) FILTER (
                WHERE mv.metric_type = 'odds_ratio'
                  AND mv.source = 'dedicated'
                  AND COALESCE(mv.metric_model_scope, '') NOT IN ('clinical_or_comparator_model', 'pgs_plus_clinical_predictors')
                  AND NOT COALESCE(pp.applicability_requires_match, FALSE)
            ),
            FALSE
        ) AS odds_ratio_personalization_allowed,

        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type = 'hazard_ratio'
              AND mv.source = 'dedicated'
        ) AS hazard_ratio_estimate,
        MAX(mv.ci_lower) FILTER (
            WHERE mv.metric_type = 'hazard_ratio'
              AND mv.source = 'dedicated'
        ) AS hazard_ratio_ci_lower,
        MAX(mv.ci_upper) FILTER (
            WHERE mv.metric_type = 'hazard_ratio'
              AND mv.source = 'dedicated'
        ) AS hazard_ratio_ci_upper,
        MAX(mv.effect_basis) FILTER (
            WHERE mv.metric_type = 'hazard_ratio'
              AND mv.source = 'dedicated'
        ) AS hazard_ratio_basis,
        COALESCE(
            BOOL_OR(mv.personalization_allowed) FILTER (
                WHERE mv.metric_type = 'hazard_ratio'
                  AND mv.source = 'dedicated'
                  AND COALESCE(mv.metric_model_scope, '') NOT IN ('clinical_or_comparator_model', 'pgs_plus_clinical_predictors')
                  AND NOT COALESCE(pp.applicability_requires_match, FALSE)
            ),
            FALSE
        ) AS hazard_ratio_personalization_allowed,

        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type = 'beta'
              AND mv.source = 'dedicated'
        ) AS beta_estimate,
        MAX(mv.ci_lower) FILTER (
            WHERE mv.metric_type = 'beta'
              AND mv.source = 'dedicated'
        ) AS beta_ci_lower,
        MAX(mv.ci_upper) FILTER (
            WHERE mv.metric_type = 'beta'
              AND mv.source = 'dedicated'
        ) AS beta_ci_upper,
        MAX(mv.effect_basis) FILTER (
            WHERE mv.metric_type = 'beta'
              AND mv.source = 'dedicated'
        ) AS beta_basis,
        COALESCE(
            BOOL_OR(TRUE) FILTER (
                WHERE mv.metric_type = 'beta'
                  AND mv.source = 'dedicated'
                  AND COALESCE(mv.metric_model_scope, '') NOT IN ('clinical_or_comparator_model', 'pgs_plus_clinical_predictors')
                  AND NOT COALESCE(pp.applicability_requires_match, FALSE)
            ),
            FALSE
        ) AS beta_personalization_allowed,

        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type = 'auroc'
              AND mv.source = 'dedicated'
        ) AS auroc_estimate,
        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type = 'c_index'
              AND mv.source = 'dedicated'
        ) AS c_index_estimate,

        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'r2') AS r2_value,
        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'incremental_r2') AS incremental_r2_value,
        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'partial_r2') AS partial_r2_value,
        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'pseudo_r2') AS pseudo_r2_value,
        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'adjusted_r2') AS adjusted_r2_value,
        MAX(mv.estimate) FILTER (WHERE mv.metric_type = 'liability_r2') AS liability_r2_value,

        STRING_AGG(DISTINCT mv.raw_metric, ' || ') FILTER (
            WHERE mv.metric_type = 'group_odds_ratio'
        ) AS group_odds_ratio_raw,
        STRING_AGG(DISTINCT mv.raw_metric, ' || ') FILTER (
            WHERE mv.metric_type = 'group_hazard_ratio'
        ) AS group_hazard_ratio_raw,
        STRING_AGG(DISTINCT mv.raw_metric, ' || ') FILTER (
            WHERE mv.metric_type = 'relative_risk'
        ) AS relative_risk_raw,
        STRING_AGG(DISTINCT mv.raw_metric, ' || ') FILTER (
            WHERE mv.metric_type = 'relative_risk_ratio'
        ) AS relative_risk_ratio_raw,
        STRING_AGG(DISTINCT mv.raw_metric, ' || ') FILTER (
            WHERE mv.metric_type = 'incidence_rate_ratio'
        ) AS incidence_rate_ratio_raw,

        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type IN ('auroc', 'auroc_other')
              AND mv.preferred_for_primary_quality
              AND mv.quality_plot_role LIKE 'primary_pgs_quality%'
              AND NOT COALESCE(pp.applicability_requires_match, FALSE)
        ) AS primary_quality_auroc,

        MAX(mv.estimate) FILTER (
            WHERE mv.metric_type IN ('c_index', 'cindex_other')
              AND mv.preferred_for_primary_quality
              AND mv.quality_plot_role LIKE 'primary_pgs_quality%'
              AND NOT COALESCE(pp.applicability_requires_match, FALSE)
        ) AS primary_quality_c_index

    FROM data_libraries.pgs_performance_metric_values AS mv
    WHERE mv.ppm_id = pp.ppm_id
) AS metrics ON TRUE;

COMMIT;
