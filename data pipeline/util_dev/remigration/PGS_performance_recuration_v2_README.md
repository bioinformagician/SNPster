# SNPster PGS performance metadata — second-pass curation

## Why this pass exists

The first curation correctly normalized PGS Catalog performance metrics into a long table, but it treated a reported AUROC, C-index or R² largely as a property of the evaluated PGS. That is not always true.

A performance metric can describe:

1. the PGS alone;
2. the PGS plus basic adjustment variables such as age/sex/PCs;
3. the PGS plus clinical predictors such as BMI, glucose, blood pressure or family history;
4. a clinical/base/comparator model rather than the PGS;
5. the incremental contribution of the PGS after covariates;
6. a restricted or risk-enriched evaluation population.

The second-pass curation makes those distinctions explicit.

## Main additions

Each normalized metric row now receives:

- `metric_family`
- `metric_model_scope`
- `metric_model_scope_confidence`
- `metric_model_scope_reason`
- `covariate_class`
- `clinical_predictors_present`
- `clinical_predictor_terms`
- `population_scope_hint`
- `population_scope_confidence`
- `population_scope_reason`
- `r2_scale_basis`
- `r2_estimate_canonical_fraction`
- `r2_estimate_display_percent`
- `quality_plot_role`
- `quality_context_warning`
- `preferred_for_primary_quality`

## Model-scope categories

- `pgs_only`: explicit PGS/PRS-only or no-covariate wording in the metric itself.
- `pgs_only_likely`: no covariates reported, but no explicit PGS-only wording.
- `pgs_plus_basic_adjustment`: age/sex/ancestry PCs/array/batch/site and similar adjustment variables.
- `pgs_plus_clinical_predictors`: model includes non-PGS clinical predictors such as BMI, blood pressure, glucose, lipids, family history, smoking, renal disease, etc.
- `pgs_plus_covariates`: covariates are present but not confidently classifiable.
- `incremental_pgs_contribution`: incremental R², delta AUROC/C-index, or explicit “variance explained by PGS alone” wording.
- `clinical_or_comparator_model`: explicit clinical/base/covariate-only comparator wording.

## Preferred metric rule

For the main SNPster quality visualization, explicit PGS-specific metrics should outrank adjusted/full-model metrics when both are available for the same PPM and quality family.

Example: `PPM008934` contains:

- dedicated AUROC = `0.99354` for the adjusted/full model;
- `PGS AUROC (no covariates) = 0.47826`;
- total/full-model R² = `0.5357`;
- `PGS R2 (no covariates) = 0.00227`.

The new curation marks the PGS-only AUROC and PGS-only R² as preferred for the primary quality figure. The full-model values remain available as contextual evidence.

This situation is common: 1,189 PPM/metric-family combinations contain an explicit PGS-specific alternative that should supersede an adjusted/dedicated value for the primary PGS-quality display.

Among dedicated AUROC/C-index values >= 0.85, 94 have a different explicit PGS-specific metric in the same PPM. The median drop from the high dedicated value to the PGS-specific value is about 0.32 AUROC/C-index units.

## R² scale normalization

The previous report heuristic (`<=1` means fraction, `>1` means already percent) was insufficient because strings such as `ΔR2 (%) = 0.96` are explicitly percentage points even though the value is below 1.

The new curation uses the raw metric wording:

- explicit `%`, `percent`, `R2*100`, etc. -> percentage points;
- no marker and value within [-1, 1] -> fraction, displayed as ×100 percent;
- no marker and absolute value >1 but <=100 -> percentage points inferred from mathematical range, because an R² fraction cannot exceed 1;
- values outside plausible fraction/percentage ranges -> ambiguous and should not be auto-plotted.

Current counts across R²-family metric rows:

- fraction assumed: 20,224
- explicit percentage points: 316
- percentage points inferred from range: 55

## Evaluation-population scope

This pass uses the Performance Metrics sheet (`Reported Trait` + `Other Relevant Information`) to detect explicit restrictions such as:

- cases only;
- family-history/high-risk subgroups;
- disease-restricted cohorts;
- demographic subgroups.

This is intentionally incomplete. A full evaluation-population curation should also inspect the `Evaluation Sample Sets` sheet/PSS metadata because risk enrichment can be documented there rather than in the performance row.

Therefore `general_or_unspecified` means “no restriction detected in the performance-row text”, **not** “proven general-population validation.”

## Recommended report behavior

### Primary quality figure

Use only the preferred metric per PPM/family.

- PGS-only AUROC/C-index -> normal primary point.
- PGS-only standard R² -> primary bar.
- basic-adjusted metrics -> allowed, but explicitly labelled `adjusted`.
- clinical-predictor/full-model metrics -> visually secondary/contextual, not presented as PGS-only quality.
- clinical/comparator models -> exclude from primary PGS-quality figure.
- restricted/risk-enriched evaluations -> badge/warning on the metric.

### Evidence details

Keep every parsed metric visible, including full-model values, incremental metrics, correlation metrics, group contrasts and comparator metrics. The purpose is to contextualize rather than discard evidence.

## Remaining curation work before final DB/report integration

1. Curate the `Evaluation Sample Sets` sheet for population restrictions, recruitment enrichment, case/control counts and reliability flags.
2. Add a sample-size/reliability layer (e.g. very small case counts should not visually compete with well-powered validations).
3. Persist model-scope and scale fields in the normalized DB table or a linked context table.
4. Update the report to select `preferred_for_primary_quality = true` rather than dedicated AUROC/C-index/R² columns by default.
