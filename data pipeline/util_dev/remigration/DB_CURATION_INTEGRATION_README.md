# SNPster DB curation integration

This changes the architecture so **report generation has no CSV dependency**.
The curated CSVs are used only as ETL/seed inputs when the global metadata is
updated. Once seeded, the QMD and scoring setup read PostgreSQL only.

## Files

- `001_add_pgs_report_curation.sql`
  - one-time schema migration
  - adds consumer report taxonomy, outcome mappings and normalized metric tables
  - adds two DB views used by the pipeline/report
- `seed_pgs_report_curation.py`
  - imports the current curated score list, outcome mapping and normalized metric
    rows into PostgreSQL
- `refresh_pgs_metadata_db.py`
  - replacement for the existing PGS Catalog metadata loader
  - preserves `pgscatalog_data`, publications, ontology rows, user PRS results and
    report curation instead of `TRUNCATE ... CASCADE` on parent tables
  - refreshes replaceable child metadata tables
  - can refresh normalized metric rows from `pgs_performance_metrics_long.csv`
- `setup_pipeline_db_curated.py`
  - replacement for the existing scoring/setup script
  - `setup_pgs_reports()` now queries `data_libraries.report_score_catalog_v`
    rather than reading a categorized-score CSV
- `SNPster_full_report_db_only.qmd`
  - report version that reads the curated score taxonomy and structured
    performance metadata entirely from PostgreSQL

## One-time migration order

### 1. Back up PostgreSQL

The current metadata loader uses destructive truncation, so take a backup before
changing the refresh workflow.

### 2. Run the SQL migration

```bash
psql -h 127.0.0.1 -p 5433 -U postgres -d snpster_db \
  -f 001_add_pgs_report_curation.sql
```

### 3. Put the curation seed files in the metadata directory

By default the seed script looks in:

```text
/home/frederik/github_projects/SNPster/data pipeline/reporting_module/data/
```

for:

```text
pgs_reports_shop_curated.csv
evaluation_outcome_mapping.csv
pgs_performance_metrics_long.csv
pgs_metric_overview.csv
```

You can instead set `PGS_CURATION_DIR` or the individual path environment
variables documented in `seed_pgs_report_curation.py`.

### 4. Seed the curation tables

```bash
python seed_pgs_report_curation.py
```

The seed is transactional. It loads:

```text
report_categories
report_traits
report_trait_scores
evaluation_outcome_mappings
performance_metric_types
pgs_performance_metric_values
```

### 5. Replace the two existing update/setup scripts

Use:

```text
refresh_pgs_metadata_db.py
setup_pipeline_db_curated.py
```

instead of the previous metadata-loader and `setup_pgs_reports()` scripts.

### 6. Rebuild the operational PGS scoring registry

Running `setup_pipeline_db_curated.py` makes `setup_pgs_reports()` query the DB
curation and rebuild:

```text
snpster_users.pgs_reports_shop
```

`pgs_reports_shop` remains the operational table containing scoring-file paths;
consumer names/categories live in `data_libraries.report_*`.

## Normal update workflow afterwards

For a new PGS Catalog metadata release:

1. replace/update `pgs_all_metadata.xlsx`
2. run `refresh_pgs_metadata_db.py`
3. if the curation/metric parser output changed, regenerate the curated CSVs and
   rerun `seed_pgs_report_curation.py`
4. run `setup_pipeline_db_curated.py` if the curated score membership changed
5. render reports with `SNPster_full_report_db_only.qmd`

The **report itself never reads those CSVs**.

## Why the raw metadata loader changed

The old loader truncated `pgscatalog_data` with `CASCADE`. Because user PRS
results and the new curation tables reference PGS IDs, that would also remove
those dependent rows. The replacement loader therefore:

- upserts `pgscatalog_data`
- upserts `pgs_publications`
- upserts `ontology_mappings`
- truncates/reloads only replaceable child metadata such as performance rows,
  evaluation sample sets and score-development samples

This preserves user-specific scoring results and the manually curated taxonomy.

## Useful verification queries

```sql
SELECT COUNT(*) FROM data_libraries.report_categories;
SELECT COUNT(*) FROM data_libraries.report_traits WHERE is_active;
SELECT COUNT(*) FROM data_libraries.report_trait_scores WHERE include_in_report;
SELECT mapping_status, COUNT(*)
FROM data_libraries.evaluation_outcome_mappings
GROUP BY mapping_status
ORDER BY mapping_status;
SELECT COUNT(*) FROM data_libraries.pgs_performance_metric_values;
```

Inspect the consumer scoring registry:

```sql
SELECT category, reported_trait, COUNT(*) AS pgs_count
FROM data_libraries.report_score_catalog_v
GROUP BY category, reported_trait
ORDER BY category, reported_trait;
```

Inspect routed performance evidence:

```sql
SELECT
    ppm_id,
    pgs_id,
    reported_trait AS raw_evaluation_outcome,
    evaluation_outcome_clean,
    evaluation_target_trait,
    odds_ratio_estimate,
    odds_ratio_basis,
    odds_ratio_personalization_allowed,
    r2_value,
    incremental_r2_value
FROM data_libraries.pgs_performance_curated_v
WHERE evaluation_target_trait IS NOT NULL
LIMIT 100;
```
