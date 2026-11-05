SNPster PGS performance context v3 integration
==============================================

Purpose
-------
This update moves the v3 performance-metric interpretation context into PostgreSQL
and makes the Quarto report consume it directly.

The conservative applicability rule is intentional:
- Sample Age mean/range is stored as cohort context only. It does NOT by itself
  make an evaluation age-restricted.
- Explicit age restrictions (for example, "before age 75", "aged 40-54"),
  explicit disease/carrier/high-risk/case-only restrictions, and explicit sex
  restrictions set applicability_requires_match = TRUE.
- Until a user-context matcher can confirm the required age/condition/sex,
  restricted PPMs are retained for transparency but excluded from personalized
  OR/HR/Beta effects, primary AUROC/C-index quality, and main R2 evidence.
- Clinical/comparator models and PGS+clinical-predictor model metrics are not
  treated as PGS-only personalized/primary-quality evidence.

Files
-----
1. recurate_pgs_performance_context_v3_applicability.py
   Reads Performance Metrics + Evaluation Sample Sets from the PGS Catalog XLSX
   and the first-pass normalized metric CSV. Produces:
     - pgs_performance_context_curation_v3.csv
     - pgs_performance_metrics_long_context_v3.csv
     - pgs_performance_context_audit_v3.csv

2. 002_add_pgs_performance_context_v3.sql
   Migration for an existing SNPster database.

3. db_structure_v3.sql
   Full updated schema for a fresh database rebuild.

4. refresh_pgs_metadata_db_v3.py
   Refreshes PGS metadata and now retains PSS Sample Age, phenotype definitions,
   ancestry detail, follow-up and additional sample/cohort context.

5. seed_pgs_report_curation_v3.py
   Loads the v3 PPM context and metric-level semantics into PostgreSQL.

6. snpster_polygenic_score_report_v3.qmd
   Updated report that consumes the curated DB fields.

Recommended order for an existing database
------------------------------------------
A. Produce the v3 curation CSVs (after your existing first-pass metric parser):

   python recurate_pgs_performance_context_v3_applicability.py \
     --metadata /path/to/pgs_all_metadata.xlsx \
     --long /path/to/pgs_performance_metrics_long.csv \
     --output-dir /path/to/curation_dir

B. Apply the database migration:

   psql -h 127.0.0.1 -p 5433 -U postgres -d snpster_db \
     -f 002_add_pgs_performance_context_v3.sql

C. Refresh the PGS Catalog metadata using refresh_pgs_metadata_db_v3.py.
   This populates the newly retained Evaluation Sample Set context fields.

D. Run seed_pgs_report_curation_v3.py AFTER the metadata refresh.
   Point these environment variables at the v3 files if they are not in the
   script's default curation directory:

   PGS_METRIC_LONG_FILE=/path/to/pgs_performance_metrics_long_context_v3.csv
   PGS_PPM_CONTEXT_FILE=/path/to/pgs_performance_context_curation_v3.csv

   Important: refresh_pgs_metadata_db_v3.py reloads performance metadata, so if
   you rerun it later, rerun seed_pgs_report_curation_v3.py afterwards.

E. Render snpster_polygenic_score_report_v3.qmd so the DB queries execute against
   the newly seeded data.

Fresh database
--------------
Use db_structure_v3.sql instead of the migration, then run the metadata refresh
and v3 seed as above.

No change required
------------------
setup_pipeline_db_curated.py does not need a v3-specific change for this feature;
it already consumes the report catalog tables and does not interpret performance
metrics itself.
