# SNPster v3 — Unified ancestry-proportion selection (Ped-sim integration)

## Purpose

This package **replaces reference-only final model selection** with one objective: accuracy of the eight predicted ancestry proportions against truth, regardless of whether a person is reference-labelled 100% EUR, an F1 50/50, or a realistic multi-generation admixed descendant.

The TV loss used for every sample is `0.5 * sum(abs(predicted_source - truth_source))` over `[EUR,EAS,SAS,SEA,GME,NAM,OCE,SSA]`. Source-reference samples use source-label truth; Ped-sim offspring use recombination-derived **realized founder-label inheritance**, not nominal pedigree expectations. These are operational truth definitions, not proof that real reference people have biologically pure ancestry.

**Important:** this version uses high-density VCFs projected directly to the FLARE panel. It **does not** run consumer-array masking, harmonization, Beagle imputation or production QC. Therefore it does not reproduce the entire Lonni domain-shift problem yet.

## Selection architecture

```text
Unchanged costly v3 FLARE representation building
       |
Stage 1: existing reference-based coarse representation screen
       |
Stage 2: decoder sweep; broaden shortlist to include at least one
         candidate from every retained representation
       |
Stage 3: existing full validation; save model.joblib, predictions.tsv
       |                               |
       |          multi-generation Ped-sim descendants from held-out donors
       |                               |
       |                     FLARE + each stage-3 saved decoder
       |                               |
       +-------------------------------+
                       |
       UNIFIED SELECTION (NEW)
       - Same eight-component TV distance for all samples
       - Balanced over scenario, ancestry combination, family
       - Stratified p95 TV and off-target mass safety gates
       - Rare 2%-5% ancestry detection-recall diagnostic/gate
       - Keep Lonni regression gate, no separate pure-accuracy objective
                       |
           finalists.tsv + unified_leaderboard.tsv
```

The former `SELECT_FINALISTS` process is executed **only when `PEDIGREE_ENABLED=false`**. The original `AGGREGATE_FINAL` process still writes the reference-only leaderboard for comparison, but it no longer controls final selection when unified mode is enabled. The same model artifacts already fitted by v3 are reused; Ped-sim samples are **not** decoder training data.

## Group balancing

To stop 100 siblings from outweighing a separate pedigree, the selector averages per-sample TV within each family, averages families within each ancestry combination, then averages combinations within each scenario. Scenarios (the reference-source group plus the nine Ped-sim classes) have equal aggregate importance, rather than raw sample-count weighting. Reference-source combinations are the eight source labels.

`balanced_mean_tv` is the primary ranking variable. `balanced_p95_tv` is an **equal-scenario average of per-scenario family/combination-weighted 95th percentiles**, not the empirical 95th percentile of all pooled individuals. `balanced_off_target_mass` and optional Lonni/low-level constraints are gates, not alternative definitions of ancestry accuracy. The configuration defaults are **provisional engineering gates**, not research-validated thresholds; change them with evidence after viewing the first results.

For 2–5% realized ancestry components, recall is the proportion with predicted fraction >=1%. A model predicting 100% EUR for a true 96.875% EUR / 3.125% SEA descendant has TV=0.03125 but recall=0 for that SEA component. The gate is applied once at least 20 eligible component observations are available. If support is smaller, the diagnostic includes `minor_gate_not_evaluable=true`; it is **not** evidence of successful rare-ancestry detection.

## How to install

For this server deployment, follow [the new sweep deployment notes](./README_NEW_SWEEP.md#deploy):
Ped-sim is bundled in the derived Docker image, the existing autosomal GRCh38
map is used, and founders must be training/reference-disjoint with validated
fine/source labels. The five historical supplemental founders are excluded.
The default is the separately labeled reduced design because only two eligible
OCE donors remain. The instructions below describe the original package
integration approach, not the current separate-directory deployment.

1. **Do not overwrite files during the actively running original v3 process.** Finish or freeze that run first. Back up `main.nf`, `nextflow.config`, launchers, `config/`, `results_v3`, and `.nextflow` metadata.
2. Unzip the package and copy the `SNPster_v3_unified_selection` contents into the existing `ancestry_benchmarking` project directory, preserving your existing live `results_v3`, `work`, and `.nextflow` directories.
3. Ensure `/srv/dependencies/ancestry/ped-sim/ped-sim` is executable. The Docker launcher mounts `/srv/dependencies/ancestry/hgdp_1kg_references` and `/srv/dependencies/imputation_runner/imputer` read-only.
4. Founder preflight validates the actual requested trees before costly training. It excludes both training-mosaic donors and reference donors and checks all source/fine labels against `ancestry_role_audit.tsv`. `PEDIGREE_DESIGN=full` requires additional independent donors; `PEDIGREE_DESIGN=reduced` rejects infeasible trees and records the changed coverage.
5. The costly upstream model-building scripts are unchanged. Use recovery with `RECOVERY_SOURCE_OUTDIR` pointing to existing models/features and a fresh `RUN_OUTDIR` to avoid both rebuilding models and overwriting confounded historical pedigree scores. Configuration changes can otherwise invalidate upstream caches.

Smoke test (**reduced chromosomes are not a genome-wide precision benchmark**):

```bash
MAX_MODELS=0 MAX_STAGE2_CANDIDATES=3000 \
PEDIGREE_ENABLED=true PEDIGREE_FAMILIES=9 \
PEDIGREE_BATCH_SIZE=3 PEDIGREE_CHROMOSOMES=1,2 \
./run_in_container_v3.sh
```

Full run:

```bash
MAX_MODELS=0 MAX_STAGE2_CANDIDATES=3000 \
PEDIGREE_ENABLED=true PEDIGREE_FAMILIES=90 \
PEDIGREE_BATCH_SIZE=10 \
PEDIGREE_CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
./run_in_container_v3.sh
```

`MAX_STAGE2_CANDIDATES=3000` is a **candidate-search cap**, not the number of finalists. Keep your existing stage-2 settings if resuming a finished search and minimizing rework is more important than expanding the search.

Legacy comparison:

```bash
PEDIGREE_ENABLED=false MAX_MODELS=0 MAX_STAGE2_CANDIDATES=3000 ./run_in_container_v3.sh
```

Environment variables: `PEDIGREE_FAMILIES` (default 90), `PEDIGREE_BATCH_SIZE` (default 10), `PEDIGREE_SEED` (20261008), `PEDIGREE_DESIGN` (default `reduced`), `PEDIGREE_ALLOW_REUSE`, `PEDIGREE_KEEP_HIGH_DENSITY`, and `PEDIGREE_CHROMOSOMES`. `PEDIGREE_FAMILIES` must be at least 9 in unified mode to include every required pedigree class. The generator cycles through nine classes. Reduced and full designs must not share a leaderboard. The default output is `results_v3_sweep_100_disjoint_reduced/`; substitute the selected `RUN_OUTDIR` for the historical `results_v3/` examples below.

## Outputs

- `results_v3/final/unified_leaderboard.tsv` — all stage-3 candidates ranked by balanced mean TV with gates, scenario count, sample counts and minor-ancestry recall.
- `results_v3/final/unified_diagnostics.json` — per-candidate per-scenario / ancestry-combination metrics, family counts, per-source MAE, gate outcomes, rare-ancestry support.
- `results_v3/final/finalists.tsv` and `finalists/*.json` — models chosen by **unified** accuracy. If no candidate passes gates, selection fails instead of silently promoting one.
- `results_v3/final/leaderboard.tsv` — old reference-only comparison, **not** the new selection leaderboard.
- `results_v3/pedigree/evaluations/*.pedigree_predictions.tsv` — per-person predictions and realized ancestry fractions.
- `results_v3/pedigree/truth/all_samples_ancestry.tsv` — realized inheritance vs pedigree expectation.
- `results_v3/pedigree/inputs/*.summary.json` — source availability and FLARE reference-overlap audit.

## Limitations to resolve for a scientific final benchmark

1. **Early representation screening remains reference-driven.** Retaining a candidate per surviving representation helps prevent stage-2 pruning, but a representation lost in Stage 1 can still be superior on admixture. A fully selection-unbiased architecture would also include development Ped-sim data in Stage 1 or retain every successful representation in that stage, which would substantially expand the stage-2 search.
2. Ped-sim families can reuse held-out donors across pedigrees if `PEDIGREE_ALLOW_REUSE=true` (the default); these are therefore **development/selection simulations, not an independent untouched final test**. A strict final test needs new source/founder-disjoint pools and preferably fine-population splits.
3. Founders are drawn from `source_validation_metadata.tsv`; the purity of their broad source label is an assumption, not independently verified ancestry.
4. Selection is currently from **direct high-density genotypes**, not SNPster's production microarray -> Beagle -> QC domain. Production-domain validation should ultimately become the required selection set.
5. The corrected founder preflight, regression tests and offline Nextflow preview have passed on the server. Controlled genome-wide founder inference confirms the reference-overlap failure. A corrected full offspring benchmark has not yet been run; its performance cannot be inferred from the historical scores.

## The most important files changed

- `main.nf` — new native `PDB_SELECT_UNIFIED_FINALISTS` DAG stage and pre-finalist Ped-sim branch.
- `pedigree_benchmark/scripts/select_unified_ancestry.py` — unified metrics, balancing, gates and finalist selection.
- `pedigree_benchmark/scripts/select_stage2_diverse.py` — ensure each retained representation reaches final Ped-sim evaluation.
- `config/unified_selection.json` — versioned scoring/gate config.
- `nextflow.config`, `run_benchmark_v3.sh`, `run_in_container_v3.sh`, `scripts/config_fingerprint.sh` — resources, activation, fingerprints, caching.

All other previously supplied v3 and pedigree files are included so this is a complete project overlay.
