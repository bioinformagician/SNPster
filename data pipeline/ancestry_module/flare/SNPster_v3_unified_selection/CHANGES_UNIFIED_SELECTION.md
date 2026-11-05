# Unified ancestry selection — change summary

## Files added
- `pedigree_benchmark/scripts/select_unified_ancestry.py` — one loss for reference + Ped-sim, macro-balancing by scenario, combination, family; p95 and 2–5% recovery safeguards.
- `pedigree_benchmark/scripts/select_stage2_diverse.py` — prevents Stage 2 from discarding every candidate from a surviving representation.
- `config/unified_selection.json` — thresholds and class requirements.
- `tests/test_unified_selection.py`, `tests/test_stage2_diverse.py` — unit tests.

## Files modified
- `main.nf`: Ped-sim executes after Stage 2 shortlist, scores *all Stage 3 candidates*, then the new `PDB_SELECT_UNIFIED_FINALISTS` process writes actual finalists. The old selector is used only when `PEDIGREE_ENABLED=false`.
- `nextflow.config`: resources for the unified selector.
- `run_benchmark_v3.sh`, `run_in_container_v3.sh`: unified branch enabled by default, 90 pedigree families, batch size 10, separate unified-config fingerprint.
- `scripts/config_fingerprint.sh`: excludes the new unified selection config so the legacy upstream cache fingerprint stays identical.
- `README_INTEGRATION.md`: instructions and limitations.

## What remains unchanged
- Model hyperparameter search JSON files.
- Existing expensive model construction/FLARE scripts and `config/params.json`.
- Decoder architecture and persisted Stage 3 `model.joblib` artifacts.
- Original reference-only leaderboard (kept for comparison; not a final selector).

## Important limitation
Stage 1 still prunes representations using reference-only validation. The main **final** selection is unified; the entire search is not yet fully unbiased to admixture. Also these samples bypass microarray/Beagle/QC and are not an independent untouched final test.
