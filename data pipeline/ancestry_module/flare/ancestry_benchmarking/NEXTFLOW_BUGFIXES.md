# Nextflow benchmark fixes and stage-2 recovery

Date: 2026-10-09

This document explains the workflow fixes and the recovery that allowed the
ancestry benchmark to finish. It separates confirmed failures from causes that
were not established.

## 1. Prevent copying a refinement config onto itself

**Problem:** `REFINE_MODEL` needs to publish a representation-specific config
alongside its final model. If the staged input already has the output filename,
the final `cp` can fail with an "are the same file" error, even after the expensive
model work has succeeded.

**Fix:** In [main.nf](./main.nf), stage the config input as
`refinement_input.config.json`, then copy it to `${rep_id}.config.json`.

**Effect:** Input and output filenames are distinct. The copy works both when
FLARE refinement runs and when refinement is disabled and the initial model is
copied directly to the final model.

## 2. Replace large candidate-directory staging with a manifest

**Observed failure:** The earlier run completed all **27,121 stage-2 fits**, then
`SELECT_STAGE2_CANDIDATES` exited with status **139** and no command output.
Status 139 indicates a segmentation fault; the exact native component responsible
was not established.

The old selection task staged thousands of candidate directories as individual
symlinks. Its generated `.command.run` was approximately **5 MB**.

**Fix:**

- [main.nf](./main.nf) collects result-directory paths into
  `stage2_result_dirs.txt` and stages that single manifest for selection.
- [select_candidates.py](./scripts/select_candidates.py) accepts `--manifest`
  as an alternative to `--root`.
- Each manifest entry is checked for `metrics.json`; a missing file raises an
  explicit error identifying the manifest line and path.

**Effect:** Selection reads the completed result directories directly rather than
requiring Nextflow to stage tens of thousands of directory inputs. The recovery
selection task succeeded with this approach. This verifies the workaround, not
the precise cause of the original segmentation fault.

## 3. Add a recovery entry point at stage-2 selection

**Problem:** Restarting with ordinary `-resume` did not continue directly at the
failed selection task.

The previous run had **55 successful initial-model builds and 18 failed builds**.
The failed builds were ignored, allowing the workflow to continue with successful
models. Nextflow reused the 55 successes on resume, but attempted the 18 failures
again, including expensive Bayesian mixture fits that had failed to converge.

This is expected Nextflow behavior: `-resume` reuses valid successful task caches;
it is not a "start at this process" instruction, and ignored failures are not
successful cached results.

**Fix:** Add the optional `--resume_stage2_manifest` parameter in
[main.nf](./main.nf). When supplied, the workflow:

1. Reads the existing stage-2 result manifest.
2. Loads the published `selection/selected_representations.tsv`.
3. Loads selected representations' final models and config metadata.
4. Loads their existing development feature tables.
5. Starts at `SELECT_STAGE2_CANDIDATES` and continues through diagnostics, full
   evaluation, aggregation and finalist selection.

Required files are resolved with `checkIfExists: true`; recovered model metadata
must match its representation ID. The ordinary workflow remains available when
the recovery parameter is absent.

**Effect:** Recovery bypasses initial model building, refinement, development
inference, representation screening, representation selection, candidate
generation and stage-2 fitting. It does not change convergence criteria or pretend
that failed models succeeded.

## 4. Preserve valid evaluation behavior

The related evaluation safeguards in
[fit_evaluate_candidate.py](./scripts/fit_evaluate_candidate.py) ensure that a
successful task does not silently produce incomplete or differently trained
results:

- Missing features for an enabled dataset fail explicitly.
- Missing truth data for the required selection role fail explicitly.
- Estimators that cannot accept required sample weights are rejected.
- Weighted fitting errors propagate instead of falling back to unweighted fits.

These safeguards protect result correctness; they are not a fix for the native
segmentation fault.

## 5. Understand fingerprint invalidation

[config_fingerprint.sh](./scripts/config_fingerprint.sh) hashes the search space,
top-level JSON/TSV configs and Python helper scripts. The fingerprint is included
in downstream task scripts so changes invalidate affected task caches.

The failed selection run used fingerprint `9301634...`; the recovery used
`1c5d747...`. This change can invalidate downstream screening, fitting, selection
and aggregation tasks in an ordinary resume. It was **not** the reason the
18 initial-model failures were retried: those tasks do not include this fingerprint.

The recovery entry point deliberately reuses existing stage-2 metrics. It is
appropriate for a selection/downstream fix with unchanged scientific inputs,
not for reusing fits after changing their datasets, training configuration or
evaluation criteria.

## 6. Recovery performed and verified

The recovery on 2026-10-09:

1. Reconstructed the exact original selection input set from the failed task's
   staging commands, rather than mixing published results from different sweeps.
2. Verified all **27,121** directories had nonempty `metrics.json` and
   `candidate.json`, and their parent fitting tasks had exit code **0**.
3. Saved the input paths in
   `results_v3/selection/recovery_stage2_result_dirs.txt`.
4. Previewed the recovery workflow and confirmed only selection and later
   processes were included.
5. Stopped the retrying container and retained it as
   `snpster-ancestry-v3-full-before-stage2-recovery`.
6. Launched the recovery as `snpster-ancestry-v3-full` with the same workspace,
   read-only reference mounts, user, 20-CPU limit and 180-GiB memory limit.

The earlier temporary-directory startup error was not reproduced during this
recovery. Its exact cause was not established; no permissions fix is claimed.

### Verified outcome

| Process | Result |
|---|---|
| `SELECT_STAGE2_CANDIDATES` | 1/1 succeeded |
| `RUN_REPRESENTATION_DIAGNOSTIC` | No tasks: all diagnostic datasets were disabled |
| `FULL_EVALUATE_CANDIDATE` | 17/17 succeeded |
| `AGGREGATE_FINAL` | 1/1 succeeded |
| `SELECT_FINALISTS` | 1/1 succeeded |

The container exited with **code 0**. Nextflow recorded **20 successful tasks,
0 failed tasks and 0 ignored tasks** in this recovery run. The final leaderboard
contains 17 candidates, and the finalist table contains 5 finalists.

Outputs are in [results_v3/final/](./results_v3/final/), including
[leaderboard.html](./results_v3/final/leaderboard.html),
[leaderboard.tsv](./results_v3/final/leaderboard.tsv) and
[finalists.tsv](./results_v3/final/finalists.tsv).

Both normal and recovery workflow paths were validated with Nextflow `-preview`.
The recovery path was also validated by the complete successful execution above.

## Recovery command

Inside the configured benchmark container, with `/workspace` as the working
directory:

```bash
fingerprint=$(bash scripts/config_fingerprint.sh config/search_space.json)

nextflow run main.nf \
  -resume 24289aff-d45c-4714-8cc4-ad4a3b3af817 \
  -params-file config/params.json \
  --search_space /workspace/config/search_space.json \
  --outdir /workspace/results_v3 \
  --config_fingerprint "$fingerprint" \
  --max_stage2_candidates 0 \
  --resume_stage2_manifest \
    /workspace/results_v3/selection/recovery_stage2_result_dirs.txt
```

Stop any existing run in the same workspace first. Keep the manifest, original
work directories, published models and development features intact. Do not
regenerate model configs or use the model-generating wrapper for this recovery.
