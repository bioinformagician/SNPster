# SNPster ancestry strategy laboratory v3

This workflow is the high-throughput ancestry-model development harness for SNPster.
It intentionally starts **before** the current frozen FLARE K18 model: it creates many
candidate FLARE2 latent models from the same copying profiles, refines them through
FLARE, screens their latent representations, performs a large decoder/calibration
search, runs the common validation battery, and produces strategy-aware leaderboards.

There is **no ADMIXTURE execution** in this workflow.

## Why v3 exists

The frozen v1 model performed well on clean synthetic/reference-domain validation but
failed an important real SNPster case: Lonni, who is approximately first-generation
European/Filipino, acquired a large spurious SAS component. Lonni is therefore retained
as a **real-world regression gate**, while not being used as supervised training truth.

v3 broadens the search upstream instead of only tuning the ExtraTrees decoder.

## Model-building strategies

The default `config/search_space.json` searches 76 newly trained FLARE models plus the
exact frozen v1 K18 baseline. Candidate initial FLARE models are created using:

- diagonal Gaussian mixture models;
- tied-covariance Gaussian mixture models;
- spherical Gaussian mixture models;
- Bayesian Gaussian mixtures with a Dirichlet-process prior;
- Bayesian Gaussian mixtures with a finite Dirichlet-distribution prior;
- K-means clustering;
- PCA + Gaussian-mixture clustering;
- PCA + K-means clustering;
- selected unrefined GMM controls;
- selected `update-p=true, em=false` refinement controls.

All surviving models are then evaluated through the same frozen-inference path
(`update-p=false`, `em=false`) on development/validation users.

The optional `config/search_space_big.json` expands this further to roughly 522 FLARE
models, including MiniBatch K-means, more K values, two panel-cm settings, additional
Bayesian priors, more PCA dimensions, and additional seeds. It is intended for an
unattended large search, not the first smoke run.

## Decoder / calibration strategies

For retained FLARE representations, stage 2 can test:

- ExtraTrees;
- Random Forest;
- Histogram Gradient Boosting;
- MLP;
- Ridge;
- source-centroid NNLS;
- no calibration;
- Ridge calibration;
- multiple pure/mosaic training weights.

The exact frozen v1 ExtraTrees + Ridge artifacts are always carried as a regression
baseline.

## Successive search design

The workflow avoids blindly running the full Cartesian product:

1. **PANEL_PROBS** – FLARE copying profiles are generated once per panel configuration.
2. **BUILD_INITIAL_MODEL** – clustering strategy/hyperparameters create an initial FLARE model.
3. **REFINE_MODEL** – FLARE optionally refines the model with `update-p=true`.
4. **RUN_REPRESENTATION_DEV** – each final latent model is run on the shared development datasets.
5. **Stage 1 screen** – a common ExtraTrees + Ridge decoder ranks latent representations.
6. **Diversity retention** – at least the best representation from each model-building strategy is preserved before filling remaining slots by overall score.
7. **Stage 2 decoder search** – many decoder/leaf/tree/calibration/weight combinations are tested.
8. **Stage 3 full evaluation** – shortlisted candidates are re-evaluated and their full model artifacts are saved.
9. **Finalists** – the best candidates are written to `results_v3/final/`.

Stage 2 result directories are passed to the shortlist selector through a single manifest
file so large searches do not generate oversized Nextflow input-staging scripts.

One pathological GMM/Bayesian/FLARE candidate is allowed to fail independently; it does
not terminate the whole search.

## Current validation datasets

Enabled now:

- 307 pure-source training individuals;
- 1,008 actual-FLARE training mosaics;
- 159 untouched pure-source validation individuals;
- Lonni real-user regression benchmark.

Present but deliberately disabled until their metadata adapters are normalized:

- matched low-level ancestry benchmark;
- same-source null benchmark.

The previously inspected 672-genome synthetic final test remains disabled for tuning.
It should not silently become model-selection data.

## Lonni regression gate

Lonni is currently treated as approximately European + Filipino. The default hard gate
is deliberately broad rather than fitting her to an exact 50/50 number:

- EUR: 35–65%
- EAS + SEA: 35–65%
- SAS: <=15%
- GME + NAM + OCE + SSA: <=12%

This catches the known v1 failure without turning one individual into the objective
function.

## Directory layout

```text
ancestry_benchmarking_v3/
├── main.nf
├── nextflow.config
├── run_benchmark_v3.sh
├── run_in_container_v3.sh
├── requirements.txt
├── config/
│   ├── search_space.json
│   ├── search_space_big.json
│   ├── params.json
│   ├── dev_datasets.tsv
│   ├── diagnostic_datasets.tsv
│   ├── full_datasets.tsv
│   ├── fixed_models.tsv
│   ├── real_benchmarks.tsv
│   └── stage*_criteria.json
└── scripts/
    ├── create_model_file_strategy.py
    ├── generate_model_configs.py
    ├── generate_decoder_candidates.py
    ├── fit_evaluate_candidate.py
    ├── select_representations.py
    ├── select_candidates.py
    ├── aggregate_results_v3.py
    ├── estimate_search_size.py
    ├── extract_global_anc.py
    ├── make_screen_candidate.py
    └── preflight_v3.py
```

## Docker

Reuse the image that already worked for the previous benchmark:

```text
snpster-ancestry-benchmark:v1
```

The benchmark folder is bind-mounted to `/workspace` and the existing ancestry data are
mounted read-only from `/srv/dependencies/ancestry/flare2_test` and `/srv/imputed`.
No new production image is required.

## Inspect search size

Default broad search:

```bash
SEARCH_SPACE=config/search_space.json python3 scripts/estimate_search_size.py
```

Large strategy search:

```bash
SEARCH_SPACE=config/search_space_big.json python3 scripts/estimate_search_size.py
```

The large profile contains about 522 trainable FLARE models. Its full decoder Cartesian
product is intentionally enormous, so use `MAX_STAGE2_CANDIDATES` as a deterministic
subsample cap unless you explicitly want the entire decoder grid.

## Recommended first smoke run

```bash
./run_in_container_v3.sh
```

Defaults:

```text
MAX_MODELS=8
MAX_STAGE2_CANDIDATES=80
SEARCH_SPACE=config/search_space.json
```

The model cap is strategy-aware, so the smoke run spans several distinct model-building
approaches rather than taking the first eight ordinary GMM configurations.

## Full default strategy run

```bash
MAX_MODELS=0 \
MAX_STAGE2_CANDIDATES=3000 \
./run_in_container_v3.sh
```

This trains all 76 default FLARE representations, then samples 3,000 decoder/calibration
candidates across the retained representations. Because selection is deterministic and
round-robin by representation, the cap remains broad rather than being dominated by one
representation.

If compute is not a concern and you truly want the complete default stage-2 Cartesian
product:

```bash
MAX_MODELS=0 \
MAX_STAGE2_CANDIDATES=0 \
./run_in_container_v3.sh
```

That is roughly 30,000 stage-2 candidates and will be substantially slower.

## Resuming after a workflow fix

See [Nextflow fixes and stage-2 recovery](./NEXTFLOW_BUGFIXES.md) for the failure
history, implemented fixes and verified recovery results.

Stop the existing Nextflow run before starting another run in the same workspace.
Reuse its exact `nextflow run main.nf` command with `-resume` (or
`-resume <run-name>`), keeping the same parameters and generated model configs.
The launch command and run name are recorded near the start of `.nextflow.log`.
Avoid regenerating the configs for a script-only fix: changed file timestamps can
invalidate otherwise reusable task caches.

Nextflow reuses successful tasks whose inputs and process definitions are unchanged.
Failed or interrupted tasks are retried; changing a process definition also
invalidates its previous task hashes. In particular, non-converged Bayesian mixture
builds will be retried too, without changing their convergence settings.

For a helper-script-only resume, compute the downstream fingerprint without
regenerating model configs:

```bash
CONFIG_FINGERPRINT="$(bash scripts/config_fingerprint.sh config/search_space.json)"
```

Pass it as `--config_fingerprint "$CONFIG_FINGERPRINT"` in the original Nextflow
command. The fingerprint includes config files and Python helpers, invalidating
screening, candidate fitting, selection and aggregation without changing upstream
FLARE/model-building task definitions.

### Continue directly from completed stage-2 fits

If stage-2 fitting finished and candidate selection failed, use
`--resume_stage2_manifest /workspace/path/to/stage2_result_dirs.txt` with the
original launch command. The manifest must contain one absolute, container-visible
result-directory path per line, each containing `metrics.json` and `candidate.json`.
Use the exact result set from the failed selection task, not a mix of prior sweeps.

This mode bypasses model building, refinement, development inference, representation
screening and stage-2 fitting, including previously ignored model failures. It starts
at `SELECT_STAGE2_CANDIDATES`, then runs diagnostics, full evaluation and final
aggregation. It requires the original `--outdir` with
`selection/selected_representations.tsv`, the selected representations' final models
and config metadata under `models/`, and their development feature tables under
`features/dev/`. Missing required artifacts fail explicitly.

Keep the original search space and dataset/criteria configs: this is a recovery
entry point for a selection/downstream fix, not a way to reuse fits after changing
their scientific inputs. The new fingerprint applies to downstream tasks; existing
stage-2 metrics are deliberately reused. Do not regenerate configs or run the normal
model-generating wrapper for this recovery.

`REFINE_MODEL` stages its config as `refinement_input.config.json`, then copies it
to the representation-specific output name. The distinct input name prevents
`cp` from failing with "are the same file" after successful FLARE refinement,
including when refinement is skipped by configuration.

Candidate evaluation requires features for every enabled dataset and truth metrics
for the configured selection role. Missing data fail the task rather than producing
an artificially low selection score. Weighted training errors propagate; estimators
that do not support required sample weights are rejected instead of silently fitting
an unweighted model.

## Large model-strategy run

A practical large run is:

```bash
SEARCH_SPACE=config/search_space_big.json \
MAX_MODELS=0 \
MAX_STAGE2_CANDIDATES=5000 \
./run_in_container_v3.sh
```

The large profile trains about 522 FLARE models before stage-1 screening. The full
stage-2 Cartesian product would exceed 100,000 candidates, so a 5,000–10,000 candidate
cap is much more sensible for the first large run.

## Output

```text
results_v3/
├── panel_probs/
├── models/
├── features/
├── stage1_screen/
├── stage2_candidates/
├── stage3_full/
├── selection/
└── final/
    ├── leaderboard.tsv
    ├── leaderboard.html
    ├── all_metrics.json
    ├── strategy_summary.tsv
    ├── model_family_summary.tsv
    ├── finalists.tsv
    └── finalists/
```

`leaderboard.tsv` includes the upstream strategy and downstream tuning parameters, e.g.
model family, K, covariance, PCA dimension, Bayesian prior, refinement mode, decoder,
leaf size/tree count, calibration and training weights.

`strategy_summary.tsv` compares model-building strategies directly; `model_family_summary.tsv`
aggregates Gaussian-mixture, Bayesian-mixture, K-means and fixed-baseline candidates.

## Scientific interpretation

The purpose of the broad search is not to declare the configuration with the smallest
single development metric the production winner. The automated search should identify a
small, robust shortlist. Those finalists still need the production-domain benchmark:

```text
known high-density truth
    -> microarray mask
    -> SNPster harmonization
    -> Beagle imputation
    -> SNPster GP/DS QC
    -> ancestry model
```

That domain-matched benchmark is the next major dataset factory to add. It is especially
important because Lonni demonstrated that clean-reference latent states can change their
meaning after the real SNPster imputation path.
