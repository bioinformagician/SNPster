# Provisional post-sweep experiment

**Full-resource release verified, 2026-10-25:** the historical releases below
cover an inherited HGDP/1000 Genomes subset. Their replacement is the
[full-resource standardized dataset](README.md#full-resource-standardized-dataset):
one identity ledger, biological-only family ownership, assay-specific panels,
G0-G3 pedigrees and separately unscored affinity/relatedness-stress records.
Use [master_resources_v1_20261025h/full/release](datasets/master_resources_v1_20261025h/full/release)
for the next comparison. All 4,957 resource identities and 441 new descendants
are retained, with 22-autosome/four-panel verification and all five model-input
interfaces passing. All 1,023 old omissions are recovered: 826 in primary roles
and 197 in unscored related-reference stress. The
[observed compact config](datasets/master_resources_v1_20261025h/full/release/panels/observed/configs/compact.json)
is ready, with separate configs for the other assay panels.
Do not move reference relatives into primary test while leaving their families
in fitting, add overlapping panel counts, or treat reused donors as a fresh
external test. Dataset preparation does not authorize model fitting.

**Historical expanded benchmark, verified 2026-10-25:**
[pedigree_observed_v1_20261025_v3/full/release](datasets/pedigree_observed_v1_20261025_v3/full/release)
remains a verified HGDP/1000 Genomes-only development subset, superseded by the
full-resource release above. Its design has training/calibration/validation/test totals of
420/105/210/210, including respectively 224/56/112/112 pure founders followed
by three halving generations. All descendants' phased observed-panel VCFs,
founder-haplotype inheritance and realized ancestry are retained. There are
1,765 separate reference representatives. Both mask input interfaces and all 22
chromosome inheritance audits pass, with no cross-role known/flagged relationship.
The [compact config](datasets/pedigree_observed_v1_20261025_v3/full/release/configs/compact.json)
is ready; **no model sweep is authorized or running as part of this preparation**.
Validation/test each contain 14 terminal pedigree families, not 210 unrelated
people, and OCE/NAM founder diversity remains constrained. Use generation-specific
and family-aware error reports; reused donors still are not a fresh external test.
See the [persistent benchmark contract](README.md#persistent-three-generation-benchmark).

**Dataset rebuild, 2026-10-25:** the approved next step is a reusable,
whole-biological-family/kinship-component cohort, **not another sweep yet**.
[run_kinship_dataset.sh](run_kinship_dataset.sh) exports observed genotypes,
local/whole-genome ancestry, immediate parents and founder-donor metadata,
plus full/compact/smoke input configs. The
[dataset contract](README.md#reusable-kinship-aware-dataset) describes the
schema, relatedness-screen threshold and independent verification gates.
Only a successfully verified new release should feed the next comparison.
That release is now available at
[kinship_observed_v1_20261025_v2/release](datasets/kinship_observed_v1_20261025_v2/release).
Its [compact config](datasets/kinship_observed_v1_20261025_v2/release/configs/compact.json)
uses five base settings and native interpretation only. Both pipeline input
routes have been exercised, but **no new model comparison has run**. The
genotype/pedigree screen found no retained cross-role known/flagged relationship
at its stated boundary. Four/two mixed validation/test people and prior use of
the underlying donors still limit model-selection conclusions.

**Blocker established, 2026-10-24:** the completed clean observed experiment has
biological relatives crossing reference/validation roles despite disjoint
declared donor IDs. Controlled native FLARE experiments identify this as the
trigger for the catastrophic SSA-to-NAM and approximately half-genome KHV
misassignments, through ancestry-specific copying/switch-rate inference.
See the [mechanism and matched exclusions](README.md#established-flare-failure-mechanism-2026-10-24).
Do not finalize this draft using that leaderboard, rerun the unchanged historical
role planner, or attempt to calibrate away these failures. Rebuild a
kinship-aware cohort and repeat a compact, locked native/control comparison
before expanding hyperparameters. The GME source-separability problem remains;
this is not a claim that all FLARE errors were repaired.

**Status update, 2026-10-23:** the user cancelled the 96-setting sweep and requested
a clean standardized launch without model handover. The observed-only replacement
uses the eight provisional base settings and five heads below, but only the two
harmonized masks. Its actual five-role cohort and generated configs are produced
by [run_observed.sh](run_observed.sh); see the
[clean-launch documentation](README.md#clean-observed-only-launch).
It must pass a real-cohort chr22 smoke and independent audit before the full stage.

The [original four-route draft](post_sweep_v1.draft.json) remains a proposal for
extending the experiment to passing imputed inputs, not the active launch config.
The original dataset paths in that draft are not materialized. Do not launch it
unchanged or infer that the interrupted sweep supplied tuned winners. The
replacement's generated config and input verification are retained in its
dedicated run directory.

## Research-led shortlist

| Priority | Comparison | Question |
| --- | --- | --- |
| Primary | Native labelled FLARE | Does a well-chosen local-ancestry model already give the best physical ancestry percentages? |
| Primary | FLARE + constrained simplex map / Ridge | Can a small, interpretable correction remove systematic source confusion without overfitting? |
| Secondary | FLARE + shallow RandomForest / XGBoost | Are there reproducible nonlinear corrections beyond the linear controls? |
| Reference policy | Full versus fine-population-capped FLARE; fixed versus EM generations | Does more balanced representation or target-dependent fitting improve mixed-person transfer? |
| Research comparator | Genuine Gnomix fast LR + CRF | Does its window/smoothing representation outperform FLARE on the same people and features? |
| Controls | Simplex-constrained reference-frequency fit; PCA + Ridge | Is local ancestry adding value over simpler phase-free global inference? |

These are hypotheses, not claims that a paper or synthetic smoke established a
winner. In this implementation, **FLARE + XGBoost means a decoder of the named,
physically aggregated FLARE fractions**. It is not SNP-level XGBoost, a local
posterior-feature model, or an ensemble concatenating multiple base models.
All decoders fit on calibration people who never train the base model.

The draft contains:

- Four provisional FLARE settings: full/cap16 references x fixed/EM, initialized
  at ten generations. `update-p=false` follows the original sweep. These are
  anchors to revise, not four asserted winners.
- Two Gnomix settings sharing the same 0.5-cM LR base but actually fitting
  CRF `c2=0.1` versus `c2=1`. The new pipeline applies this parameter explicitly.
- Frequency and PCA/Ridge controls.
- Native, constrained simplex, Ridge, shallow RF and shallow XGBoost outputs.
  Native labelled FLARE does not need calibration by definition; the learned
  heads must demonstrate added value.

RF is limited to depth 3 and leaf size 10; XGBoost uses depth 2, 100 trees and
learning rate 0.03. These are deliberately small starting hypotheses, not tuned
optima. Check calibration cohort diversity/size before retaining nonlinear heads;
with too few distinct families, remove them rather than fit an impressive-looking
in-sample decoder.

The current workflow applies all five interpretations to all eight bases:
**40 base/interpretation combinations per route, not just eight models**.
Four routes (two stages x two array masks) give 160 route-specific output sets.
This is not another 96-setting base-model grid. At finalization, remove redundant
FLARE settings and do not expand the head/hyperparameter grid repeatedly using
the same validation people.

### Important qualification of the interrupted sweep's Gnomix labels

The interrupted observed expansion passes `crf_c2` to
`model_family_gnomix.fit_chromosome`, but its `initialize_model` does not set the
underlying CRF parameter before training. A read-only check against the pinned
upstream runtime initialized models with requested values 0.01, 0.1, 1 and 10:
all returned effective `smooth.model.CRF.get_params()["c2"] == None`, meaning
the library default, not the requested value.

Consequently **do not select a Gnomix L2 winner from those candidate names**.
Any future analysis of its retained outputs must inspect fitted parameters and
prediction aliases first. No legacy L2 winner was used in the clean launch.
The new pipeline's explicit L2 contrast is a new comparison, not a replication
of an assumed legacy L2 result. The legacy producer was not patched in place;
the run itself was stopped by explicit user request.

### Deferred rather than disguised as supported candidates

- A donor-purged classifier/frequency-feature stack plus boosted decoder is a
  promising phase-free follow-up from the historical experiments. It needs a
  separate feature/OOF adapter; plain frequency + XGBoost is not that model.
- New FLARE2 copying-profile clusters plus a learned source map are interesting,
  but the current adapter does not train that representation. PCA/GMM is not an
  equivalent substitute. Do not add a fictitious FLARE `k=1234` option.
- PCA/RF and PCA/XGBoost global controls are available but are not default
  priorities: richer PCA/tree features have previously failed transfer tests.
- Cross-model stacking and concatenated local-posterior summaries are not
  implemented by this config and must not be inferred from its decoder names.

Primary-source links and historical lessons are in [RESEARCH.md](RESEARCH.md).
FLARE's Apache-2.0 software license motivates its priority; it does not clear
data/dependency licensing. Gnomix remains research-only in separate rankings.

## Dataset and input-stage contract

The draft uses the sweep's **seven-source order**:
EUR, EAS_SEA, SAS, GME, NAM, OCE, SSA. Its 22 physical spans are copied from the
retained observed benchmark cohort, not guessed from contig lengths.
Preserve these boundaries for directly comparable truth, or explicitly freeze
new boundaries and regenerate truth before fitting. An eight-source EAS/SEA
experiment requires authoritative eight-source labels and a separate config;
never manufacture that split from seven-source predictions.

Existing FLARE/Beagle jars, Gnomix checkout and genetic-map paths are real local
assets. The following **planned input bundle does not yet exist**:

```text
inputs/post_sweep_v1/
  samples.tsv
  realized_segments.tsv
  reference/chr{1..22}.vcf.gz
  training/chr{1..22}.vcf.gz
  query/lonni/harmonized/chr{1..22}.vcf.gz
  query/frederik/harmonized/chr{1..22}.vcf.gz
  query/lonni/imputed/chr{1..22}.vcf.gz
  query/frederik/imputed/chr{1..22}.vcf.gz
```

Those are explicit staging destinations in the draft JSON, not working inputs.
Materializing them requires a separately orchestrated preparation step after
the source cohort and role split are agreed. Do not fill them with synthetic
smoke data or silently copy the already inspected 14 controls into a supposedly
independent test set.

Use disjoint reference, training, calibration, validation and test roles with
canonical founder donors and genetically audited family components. Include
pure and mixed people and multiple mixture scenarios. Reserve the external
holdout from model/decoder selection. `independent_test=false` is intentional
until the provenance actually supports a stronger statement.

The reference and labelled training VCFs must retain the eligible **dense**
source panel so that the imputed route can use imputed SNPs. Copying the old
array-only reference/training panels would intersect the imputed route back
down to observed markers and defeat this comparison. Log the retained marker
sets/counts; different stage-native panels do not establish a matched-panel
causal effect of imputation.

Only already-verified production imputation output can populate the imputed
routes. The earlier failed 0.7711-concordance run is **not** a passing source.
Do not relax the existing genotype gates to make this config run. If no passing
imputed cohort is available, explicitly finalize an observed-only config first
and defer imputation; do not relabel phase-only inputs as imputed.

## Finalization after the current sweep

1. Require complete expected chromosome/mask/candidate coverage and verify
   published predictions; do not shortlist from partial task output.
2. Summarize all-autosome, physical-span TV, equal-person mask-averaged loss,
   pure/mixed/scenario behavior, runtime and memory. Avoid selecting a candidate
   that merely improves pure-reference scores while harming mixed people.
3. Use these development results to revise the four FLARE slots: generations,
   EM policy and reference cap. Retain an untuned native anchor plus useful,
   nonredundant alternatives. Deduplicate effective settings; `array=true/false`
   is not an extra model with explicit zero internal MAC/MAF filters.
4. Audit the Gnomix caveat above. Select its window from reliable evidence and
   test explicit CRF L2 values in the new pipeline. Update IDs to match actual
   parameters. Do not claim the old sweep selected L2.
5. Freeze cohort roles/donors, common feature policy, source mapping, exact input
   provenance and the finite decoder budget **before** consulting the new
   validation/test outcomes. The current sweep's inspected controls remain
   development evidence.
6. Materialize/audit the bundle above, or rebind the config to a verified
   equivalent bundle. Verify all 134 cohort files (132 chromosome VCFs plus
   two manifests), phased training truth, span coverage and imputation gates.
7. Save a new `post_sweep_v1.final.json`, freeze its hashes, and run a small
   **real-cohort chr22 smoke** through the same containerized Nextflow workflow.
   The earlier synthetic smoke is not a substitute for this input-level check.
8. Only after that passes, launch the all-autosome final config. Select on
   validation only, then report the locked finalists on the reserved test role.
   Inspect transfer across both input stages and masks; do not automatically
   promote the winner to production.

## Resource plan, not a launch instruction

Once the current sweep has exited and released resources, use up to 22 task
CPUs and a 180-GiB task allocation budget, initially 8 GiB per one-core task,
plus the 2-GiB controller and host headroom. Measure a real dense-panel smoke
first: some reference/PCA/Gnomix tasks may need more memory, in which case reduce
concurrency rather than exceed the budget. No process in this draft is started
automatically when the existing sweep completes.

After approval and successful input smoke, invoke the existing
[launcher](launch.sh) with the **final**, not draft, config:

```bash
SNPSTER_CPUS=22 SNPSTER_MEMORY_GB=180 SNPSTER_TASK_MEMORY_GB=8 \
  bash data_analysis/ancestry_training/launch.sh \
  "$PWD/data_analysis/ancestry_training/post_sweep_v1.final.json" \
  "$PWD/data_analysis/ancestry_training/results_post_sweep_v1"
```

The launcher freezes sources/config, prints the real controller/log command,
saves graphical progress/report/trace evidence and independently verifies
published results. Preserve failed work for resume and clean only the exact
successful cache after verification, as described in [README.md](README.md).
