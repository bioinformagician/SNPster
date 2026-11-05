# Standard ancestry training and evaluation

This is a new, isolated, configuration-driven Nextflow DSL2 workflow. It does not
change production preprocessing,
production ancestry, or PGS. No large experiment is launched by installing it.

The historical [kinship-aware](#reusable-kinship-aware-dataset) and
[three-generation](#persistent-three-generation-benchmark) releases are verified
HGDP/1000 Genomes development subsets, **not a consolidation of every reference
resource**. The full-resource build below replaces that scope assumption.
Use a release for mod
The [original post-sweep draft](post_sweep_v1.draft.json) is not the
active input bundle; its [research-led shortlist and finalization gates](POST_SWEEP_PLAN.md)
remain provisional rather than tuned model winners.

## Full-resource standardized dataset

[master_resources.json](master_resources.json) and
[run_master_resources.sh](run_master_resources.sh) build a reusable dataset,
**not another model sweep**. The local inventory is 4,957 sample identities:
4,091 HGDP/1000 Genomes, 196 GSE77508, 381 GSE73996, 218 Kazakh2025 and
71 locally materialized SGDP records. This is not a claim to have downloaded
every public SGDP person, nor that cross-resource records are distinct people.
The reused SGDP Chane recall control is counted only once.

The inventory retains 3,984 source-labelled and 973 affinity-only records.
Source labels are operational population annotations, not measured individual
purity. Affinity-only people have `null` additive ancestry truth; a Central Asian
population label is not converted into invented seven-source percentages.

The build deliberately separates three concerns:

1. **One complete identity ledger and biological-family split.** Published
   biological families and all-real-person genome-wide KING edges define
   components. Old simulated-pedigree donor bridges do not. Every screened
   pair must have at least 10,000 jointly called, LD-pruned SNPs; an underpowered
   screen fails explicitly instead of declaring people unrelated.
2. **Assay-specific model panels.** The rich observed-mask, GSE77508, GSE73996
   and SGDP panels retain their own matching markers. The old four-way
   43,722-SNP universal intersection is an identity-screen backbone, not the
   ancestry-model backbone. Native projected GT archives preserve missingness
   and source fingerprints; full original VCFs remain at recorded, hashed
   locations. Complete-case model views may remove sites, never fill GT.
3. **Honest evaluation tracks.** All members of a biological component share
   primary role ownership. Relatives in held-out families can contribute
   evaluation data, with family-aware uncertainty. Additional relatives of
   fitting references and high-KING duplicate/MZ candidates remain available
   in explicitly unscored stress tracks, never as independent test people.

The 224/56/112/112 HGDP/1000 Genomes founder design and three halving generations
are preserved. Additional real relatives and source-/platform-stratified
holdouts can increase pure-role counts; the pure counts therefore need not
equal the pedigree seed counts. The same founder trees and meiosis seeds are
replayed per panel, and their realized physical ancestry must agree.
External unphased targets are phase-only processed by role against the fitting
HGDP reference, with unchanged diploid calls and no added variants. Supplied
HGDP founder homologs are preserved. GSE77508's provided phased GTs do not have
independently certified original per-call pre-phase missingness; this limitation
is retained rather than silently advertised as pristine observed truth.

```bash
bash data_analysis/ancestry_training/run_master_resources.sh \
  /absolute/path/to/master_resources.json \
  /absolute/path/to/new_dedicated_run
```

The launcher freezes producers and inputs, uses pinned Nextflow 26.04.4 inside
named Docker controllers, and first runs real all-resource chr22 and larger
chr4 genotype/memory-scale smokes. In the full run, all four chr22 pedigree/export audits must pass before
the remaining model-panel chromosomes start. Defaults are 12 task cores and
96 GiB task memory, plus a 2-GiB controller. Most tasks use one hard-limited
core/8 GiB; genome-wide KING uses four hard-limited cores/32 GiB. Override the
aggregate budget with `SNPSTER_CPUS` and `SNPSTER_MEMORY_GB` (up to 22/178).
Live graphical Docker logs, source/input hashes, DAG, trace, reports and
failure evidence are preserved. Only independently verified successful task
work/cache is eligible for run-scoped cleanup.

For a new split of the same prepared people and source GTs, optionally set
`SNPSTER_PREPARED_RESOURCES=/absolute/path/to/prepared/full/release`.
The import task checks all source/archive hashes, all 22 independent source-GT
audits, compatible complete all-pair evidence and the frozen producer manifest
before reusing **only** preparation. It never imports old role assignments or
descendants and does not assume the original whole workflow completed.

```text
full/release/
  inventory/inventory.json            # Every local resource ID, provenance and label basis
  archive/chr*/projection/            # Source-specific GT/missingness and identity backbone
  kinship/                           # All-pair coverage evidence and related/duplicate flags
  cohorts/                           # Global family ownership and per-panel founder plans
  panels/{observed,gse77508,gse73996,sgdp}/
    input/chr*/base/                  # Complete matched panel and phase-only checks
    genotypes/all/chr*/master/        # Persistent G0-G3 phased VCFs and inherited ancestry
    genotypes/pipeline/chr*/bundle/   # Ready reference/training and unphased query VCFs
    metadata/                        # Samples, segments and parents/ancestry
    configs/{dataset,compact,sweep,smoke}.json
    verification/chromosomes/         # Exact independent founder-haplotype reconstruction
  unscored/chr*/cases/                # Affinity and relatedness/replicate stress VCFs
  master/{dataset,individuals}.json    # Whole-release contract and parent/ancestry ledger
  master/individuals.tsv
  verification/independent.json       # Durable hashes and real pipeline-interface verification
```

Runnable configurations use local absolute genotype/tool paths; descriptors
use release-relative paths. The compact panel configs compare native FLARE,
genuine Gnomix and controls; full configs keep provisional hyperparameters and
optional interpreters. These are inputs, not fitted results. Reused donors are
not a fresh independent test. Fixed production-imputation-reference membership
is explicit, so this observed-stage dataset does not certify independent
imputed-stage accuracy.

### Verified full-resource release

**Verified, 2026-10-25:**
[master_resources_v1_20261025h/full/release](datasets/master_resources_v1_20261025h/full/release).
Start with the [global dataset descriptor](datasets/master_resources_v1_20261025h/full/release/master/dataset.json),
[ancestry/parent ledger](datasets/master_resources_v1_20261025h/full/release/master/individuals.json)
or [tabular identity ledger](datasets/master_resources_v1_20261025h/full/release/master/individuals.tsv).
The ledger covers all 4,957 resource identities and 441 new descendants:
5,398 entries, **not 5,398 independent biological people**.

These are global identity-role counts, not a single assay's sample sizes:

| Primary role | Real source-labelled records | New descendants | Total |
| --- | ---: | ---: | ---: |
| Reference | 2,015 | 0 | 2,015 |
| Train | 779 | 196 | 975 |
| Calibration | 252 | 49 | 301 |
| Validation | 389 | 98 | 487 |
| Test | 331 | 98 | 429 |

The other 1,191 real records are retained separately: 973 affinity-only,
210 related-reference stress and eight high-KING identity-replicate/MZ
candidates. They are not primary percentage-accuracy cases. The source labels
in the table are operational annotations, not measured biological purity.
The 441 descendants comprise 252 G1, 126 G2 and 63 G3 people. Their founder
seeds remain 224/56/112/112 for train/calibration/validation/test; extra real
records are not additional simulated seed founders.

All **1,023 previously omitted HGDP/1000 Genomes people are retained**.
Their new primary assignments are 46 reference, 448 train, 79 calibration,
155 validation and 98 test (826 total); the remaining 197 are retained as
related-reference stress cases. No old simulated donor bridge was treated
as a biological relationship.

Assay panels overlap in their HGDP people and descendants; do not sum their
counts or pretend the global counts describe every model input. For example,
the observed panel has 1,908 reference, 865 train, 185 calibration, 362
validation and 310 test people. Its
[compact comparison](datasets/master_resources_v1_20261025h/full/release/panels/observed/configs/compact.json)
and [full provisional sweep](datasets/master_resources_v1_20261025h/full/release/panels/observed/configs/sweep.json)
are ready. Equivalent configurations and persistent G0-G3 VCFs are available
under each [assay panel](datasets/master_resources_v1_20261025h/full/release/panels).
The configs compare FLARE, genuine Gnomix and controls; no models have been
fitted as part of this release.

All 382 full-stage tasks passed, followed by
[independent published-output verification](datasets/master_resources_v1_20261025h/full/release/verification/independent.json).
The verification covers all 22 autosomes in all four panels, 12,283,446
real-record KING pairs (minimum 32,277 jointly called LD-pruned markers),
1,294 flagged relationships with no cross-role family ownership, and
23,001,219,710 independently reconstructed founder-haplotype alleles.
The same descendants have identical realized physical ancestry across panels.
All five actual pipeline input interfaces passed, including both observed masks.
Scientific files are copied, not links into disposable Nextflow work.

Use calibration only for fitted output interpreters, validation for model
selection and test for the locked comparison. Report generation-specific,
family-aware errors: these larger pools still are not independent individuals,
and the reused donors do not constitute a fresh external test.

### Guided GME and eight-source follow-up

[The bounded follow-up specification](comparisons/guided8_observed_20261030/specification.json)
uses the completed staged genome comparison as its predecessor. Its original
GME reference has **65 people**, spanning Druze, Bedouin, Palestinian and
Mozabite, rather than only the four pure GME people in the small validation
subset. Those four validation people have mean genome-wide native TV errors
of 44.04 pp with full-reference EM, 31.98 pp with fixed generation 2 and
53.02 pp with cap16 EM. This implicates model/reference representation as well
as sample scarcity; it does not isolate a causal mechanism or establish
measured individual biological admixture.

The new screen is deliberately limited to **chr20 and chr22**, both observed
masks, **five FLARE configurations** and **170 development queries**
(43 calibration, 127 validation). It retains the preceding 96 development
people and adds the available pure GME and SEA calibration people, all pure
GME validation people, and all EAS/SEA validation people. There are now
20 pure GME, 40 pure EAS and 12 pure SEA validation people, still belonging
to only 14 validation family components. All 1,908 original reference and
865 training identities remain in the input; test people are excluded.
The cap65 candidate alone deliberately subsets its fitting reference.

The five candidates are:

- The previous winning full-reference EM configuration, with eight sources.
- Full-reference fixed generation 2, the previous challenger.
- EM with 280 rather than 100 IBS states.
- EM allowing `update-p=true` to adapt native copying parameters.
- EM with at most 65 reference people per source, rotating fine populations
  deterministically. All 65 GME reference people and all smaller-source
  reference sets are retained; reference people are never duplicated or
  borrowed from held-out roles.

This is **3,400 singleton FLARE calls**, with a hard maximum of 4,000 and
**no automatic genome-wide expansion**. Only native interpretation is fitted;
three pure SEA calibration people do not justify another learned-decoder grid.
The reference-copying matrix in a native FLARE model is a diagnostic
parameter, not the person's ancestry percentages.

[guided_followup.py](guided_followup.py) restores EAS/SEA using the existing
audited original label authority, not a new geographical guess: CDX and Dai
remain EAS; KHV, Cambodian and Lahu remain SEA. Every descendant interval
is assigned its actual founder donor's original label using the retained
founder-haplotype transmission table. The transmission and authority hashes
must match the independently verified resource release; merging reconstructed
EAS/SEA segments must reproduce the original seven-source physical truth
exactly. Original reference/training identities, genotype files, marker filters,
spans, seed and role ownership are unchanged.

The old merged model itself **cannot** discriminate EAS and SEA. This screen
refits the winning *configuration* with eight sources; it does not relabel old
seven-component weights. Output reporting includes eight-source and matched
old-cohort merged-seven TV, pure-source confusion, separate and class-balanced
EAS/SEA comparisons, full eight-class decisions and EAS/SEA fraction errors
in people with at least 2% realized combined EAS/SEA truth. Reporting separately
reconstructs truth from the original founder transmissions before accepting
the discrimination report. These remain partial-chromosome development scores,
not independent test or genome-wide accuracy.

[guided_followup.nf](guided_followup.nf) reuses the standard fitting,
inference, aggregation and validation processes, with cacheable **16-person
phase chunks** and **32-person singleton FLARE chunks**. The original 96
people's matching singleton phases are reused after exact input/hash/key/
diploid-GT checks; the 74 added people are phased without imputation.
No merged-source FLARE prediction is reused as eight-source evidence.

```bash
SNPSTER_CPUS=22 SNPSTER_MEMORY_GB=160 \
SNPSTER_TASK_MEMORY_GB=6 SNPSTER_PCA_MEMORY_GB=16 \
SNPSTER_JAVA_PROFILE=c1 \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/path/to/comparisons/guided8_observed_20261030/specification.json \
  /absolute/path/to/new_dedicated_run --guided-followup

# After interruption, use the same specification and run path and add --resume.
```

The graphical named controller saves a DAG/report/timeline/log/trace for
planning, the screen and discrimination reporting. Its independent native
posterior/physical-truth audit gates the report and exact successful-stage
cache cleanup. Completed output and audit evidence remain durable; unsuccessful
runs and the predecessor are retained.

### Final closure protocol

[The locked final protocol](comparisons/final_observed_20261030/protocol.json)
limits the final genome-wide comparison to **three finalists on all 22
autosomes**. Keep eight-source truth and predictions unchanged; choose the
winner using mean person-level validation TV after summing EAS+SEA in both
truth and prediction **only during metric calculation**. Average the two
observed masks within each person, and report full eight-source accuracy,
pure/mixed source errors, off-target ancestry, EAS/SEA discrimination and
family-aware uncertainty separately.

There is one last bounded optimization hypothesis: shrink native
reference-estimated copying probabilities modestly toward their *named*
source/reference-panel identity, without changing reference people, genotypes,
phasing, miscopy probabilities or IBD rates. Its chr22 probe uses only the
43 calibration people, both masks and strengths 0, 0.25 and 0.5 (258
singleton calls). The zero-strength serialization control must match the
existing native result within 0.1 TV percentage points. Initial target
ancestry proportions and generations are reset, never borrowed from old
target-fitted outputs. Improvement and regression gates are frozen before
scoring; a failed hypothesis is an explicit valid outcome, not a reason for
another sweep.

The planned finalists are native full-reference EM100, fixed generation2
and one accepted regularized candidate, or EM280 if regularization fails.
Creating this protocol alone does **not** start the final sweep.
Lock one winner on validation and then score only that winner once on the
original untouched test role. Do not tune again, select a different model
from test results or present the repeatedly used validation leaderboard as
independent test accuracy.

The final copying probe rejected both strengths: calibration merged TV
increased from 17.76 pp to 21.94/25.52 pp, and pure GME TV increased from
41.02 pp to 53.33/60.79 pp. The zero-strength control reproduced all
86 existing person/mask predictions exactly. Consequently the final
three candidates are **native EM100, EM280 and fixed generation2**;
neither copying regularization, target copying updates nor balanced
reference capping is advanced.

The optional config boolean `combine_eas_sea_for_metrics` activates merged
scoring; it is false by default. It requires both separate source axes,
does not alter any predictions/truth and preserves separate-source TV
in prediction evidence alongside the selected metric.

[launch_final.sh](launch_final.sh) executes the final
[Nextflow metadata stages](final_comparison.nf) and the standard sharded
scientific workflow. Audited chr20/22 native predictions and phases are
reused, together with the original 96-person phases on the other autosomes.
The genome stage requires 20,400 new calls rather than 22,440.
After its audit, a durable winner lock gates a **single-model** test of the
310 original test people. Development predictions/phases are reused in
that stage; only the test people require new native calls (13,640).
The final report includes both metrics, source precision, EAS/SEA
discrimination and explicit limits. No production promotion is automatic.
Exact successful work caches are cleaned only after both audits and the
final report have passed.

```bash
bash data_analysis/ancestry_training/launch_final.sh \
  /absolute/path/to/comparisons/final_observed_20261030/protocol.json \
  /absolute/path/to/new_dedicated_final_run

# To resume an interrupted run, use the same paths and add --resume.
```

### Prepared fine-population ablation

The subsequently requested [population experiment](comparisons/population_observed_20261031/specification.json)
is **prepared, not automatically launched**. It is separate development work,
not an extension of the final locked test. Its single candidate changes only
the native FLARE reference-panel labels, matching the predecessor's
validation-locked winner's settings. It uses the original **68 dataset-qualified
reference populations** among all 1,908 references, not countries or
nationalities. Fine-population probabilities are retained and **summed**, not
hard-assigned, through an original-reference-only one-hot map into the same
eight sources. No calibration or evaluation labels fit this map.

The scope is chr20+22, both observed masks and the same 170 development people:
**680 new candidate singleton calls**, with 680 matched baseline calls and
all phases reused after hash, command, marker and genotype checks. Training,
truth, reference identities and source taxonomy remain unchanged. Preserve
EAS/SEA separately in outputs; combine them only for the primary TV metric.
Also report separate-source TV, pure/mixed source errors, absent-source calls,
EAS/SEA discrimination, per-mask results and paired family uncertainty.
There is no automatic genome-wide expansion, promotion or test stage.

Reference scarcity is explicitly reported: 19 fine populations have fewer
than ten reference people, including 1KG_CHS (one), 1KG_IBS (two) and
HGDP_Melanesian (two). This may make finer modeling less stable, not better.
There is no known country/fine-population truth for the descendants and no
claim to measure accuracy at that level.

[population_followup.py](population_followup.py) can prepare the reference
counts/map without inference. [population_followup.nf](population_followup.nf)
provides cacheable metadata/report stages; [launch_population.sh](launch_population.sh)
uses the standard sharded scientific workflow. It refuses to create/start a
run until the predecessor has completed, locked its winner and passed cleanup.
It reads only the predecessor's **development** evidence and validation lock,
never its test predictions or scores. Frozen sources/runtime, named Docker
controllers, graphical progress, execution reports, independent native-posterior
auditing and exact successful-stage cache cleanup are retained.

```bash
# Explicit future launch only, after reviewing the completed predecessor.
bash data_analysis/ancestry_training/launch_population.sh \
  /absolute/path/to/comparisons/population_observed_20261031/specification.json \
  /absolute/path/to/new_dedicated_population_run

# After interruption, use the same paths and add --resume.
```

The predecessor's locked test remains evidence for its own model only.
This additional, repeatedly used development cohort cannot give the new
candidate fresh held-out-test validation; that would require untouched data.

#### Parallel chr22 pilot

The user-authorized [parallel pilot](comparisons/population_chr22_parallel_20261031/specification.json)
is a separate version2 protocol: **chr22 only**, two masks, the same170
development people and **one fine-population EM100 candidate**. Its comparator
is the already-audited guided `full_em100`, not the unfinished final sweep's
winner. Exactly **340 new** fine-panel calls and **340 reused** broad-panel
calls are planned; phases are reused with the same strict checks.
The original post-winner chr20/22 protocol remains gated and unchanged.

Give each orchestration **11 task CPUs and80GiB**, using the same pinned
runtime and six-GiB single-CPU tasks. Nextflow controllers themselves remain
lightweight; allocating11 controller CPUs would not limit their task workers.
The launcher also sets **`-pool-size 11`** independently of the controller's
CPU quota. In the pinned Nextflow version, the default local launch pool uses
the controller's visible CPU count and waits for each child task to finish.
Without this explicit setting,11 scheduler slots can contain only two running
workers and nine submitted-but-not-started tasks. Executor capacity alone is
not evidence of actual concurrency; verify running worker containers and CPU
usage as well.
Population inference defaults to **one person per cacheable task**, so the
remaining people can use all available scheduler slots instead of being trapped
in a few serial 32-person chunks. `SNPSTER_INFERENCE_CHUNK_SIZE` accepts a positive
integer override. It is an execution setting passed explicitly on launch/resume,
not a change to the frozen scientific source, reference, seed or EM settings.
Rescheduling interrupted chunks requires a scoped stop and verified singleton
recovery first; changing chunk boundaries alone does not preserve partial work.
The existing final run must be deliberately paused/resumed with the smaller
budget without changing its frozen source, inputs, model settings or protocol.
Completed task caches are retained; interrupted chunks may restart.

[inference_recovery.py](inference_recovery.py) can preserve completed native
singletons from interrupted chunks after an exact scoped stop receipt.
It verifies frozen inputs, native commands, reference/parameter signatures,
phased genotypes, marker axes and posterior completeness before making durable
snapshots. Resume with `SNPSTER_NATIVE_RECOVERY=/absolute/path/to/native-reuse.json`
to reuse them through the standard native-reuse path without editing frozen
scientific sources or the original planning manifest. Unfinished singleton
calls are explicitly excluded and rerun.
The launcher persists the recovery path and hash for later resumes; changing
that manifest is rejected.

```bash
SNPSTER_CPUS=11 SNPSTER_MEMORY_GB=80 \
  bash data_analysis/ancestry_training/launch_population.sh \
  /absolute/path/to/comparisons/population_chr22_parallel_20261031/specification.json \
  /absolute/path/to/new_dedicated_chr22_population_run
```

This pilot never reads or influences the ongoing final winner lock or test.
No subsequent model sweep or genome-wide expansion starts automatically.
Chr22 metrics are a diagnostic preview, not genome-wide accuracy.

#### How the ancestry models are connected

The current three finalists are **standalone native FLARE models**, not
PCA/GMM clusters fed into FLARE. Common reference-qualified markers and
exact alleles are aligned first. Beagle phases each observed query without
imputation or changes to diploid genotypes; retained singleton phases can be
reused. FLARE receives those phased query haplotypes, phased reference
haplotypes, their known population labels and a genetic recombination map.

FLARE's probabilistic haplotype-copying model estimates local ancestry along
each chromosome. Similar reference haplotypes support copying states, and
the model accounts for ancestry transitions and imperfect copying. Native
reference parameters and target-dependent ancestry/generation estimates are
not fitted from evaluation truth. `states=100/280` changes the internal IBS
state budget, not the number of ancestry classes or clusters of people.
EM100 and EM280 allow native target EM; fixed2 disables EM and fixes the
generation parameter at2. All three use `update-p=false`.

The native posterior probabilities for both haplotypes are integrated over
physical chromosome length, then combined across chromosomes by physical
length. They are not hard-called before aggregation. The fine-population
candidate replaces the eight reference labels with68 original labels and
sums their soft probabilities into the eight sources using metadata alone.

Earlier comparison families are alternatives, not mandatory preceding stages:

- **Frequency/NNLS:** estimate source allele frequencies from references,
  then fit nonnegative ancestry fractions summing to one to query dosage.
- **PCA + ridge/random forest/XGBoost:** fit PCA on references, then train a
  supervised proportion-prediction head on training-only ancestry truth.
- **PCA + GMM:** fit unsupervised clusters on reference PCA coordinates and
  map cluster memberships to known sources using reference labels. This is
  the family that clusters people; it does not supply FLARE's input.
- **Gnomix:** supervised haplotype-window classifiers fitted on labelled
  references, followed by a CRF smoother fitted on training-only local truth.

Optional learned interpreters in earlier screens mapped raw model components
to ancestry using calibration-only truth. The current finalists use native
interpretation, and the fine model uses only a deterministic reference-label
sum: neither adds a learned calibration decoder. Validation ranks settings;
the locked final winner alone gets its prescribed test evaluation.

### Focused model grid and staged development

The original full-cohort chr22 smoke is superseded for new development work:
it is a chromosome benchmark, not a cheap runtime check. The reduced staged
comparison keeps **all 1,908 reference and 865 training people**, but executes
only **32 calibration and 64 validation queries**. Metadata-only seeded selection
balances pure sources, fine populations, biological families and G1-G3 tracks.
All seven calibration and 14 validation family components remain represented;
they are not 96 independent people. The test role is absent from the execution
manifest, and no test comparison starts automatically.

The stages are:

1. Screen the 16 settings below on chr22 and both observed masks, with all five
   interpretations. Reuse verified completed singleton phase and FLARE evidence
   from the explicitly stopped predecessor.
2. Independently reconstruct and audit the published screen.
3. Advance **five distinct FLARE configurations**, ranked by their best
   validation-only interpretation and mean-person TV. Record the other five
   configurations as reserves; additions require a documented scientific reason
   and a revised cost decision, not automatic CI-overlap expansion.
4. Compare those five and all six inexpensive Gnomix/control configurations on
   all 22 autosomes and the same 96 development queries, then independently audit.
   The second stage is development validation, **not a held-out test result**.
   Its mean TV is not directly comparable with the old 857-person benchmark:
   compare matched people and source/generation breakdowns instead.

This reduces the planned genome-wide FLARE calls from **377,080 to 21,120**.
The automatic development gate refuses more than 25,000 calls. Cheap Gnomix
inference does not make phasing free: both masks still receive singleton
phase-only Beagle preprocessing, shared across phase-dependent models.
Frequency fitting and PCA/Ridge are phase-free, but do not eliminate phasing
when Gnomix or FLARE is also in the same comparison.

The frozen release has seven operational sources, including **EAS_SEA**.
Neither the screen nor the genome-wide development result certifies separate
EAS and SEA estimates. Recovering an eight-source comparison requires
authoritative founder labels and recomputed inherited truth, not splitting a
merged prediction after the fact. Source-resolution claims must follow the
actual frozen truth axis.

[staged_comparison.py](staged_comparison.py) implements the cohort/cost/selection
gates. [staged.nf](staged.nf) makes planning and advancement cacheable Nextflow
tasks; model execution continues through [main.nf](main.nf).
[projection.py](projection.py) checks immutable producers, input hashes,
unchanged full reference/training axes, exact query identities, physical truth,
marker keys, genotype-preserving singleton phasing, model flags and C1 runtime.
It refits bases/heads against the new panel instead of rewriting fitted-object
hashes. Each reused FLARE target is checked against its exact phased query and
its raw posterior is reintegrated. Missing declared completed evidence fails
explicitly; incomplete singletons are enumerated as fresh work.

```bash
SNPSTER_CPUS=22 SNPSTER_MEMORY_GB=160 \
SNPSTER_TASK_MEMORY_GB=6 SNPSTER_PCA_MEMORY_GB=16 \
SNPSTER_JAVA_PROFILE=c1 \
SNPSTER_STAGE_SOURCE_RUN=/absolute/path/to/stopped_superseded_run \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/path/to/focused16/full.json \
  /absolute/path/to/new_dedicated_run --staged-development

# Resume only this exact frozen scope/source run: add --resume.
```

The named, TTY-enabled Docker controller logs the graphical Nextflow pipeline,
stage transitions, finalist list, audits and exact successful-cache cleanup.
Each stage also saves a DAG, report, timeline, plain log and trace. After both
independent scientific audits pass, only this new run's four named stage work
directories are removed, followed by durable-hash reverification. Failed and
superseded runs, original sources and shared scientific inputs are untouched.
Completed, cleaned staged runs refuse `--resume` rather than silently
recomputing the genome-wide comparison.

#### Original grid and legacy smoke/full option

The compact five-setting configuration is a diagnostic baseline, not a
comprehensive sweep. The approved
[16-setting observed comparison](comparisons/focused16_observed_20261025/full.json)
expands the verified release without changing its people, splits, marker filters
or seeds:

- **Ten FLARE settings:** fixed generations 2/5/10/20 with full versus
  fine-population-cap16 references, plus EM initialized at ten generations
  under both reference policies.
- **Four genuine Gnomix settings:** 0.25/0.5-cM windows crossed with CRF L2
  strengths 0.1/1. Window base fits are shared across smoothing strengths.
- **Two controls:** reference-frequency fitting and PCA/Ridge.
- **Five interpretations per setting:** native, simplex, Ridge, RandomForest
  and XGBoost. These are 80 model/interpretation combinations per mask, not
  80 independently trained SNP-level models.

Both lonni and frederik observed masks are included. Other assay panels remain
available for separate comparisons; this run does not silently pool their
different marker sets. Gnomix remains a research-only comparator.

[prepare_comparison.py](prepare_comparison.py) materializes this grid from the
reviewed full config and records its source hash. The optional `--smoke-first`
launcher mode exercises the **same models, interpretations and role cohorts**
on chr22, independently audits those outputs, and only then starts all 22
autosomes. Failure at either gate retains work/history and prevents promotion.

```bash
SNPSTER_CPUS=22 SNPSTER_MEMORY_GB=160 \
SNPSTER_TASK_MEMORY_GB=6 SNPSTER_PCA_MEMORY_GB=16 \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/path/to/focused16/full.json \
  /absolute/path/to/new_dedicated_run --smoke-first

# Same frozen inputs and stage histories:
# Add --resume, retaining --smoke-first.
```

For a native JVM abort in C2-compiled code, the explicit
`SNPSTER_JAVA_PROFILE=c1` recovery profile disables C2 in Java task
containers using `JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1`. It does not
change genotypes, models, seeds, scientific filters, installed software or frozen
producers, and is not an automatic retry or proof of an underlying JVM bug.
The profile selects a distinct pinned image derived from the original task image
by adding only that environment variable; see
[java_c1.Dockerfile](java_c1.Dockerfile). The launcher snapshots
[java_c1.config](java_c1.config) outside the frozen source, records its hash,
and remembers the selected profile for later resumes.
Switching back to default JIT in that same run is refused.

New launches split FLARE inference into cacheable **single-person chunks**, scheduled
by Nextflow within the existing CPU/memory budget. Every person still receives
their own FLARE invocation with the same reference, parameters and seed: this
does not pool people into a multi-target EM fit or change model semantics.
Non-FLARE inference remains vectorized. A separate collection task requires
complete, nonoverlapping person coverage, identical model/input hashes and
unchanged native evidence, then restores the original person order and output
layout before chromosome aggregation and the independent audit.
Collection starts as soon as that model/chromosome/mask's known chunk count is
complete; it does not wait for the entire inference sweep.

This avoids the old long tail in which only ten one-core tasks remained, each
looping serially over all 857 query people. It also bounds the work lost on a
failed task to one chunk instead of a whole model/chromosome/mask. It does not
remove the cost of the individual FLARE fits or guarantee continuous full CPU
utilization. Existing frozen runs and their source snapshots are not modified;
resuming an old snapshot retains its original scheduling unless its launcher
explicitly supplies an authorized chunk-size override, as the population launcher
does. Preserve verified completed singletons before rescheduling partial chunks.

An explicitly authorized migration can retain work from a stopped snapshot.
[reuse.py](reuse.py) checks unchanged scientific producers, the exact scope,
singleton flags and the pinned JVM profile, then freezes hashes for completed
preparation/phasing/base/head outputs and finished individual FLARE exports.
Set `SNPSTER_REUSE_MANIFEST=/absolute/path/to/verified-manifest.json` for the
**new** launch. The launcher freezes that manifest and its hash; a resume cannot
add or alter it. It is used only for smoke/single scope, never silently applied
to the full-autosome scope with different metadata.

Adoption copies and rechecks complete artifacts, explicitly recording
`adopted_not_refitted`. Reused FLARE singletons must match the exact phased GT
and allele keys and pass the same native-posterior checks as fresh inference.
Unfinished individuals are recomputed; changed completed evidence fails rather
than silently falling back. The final independent audit includes all copied
native evidence and does not depend on retaining the old work cache.

The changed image identity invalidates old `PHASE` and native-image `INFER`
tasks, including directory-input tasks. Gnomix keeps its separate pinned image.
Other unchanged tasks retain their caches; subsequent C1 resumes reuse completed
C1 tasks. Do not rely on a changed cache hashing mode alone to invalidate
directory-input tasks. A pinned-container directory-input cache-transition check
verified both invalidation and reuse. Exact-target checks confirmed identical genotype/posterior records
for phase-only Beagle, fixed FLARE and EM FLARE under default versus C1-only
Java. This avoids the observed C2 crash path, but does not establish whether the
intermittent native fault originated in the JVM or the host. C1 may run slower.

```bash
# Opt in once when resuming an affected run; keep the original config and path.
SNPSTER_JAVA_PROFILE=c1 SNPSTER_CPUS=22 SNPSTER_MEMORY_GB=160 \
SNPSTER_TASK_MEMORY_GB=6 SNPSTER_PCA_MEMORY_GB=16 \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/run/experiment.json /absolute/run --smoke-first --resume
```

Every task has a one-core Docker cap and a hard memory cap. PCA base tasks
receive 16 GiB in this profile; other tasks receive 6 GiB. Nextflow schedules
up to 22 task cores within the 160-GiB task budget, plus a separate 2-GiB
controller, leaving host headroom rather than promising 22 simultaneous jobs
regardless of their memory requirements.

## Reusable kinship-aware dataset

[run_kinship_dataset.sh](run_kinship_dataset.sh) prepares data only. It does not
train any models or start a sweep. It preserves the old experiments and builds
a new versioned release from their certified observed SNP genotypes.

**Verified release, 2026-10-25:**
[kinship_observed_v1_20261025_v2/release](datasets/kinship_observed_v1_20261025_v2/release),
with a [portable dataset descriptor](datasets/kinship_observed_v1_20261025_v2/release/configs/dataset.json),
[ancestry/parent metadata](datasets/kinship_observed_v1_20261025_v2/release/metadata/split/individuals.json)
and [independent verification](datasets/kinship_observed_v1_20261025_v2/release/verification/independent.json).

| Pipeline role | Pure | Mixed | Total |
| --- | ---: | ---: | ---: |
| Reference | 2,057 | 0 | 2,057 |
| Train | 362 | 71 | 433 |
| Calibration | 29 | 6 | 35 |
| Validation | 23 | 4 | 27 |
| Test | 28 | 2 | 30 |

All 70 preparation tasks passed. A separately frozen publication-only recovery
corrected an extra `configs/` directory without changing genotypes or roles.
Independent checks verified all 22 genotype chromosomes, 3,514 ancestry records
and 1,191 flagged genome-wide relative pairs on 151,709 LD-pruned SNPs, with
zero cross-role known/flagged relationships. Both actual pipeline input adapters
passed again after successful work/cache cleanup. The release retains 3,375
candidate/donor genotypes and metadata for 3,630 people/ancestor nodes; 793
reference relatives were excluded from the pipeline panel but not deleted from
the archive. All 36 regression tests pass. No models were fitted.

The small mixed validation/test counts (four/two) remain an important limitation
for any later accuracy claim. The completed run's removed task caches cannot
provide cached resume, but published data, hashes, source snapshots, reports and
failure/recovery evidence remain. The
[compact five-base/native-only config](datasets/kinship_observed_v1_20261025_v2/release/configs/compact.json)
and [full provisional comparison config](datasets/kinship_observed_v1_20261025_v2/release/configs/sweep.json)
are ready as inputs, not recorded model results.

```bash
bash data_analysis/ancestry_training/run_kinship_dataset.sh \
  /absolute/path/to/completed_observed_experiment \
  /absolute/path/to/new_dataset_run
```

The named, pinned Docker controller runs [prepare_kinship.nf](prepare_kinship.nf):
pedigree catalog, real chr22 genotype/PLINK smoke, remaining autosomes,
genome-wide KING screening, whole-component role assignment, pipeline exports
and independent genotype audits. A published-output audit then checks ancestry,
family/kinship constraints and actual `freeze`/`prepare` interfaces for both masks.
Default limits are twenty hard-one-core task slots, 176 GiB task memory and a
2-GiB controller. Source/input hashes, graph, traces, HTML report/timeline, logs
and exit statuses are preserved. Use `--resume` only for the identical frozen
source, input files and dedicated failed/incomplete run.

The release layout is:

```text
release/
  configs/dataset.json                 # Portable descriptor, schemas, paths, hashes, limitations
  configs/{sweep,compact,smoke}.json    # Local model-pipeline input configs; not executed
  metadata/split/samples.tsv           # Included samples and five pipeline roles
  metadata/split/individuals.json      # Candidate, ancestor and known biological-parent metadata
  metadata/split/individuals.tsv       # Flat per-person ancestry/parent table
  metadata/split/{split.json,exclusions.tsv}
  catalog/catalog/segments.tsv         # Physical local ancestry truth, both homologs/all autosomes
  genotypes/all/chr*/master/            # All candidate/donor phased SNPs and original mask sites
  genotypes/pipeline/chr*/bundle/       # Reference/train phased VCFs; two unphased query masks
  kinship/screen/                      # LD-pruned all-autosome data and flagged related pairs
  verification/                       # Independent chromosome and final bundle audits
```

The descriptor's paths are relative to `release/`, allowing other pipelines to
consume the data without adopting this workflow. The runnable configs have
absolute local genotype/tool/map paths; regenerate or explicitly rebind these
if relocating the release. Data are GRCh38, biallelic A/C/G/T SNP GTs, without
imputation, genotype filling or dosage substitution. Reference and training
retain the union of the original masks; each harmonized query retains its
original marker set and unchanged diploid genotype. The all-candidate VCFs also
retain excluded reference relatives so that exclusions are reversible and
future splits remain possible without reconstructing source data.

### Ancestry and pedigree metadata

`individuals.json` is keyed by sample identity and includes real people,
simulated descendants, founder aliases and their necessary ancestors:

- Immediate `father_id` / `mother_id`, distinct from `donor_id` for a simulated
  founder that copies a real donor. Available biological parentage comes from
  the hash-pinned public 1000 Genomes pedigree; simulated parentage comes from
  the original plans, not ancestry predictions.
- `expected_ancestry`: Mendelian pedigree expectation for simulations, or the
  operational source-population label for real source donors.
- `realized_ancestry`: ancestry integrated over certified physical segments,
  both homologs and all 22 autosomes. This is not assumed equal to the
  parental average in mixed descendants.
- `parent_ancestry`: parent IDs and each parent's expected/realized vector,
  including intermediate simulated parents not retained as genotype candidates.
- Genotype availability, inclusion/exclusion, pipeline role, component role,
  original fine population and stable family/kinship component ID.

Fractions use 0-1 units in the declared seven-source order; multiply by 100 for
percentages. `null` means unknown/unavailable, never zero ancestry. Parentage is
not invented for real people outside the available pedigrees. A one-hot
source-population label does not certify that an actual person's biological
ancestry is exactly pure.

For example, a consumer can read ancestry and immediate-parent ancestry without
loading any genotype files:

```python
import json
from pathlib import Path

release = Path("/absolute/path/to/dataset_run/release")
descriptor = json.loads((release / "configs/dataset.json").read_text())
metadata = json.loads((release / descriptor["individuals"]).read_text())
person = metadata["individuals"]["HG02954"]
print(person["realized_ancestry"])
print(person["father_id"], person["mother_id"], person["parent_ancestry"])
```

The TSV contains the same per-person fractions and parental IDs in flat form;
unavailable JSON values are blank TSV cells. Preserve the run's `inputs/` beside
`release/` to repeat independent verification against the original public
pedigree, and preserve `source/` and `reports/` for execution provenance.

### Split and screening contract

The graph combines known biological families, every retained simulated
parent/founder relationship, and PLINK2 KING-robust pairs at coefficient >=0.0884
(a second-degree screening boundary). The all-pairs screen uses all 22 autosomes,
complete GTs, pooled MAF >=0.05 and per-chromosome LD pruning at 200 kb / r2 0.2;
at least 10,000 retained SNPs are required. Previously absent real founder donors
are extracted from the hash-certified dense resources before screening.

Whole graph components receive one role, without reading model predictions.
The largest mixed-family component trains the models; other separate mixed
components supply calibration, validation and test. Additional pure components
are chosen deterministically to supply at least three people per source per
query role. Remaining original training components stay in training; remaining
components form the reference, retaining one representative per reference
component. Too few separate mixed families or reference-source coverage is an
explicit failure, not a fallback to a leaking split.

This is still previously studied development data. Family/second-degree
screening does not make it a fresh independent test, prove all distant
relationships absent, or eliminate ancestry/inbreeding sensitivity of KING.
`independent_test=false`, `scientific_certification=false` and
`production_promotion=false` remain explicit.

## Clean observed-only launch

The old 96-setting sweep was stopped on 2026-10-23 at the user's request. The
replacement does not adopt its model fits. [run_observed.sh](run_observed.sh)
starts one named, pinned Docker controller and runs three containerized Nextflow
workflows in sequence: input preparation, a real-cohort chr22 smoke, and the full
all-autosome experiment. The smoke's independent durable-output audit must pass
before the full experiment starts. Failures stop the sequence and retain work,
reports, source snapshots and exact-run history.
The chr22 smoke retains the real reference and training people and exercises
all model/head paths with one pure and one mixed person per query role (six
query people total). It is a mechanics check, not an accuracy leaderboard.
The subsequent full experiment uses all 65 calibration/validation/test people.

```bash
bash data_analysis/ancestry_training/run_observed.sh \
  /absolute/path/to/results_v3_observed_dense_01 \
  /absolute/path/to/new_dedicated_experiment
# To resume the same frozen source/config/code and all three stage histories:
bash data_analysis/ancestry_training/run_observed.sh \
  /absolute/path/to/results_v3_observed_dense_01 \
  /absolute/path/to/the_same_experiment --resume
```

[prepare_observed.py](prepare_observed.py) and
[prepare_observed.nf](prepare_observed.nf) materialize the draft's missing inputs
from existing certified observed GTs and physically realized ancestry segments.
The two masks are unioned for the common reference/training inputs, with an
explicit equality check for phased GT on overlapping exact keys. Each query mask
retains its own original markers. Removing phase does not change diploid GT;
the standard workflow then performs its independently checked phase-only adapter.
No imputed route or failed imputation epoch is relabelled or admitted.

### Native FLARE mechanism diagnostics

[run_flare_diagnostic.sh](run_flare_diagnostic.sh) launches a separate, frozen,
containerized [Nextflow investigation](diagnose_flare.nf) against the clean run's
published inputs, without modifying its models or using final-test people:

```bash
bash data_analysis/ancestry_training/run_flare_diagnostic.sh \
  /absolute/path/to/completed_observed_experiment \
  /absolute/path/to/new_flare_diagnostic
```

The [locked protocol](flare_diagnostic_v1.json) examines three validation outliers
and two hash-selected pure training controls on chromosome 22 and both masks.
It verifies exact genotypes, authoritative source/fine-population labels,
reference exclusions and historical input certificates before inference.
An exact native-baseline replay must pass in a real task container before the
remaining calls start. Matched contrasts cover stock versus precision export,
inherited versus adapter phase, reference-order permutation, IBS/state settings,
seed and native EM updates. Training controls retain their inherited phase;
they are not claimed to reproduce a singleton query-phasing experiment.
Previously completed theta-crossing and diagonal-copying controls are prior
negative evidence, not presumed fixes.

The controller reserves 2 GiB and by default schedules fourteen 1-core/4-GiB tasks
(56 GiB total task memory), leaving room for the original run's independent
audit. Reports, input/source hashes, native models, local/global posteriors,
commands and diagnostic comparisons are retained. Final verification separately
reconstructs every local-posterior physical integral. These are mechanistic
development tests, not independent accuracy certification or model promotion.

The launcher snapshots itself before waiting for the controller, so subsequent
workspace edits cannot corrupt the host-side log/exit-status archival step.
Use `SNPSTER_DIAGNOSTIC_PROTOCOL` to select the
[parameter/relative-exclusion follow-up](flare_diagnostic_v2.json) or the locked
[chr1](flare_diagnostic_chr1.json) and [chr12](flare_diagnostic_chr12.json)
replications or the [matched sham exclusion](flare_diagnostic_sham.json).
`SNPSTER_DIAGNOSTIC_CPUS` overrides the default fourteen task slots;
each slot still reserves 4 GiB and has a hard one-core quota. The follow-up
crosses copying matrices and ancestry-specific haplotype switch rates by
**source name**, not array index. It compares native re-estimation with unchanged
native parameters on the same filtered references.

### Established FLARE failure mechanism, 2026-10-24

Five dedicated Nextflow runs completed 368 matched calls plus five real-container
smoke replays (373 native FLARE invocations) on chr1, chr12 and chr22. Every
matched call passed independent native-posterior reconstruction.
Source/fine labels, exact target and
reference GTs and stock-versus-precision inference were consistent. Export
rounding, source-axis display, reference order, seeds, phase inheritance and
larger IBS/state settings did not explain the catastrophic SSA-to-NAM call.

There is **real biological relatedness across the declared roles**:

- Validation person HG02954 has both parents, HG02952 and HG02953, in the
  reference, including cap16.
- Validation person HG02028 has her child HG02027 in the full reference, but
  not cap16.

These relationships are explicitly present in the public
[1000 Genomes pedigree](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped)
and independently supported by genotype-only chr22 coefficients near 0.25 with
zero or one opposite-homozygote SNP. The previously checked declared donor
components did not capture these biological families.

| chr22 case/reference | Mask | Native fixed TV (pp) | Exclude flagged relatives, fixed | Exclude flagged relatives, EM |
| --- | --- | ---: | ---: | ---: |
| HG02954 / cap16 | Lonni | 100.000 | 1.481 | 0.000 |
| HG02954 / cap16 | Frederik | 100.000 | 0.358 | 0.000 |
| HG02028 / full | Lonni | 50.344 | 2.045 | 0.000 |
| HG02028 / full | Frederik | 50.533 | 0.716 | 0.000 |

The same exclusions with the original copying/rho parameters locked also
removed the catastrophic calls. Removing only one of HG02954's two parents
left approximately 50% error; removing both removed the catastrophic component.
These tests do not change the target GT or phase and do not inject target
ancestry into inference.
The effect replicated on chr1 and chr12, on both masks: all eight
locked-parameter exclusions reduced the affected cases to less than 2.04 pp
error. These are three-chromosome mechanism replications, not a corrected
all-autosome accuracy certificate.

The sham exclusion removed the same numbers of unrelated people from the
**same fine populations**, leaving the parents/child in the reference.
HG02954 still had 100 pp chr22 error on both masks, and HG02028 still had
approximately 50 pp error. This also held with the native copying/rho parameters
locked. The rescue is therefore not simply a smaller reference, altered source
counts, different fine-population representation or parameter re-estimation.

The mechanism is FLARE's ancestry-specific haplotype-switch model interacting
with exceptionally long close-relative matches. On capped Lonni chr22,
SSA's fitted rho is 118.879 and NAM's is 72.738; the copying model permits
cross-panel copying. A long inherited SSA haplotype can therefore be attributed
to the lower-switch-rate NAM ancestry. Transplanting only rho from cap16 into
the full Lonni model changed HG02954's error from 41.36 to 100 pp; transplanting
only the copying matrix did not. Equalizing rho greatly reduced the problem,
but is a diagnostic counterfactual, **not an approved universal fix**.
Target EM then reinforces the incorrect ancestry; it is not necessary for the
failure, which occurs with uniform priors and fixed generation 10.

Do not call this proof that every FLARE error is solved. HGDP00691's GME
confusion remains substantial after excluding a genotype-flagged reference
relative; indeed that relative had helped its full-reference prediction.
Broad-source separability, reference composition and operational one-hot
population labels still require evaluation after correcting the cohort.

**Do not finalize a new hyperparameter sweep on the existing role bundle.**
The historical observed-cohort planner is not a biological-kinship-aware
planner. Future cohorts must connect known biological families as well as
declared founder donors, screen reference/query relatedness independently, and
reassign/exclude affected reference people before any model fits or calibration.
The per-chromosome diagnostic threshold is not a genome-wide kinship
certification. The frozen completed run remains an execution-verified
development experiment, not a leakage-free model-selection certificate.
Genotype/axis and named-parameter-transfer regressions are in
[test_flare_diagnostic.py](test_flare_diagnostic.py).

### Clean observed cohort and resources

Selection is deterministic and does not read model performance. Shared founder
donors connect the 25 retained training families into three components (55, 5 and
4 donors); the largest remains training, and the other two become validation and
test. Their pure founders follow the same role. The 14 earlier observed query
people become calibration only. Five additional pure people per source are
reserved from the reference (three calibration, one validation, one test) before
new fitting. Connected founder IDs, sample identity and family ownership are
checked across all five roles. The resulting counts are 2,998 reference, 300
training, 35 calibration, 16 validation and 14 test people.
These declared-component checks did not exclude the biological relatives
identified above; do not interpret them as genetic family disjointness.

This is **previously studied development data**, not a fresh independent test or
a genetic-kinship certification. Validation/test each contain only one small
founder-connected pedigree group; pure/mixed results and these sampling limits
must remain explicit. Full `independent_test` remains false.

The initial experiment retains the draft's eight base configurations (four
FLARE, two genuinely different Gnomix CRF penalties, frequency and PCA/Ridge)
and five interpretation heads, but only the two harmonized routes: **80 output
sets**, not the draft's four-stage 160. Hyperparameters are provisional anchors,
not winners inferred from the interrupted sweep. The generated final configs,
role manifests, certificates and hashes are published below `prepare/results`.

Default scheduling permits up to **22 one-core tasks** and **178 GiB task memory
plus a 2-GiB controller** (180 GiB combined). Normal tasks reserve 8 GiB; PCA
bases reserve 16 GiB, so Nextflow reduces concurrency when necessary. Override
`SNPSTER_CPUS`, `SNPSTER_MEMORY_GB` (task budget) and
`SNPSTER_TASK_MEMORY_GB` explicitly when needed. Scientific VCF loading uses
compact int8 genotype blocks instead of per-genotype Python objects; physical
truth is indexed per chromosome/sample/homolog without repeated full-table scans.
Task Docker containers also have a hard one-CPU quota, not merely CPU shares.
The 70%-of-normal-task-memory Java heap is passed explicitly into each container;
an executor-side shell export alone does not propagate to sibling containers.

The clean 2026-10-23 launch passed **24 contract/regression tests** and all
**145 real-cohort chr22 smoke tasks**, including 48 native FLARE target fits.
The independent smoke audit reconstructed 16 aggregates and 72 local predictions.
The full stage then started in
`runs/observed_v1_20261023_clean`, with named controller
`snpster-ancestry-observed-20261023_181805_879906705`.
Those are mechanics/provenance checks, not full-run completion or accuracy
certification. The real smoke also exposed and fixed a seven-source rounding
reader error: FLARE's +/-0.005 bound is per exported component, not per vector.

The full stage subsequently completed on **2026-10-24**: all **1,111 tasks**
succeeded in 26h 15m. Its independent audit verified **16 aggregates and 17,160
local predictions**, and the controller exited zero after the audit.
Approximately 12 GiB of disposable preparation/smoke/full work and their exact
Nextflow cache directories were then removed. Both result manifests passed a
second hash-only verification without original scientific inputs or work mounted.
Published inputs, fitted models, predictions, summaries, frozen code and execution
reports remain intact. Cache-only resume is no longer available for this
successfully cleaned run; failed/retired development runs were not cleaned.
These successful execution and durability checks do not change the declared
development-only, non-independent status of the biological evaluation.

## What to compare

Start with labelled **FLARE** (fixed and EM-estimated generations), a constrained
allele-frequency control, and genuine **Gnomix fast** as an explicitly
research-only comparator. Add PCA/linear, PCA/RandomForest and PCA/XGBoost
controls only as declared comparisons. More features and better pure-reference
accuracy did not consistently improve the repository's mixed/pedigree or
real-sample transfer results. See [the evidence and research notes](RESEARCH.md).

| Model family | Parameters | Training semantics |
| --- | --- | --- |
| `flare` | `gen`, `em`, `update-p`, `states`, `panel-weight`, `array`; optional `model_file`, `components`, `component_map` | Freeze the reference; run each target independently. FLARE estimates target-dependent parameters during inference and saves each target's native model. `model=` alone does **not** freeze EM. |
| `gnomix` | `window_size_cM`, `context_ratio`, `logistic_c`, `crf_c2` | Pinned genuine fast architecture: reference-only window LR, separate training-only CRF. Window bases are shared across CRF L2 variants. Other architectures are deliberately not claimed. |
| `nnls` | none | Simplex-constrained least-squares fit against reference-only allele frequencies; not a local-ancestry model. |
| `pca_ridge` | `n_components`, `alpha` | Reference-only standardized dosage PCA; global regression fitted on separate training people. |
| `pca_rf` | `n_components`, `n_estimators`, `max_depth`, `min_samples_leaf` | Same reference-only PCA, training-only regression. |
| `pca_xgboost` | `n_components`, `n_estimators`, `max_depth`, `learning_rate` | Same split discipline, multi-output boosted regression. |
| `pca_gmm` | `n_components`, **`k`**, `reg_covar` | Reference-only PCA/GMM clusters; reference-labelled responsibility-weighted cluster-to-source map. A diagnostic baseline, **not** FLARE2 clustering. |

Every family optionally accepts `reference_cap`, the number of reference people
per authoritative `fine` population. Selection is deterministic by sample ID.
The marker policy is frozen **before** capping and is shared by all candidates
within a route/chromosome. Hyperparameters are allowlisted; misspellings, global
`k` passed to FLARE, invalid probability axes and invalid numeric values fail.
Do not count `array=true/false` as distinct scientific settings when explicit
FLARE MAC/MAF filters are both zero: those settings are then effective aliases.

All estimators, external programs, numerical libraries and Nextflow workers use
one CPU each. Nextflow, not Python process pools, schedules their concurrency.
FLARE Java heap is 70% of the declared task memory. The whole reference panel is
materialized in memory; choose memory appropriate to the dense panel size.
This implementation prioritizes correctness and auditable adapters, not a
streaming multi-terabyte reference backend.

## Input contract

Use [the JSON template](example.json) or generate a complete synthetic example
with [make_smoke.py](make_smoke.py). Paths must be existing **absolute host paths**
and must be identical inside the controller and sibling task containers.
The launcher snapshots the config and task implementations; resume rejects a
changed config. Scientific inputs, tools and custom model files are SHA-256
hashed and checked; Gnomix also requires the clean pinned upstream commit.

The strict top-level config contains:

- `version: 1`, `seed`, `mode: "smoke" | "full"`;
- `sources`: unique, ordered, explicitly named source groups. Seven and eight
  sources are supported as different estimands; no EAS/SEA split is invented.
- `chromosomes`, `spans`: full half-open physical spans, keyed by chromosome;
  `full` requires autosomes 1 through 22 in order. Partial runs are mechanics
  checks, never genome-wide certification.
- `reference`, `training`, `maps`: chromosome-to-file maps. Reference and
  training VCFs must already have labelled, phased diploid GT. Maps are
  headerless four-column PLINK maps: chromosome, marker ID, cM, position.
- `routes`: named input masks/stages, each with `id`, `stage`, `vcfs`.
  Query VCFs contain the calibration, validation and test people.
- `samples`, `segments`: manifests below.
- `reference_filters`: common `min_maf`, `min_mac`, `minimum_markers`.
- `models`, `interpretations`: explicit candidate IDs/families/kinds/parameters;
  no hidden grid expansion. The number of scored route configurations is
  `routes x models x interpretations`, not just the number of base models.
- `tools`: `flare_jar`, `beagle_jar` when needed, `gnomix_root` when enabled.
- `allow_research_only`: explicit Gnomix opt-in.
- `independent_test`: a provenance declaration, **not** genetic certification.
  Previously inspected development people cannot become a fresh test cohort.

Sample manifest, tab-separated:

```text
sample_id	role	family_id	donors	source	track	fine
reference_1	reference	fam_ref_1	founder_ref_1	EUR	pure	French
training_1	train	fam_train_1	founder_train_1		mixed	
calibration_1	calibration	fam_cal_1	donor_cal_1a,donor_cal_1b		mixed	
validation_1	validation	fam_val_1	donor_val_1a,donor_val_1b		mixed	
test_1	test	fam_test_1	donor_test_1a,donor_test_1b		mixed	
```

All five roles are required and disjoint. Families and the transitive union of
declared founder donors cannot cross roles. Use canonical donor IDs, include the
person's own identity for pure samples, and identify every founder for simulated
descendants. Use connected-component family IDs from your genetic-relatedness
audit, not arbitrary per-row IDs. The workflow verifies **declared** separation;
it does not replace a kinship screen or an imputation-reference exclusion audit.
`track` is a reporting subgroup, e.g. pure/mixed/scenario; it is never a feature.

Realized local truth, tab-separated:

```text
sample_id	chromosome	haplotype	start	end	source
training_1	22	0	0	20000000	EUR
training_1	22	0	20000000	50818468	SSA
training_1	22	1	0	50818468	EUR
```

Haplotypes are **0/1**, coordinates half-open, source names exact. Both copies of
every non-reference person must tile every declared chromosome span, without
gaps or overlaps. Global truth is derived from these realized lengths, not
nominal family fractions. Training VCF phase must match training segment
haplotypes; do not rephase labelled training data without a homolog/truth audit.
No truth is inferred from model predictions or reference-imputed genotypes.

### Two production handoff routes

- `harmonized`: complete biallelic SNP GT, no DS substitution or missing-call
  padding. Phase-free controls consume the same diploid calls. When FLARE or
  Gnomix is requested, a **singleton, reference-assisted phase-only** Beagle
  adapter uses `impute=false`, checks exact keys and unchanged allele counts,
  and adds no variants. Phase-only work is cached separately. It is not a
  phasing-strategy search.
- `imputed`: already-phased imputation output; no additional phasing/imputation.
  The stage name is not evidence of imputation accuracy. Existing production
  genotype-concordance gates and donor/reference audits still apply.

VCF chromosome names must be numeric 1..22. Exact CHROM/POS/REF/ALT intersections
with the phased reference and labelled training VCF define the common panel.
Only complete diploid biallelic A/C/G/T SNP GT are supported. Missing GT fail
rather than being filled. Sample IDs and marker positions must be unambiguous.
Filter selection uses reference genotypes only, never evaluation labels.
Methods share the exact panel within each route. **Across input stages the
panels may differ**; report this as stage-native evidence, not a matched-marker
causal claim about imputation.

## Interpreting anc1 / anc2 / anc3

1. **Labelled FLARE/Gnomix:** use the actual model's component names/classes,
   not a guessed positional order. Native fractions already have meaning;
   fitting a decoder is optional, not mandatory.
2. **Unlabelled supplied FLARE model:** specify its stable `components` axis.
   An explicit row-stochastic `component_map` can map that axis to `sources`.
   Otherwise raw native components remain saved and explicitly **unscored**;
   never automatically rename `anc_0` to EUR. Use a separate calibration-only
   `simplex` head to learn a nonnegative row-stochastic component/source map.
3. **Optional learned heads:** `simplex` (`l2`), `ridge` (`alpha`),
   `rf` (`n_estimators`, `max_depth`, `min_samples_leaf`) or `xgboost`
   (`n_estimators`, `max_depth`, `learning_rate`) fit on **calibration people
   only**, after base training. RF/linear/XGBoost outputs use a documented
   Euclidean simplex projection; malformed native probabilities are not clipped.
   The native track is always retained separately.

Calibration people, their families and donors never train the base model.
Their native features are genuinely out of base-training sample, so this
five-role protocol does not require a second, potentially leaky in-sample OOF
shortcut. Hyperparameters are selected on validation people only. Test labels
never train heads or select finalists. Do not repurpose inspected validation
scores to repeatedly redesign candidates and still claim a fresh test.

FLARE's pinned writer rounds each local source probability to two decimal
places. Its mass error therefore depends on the number of sources; 0.005 is
not a valid whole-vector tolerance for seven sources. Each finite component
must lie in [0,1], and its clipped +/-0.005 rounding interval must admit a
joint unit-mass vector. Only then is explicit mass renormalization used for
integration. The raw mass correction and export policy are recorded; exact
unrounded posteriors cannot be recovered. Invalid axes,
negative/nonfinite values or larger errors fail. The scorer integrates local
posteriors using physical midpoint cells over complete spans, weights copies
equally, then weights chromosomes by span. It does **not** average `.global.anc`
or count markers. Gnomix windows receive the sum of the same physical cells.

The metric is `100 * 0.5 * sum(abs(predicted - realized))` percentage points.
Average masks within person, then average people equally within input stage.
Native/learned and permissive-software/research-only leaderboards stay separate.
Test metrics are reported for validation-selected finalists only. Raw test
predictions remain auditable but must not be used as another tuning leaderboard.
Pure/mixed validation subgroup errors are emitted; relatives are not independent
evidence. No automatic production promotion or commercial eligibility claim.

### Persistent three-generation benchmark

**Verified release, 2026-10-25:**
[pedigree_observed_v1_20261025_v3/full/release](datasets/pedigree_observed_v1_20261025_v3/full/release).
Its [dataset descriptor](datasets/pedigree_observed_v1_20261025_v3/full/release/configs/dataset.json),
[compact comparison config](datasets/pedigree_observed_v1_20261025_v3/full/release/configs/compact.json)
and [independent verification](datasets/pedigree_observed_v1_20261025_v3/full/release/verification/independent.json)
are ready for consumption. All 4 smoke and 68 full scientific tasks completed;
22 independent chromosome audits reconstructed 4,482,518,860 inherited GT
alleles. The release retains 2,710 people (1,765 reference plus the 945 below),
441 new mixed descendants and 3,873 individual/ancestor metadata records.
No model fitting has run. The final post-run checker needed a narrowly scoped
recovery because the original KING table also contains older simulated people
not included in this new cohort. The corrected frozen checker validates their
original connected components rather than treating them as new active people.
The original producer and failed post-run status remain unchanged;
`verification_recovery/` preserves the successful auditor and evidence.

`pedigree_dataset.json` specifies founder counts, sources, seed and three
successive halving generations. `run_pedigree_dataset.sh` consumes an independently
verified kinship-aware release and creates a new immutable dataset, not model
fits. The default allocation is:

| Role | Pure founders (G0) | G1 | G2 | G3 | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Training | 224 | 112 | 56 | 28 | 420 |
| Calibration | 56 | 28 | 14 | 7 | 105 |
| Validation | 112 | 56 | 28 | 14 | 210 |
| Test | 112 | 56 | 28 | 14 | 210 |

The remainder supplies one real reference representative per retained kinship
component. Each pure founder is used once in one tree; each successive generation
pairs the preceding generation without reuse, siblings mating or overlapping
founder donors. Different-source first-generation parents are required; later
pairings favor different expected ancestry vectors. Generations describe pedigree
depth, not a guaranteed count of ancestry sources. All eight founder donors of a
terminal family, their known biological relatives and all descendants remain in
one role. The allocation does not inspect model predictions. Legacy connected
components are conservatively preserved, even when their connections arose from
older simulated pedigrees. Source coverage is deliberately unequal because
independent OCE/NAM founders are limited.

```bash
bash data_analysis/ancestry_training/run_pedigree_dataset.sh \
  /absolute/verified-kinship-run/release \
  /absolute/path/to/pedigree_dataset.json \
  /absolute/new-pedigree-run
```

The pinned, named Docker controller runs a real chr22 simulation, export and
independent inheritance smoke before the separate 22-autosome Nextflow build.
The launcher prints each actual `docker logs -f` command. Its default budget is
12 one-core / 8-GiB tasks, bounded by 96 GiB plus the 2-GiB controller. Sources,
the design, hashes, raw Ped-sim breakpoint evidence, trace, DAG, reports and
exit statuses persist. Failed work/history are retained.

The consumer root is `full/release/`:

- `genotypes/all/chrN/master/all.phased.vcf.gz`: indexed phased GT for every
  retained reference, pure founder and simulated descendant. These persistent
  observed-union SNP VCFs are reusable scientific data, never disposable work.
- Each master also retains `inheritance.tsv`, ancestry `segments.tsv` and
  original `pedsim.bp` / `pedsim.ids` evidence.
- `genotypes/pipeline/chrN/bundle/`: phased reference/training VCFs and
  phase-stripped, unchanged-GT query VCFs for both original array masks.
- `metadata/samples.tsv`: role, pedigree-connected family, exact founder donors
  and `pure` / `mixed_g1` / `mixed_g2` / `mixed_g3` reporting tracks.
- `metadata/individuals.json` and `.tsv`: G0--G3 identity, biological or simulated
  parents, parent ancestry, expected fractions, physically realized fractions
  and old/new relatedness components. Excluded ancestors are metadata, not extra
  evaluated people. JSON fraction values are 0--1, not percentages.
- `configs/{dataset,compact,sweep,smoke}.json`: portable descriptor and model-input
  configurations. Pipeline configs use absolute local paths; the descriptor's
  data paths are relative to the release.
- `verification/`: independent exact founder-haplotype reconstruction at every
  GT allele, Mendelian inheritance, ancestry tiling, family separation, hashes
  and actual pipeline `freeze` / both-mask `prepare` checks.

Selection additionally reports validation and selected-finalist test mean,
median and p95 person TV error, per-source MAE/RMSE/bias and separate pure/G1/G2/G3
errors. A deterministic 2,000-replicate family-component bootstrap provides a
95% mean-error interval when at least two family units exist; it is unavailable
(`null`) with only one. Test results never enter selection.

The 210-person validation/test sets each have only 14 independent simulated
terminal trees, not 210 independent families. Shared founders, chromosomes and
masks do not create independent people. Reused donors are not a fresh external
test; operational pure labels are not individual admixture measurements. Three
halving generations do not certify rare ancestry below nominal 12.5%, and sparse
OCE/NAM diversity still limits conclusions. No imputation or model fitting is
performed by this dataset build.

## Running and resuming

The launcher uses the locally pinned Docker images in
[nextflow.config](nextflow.config), the cached **Nextflow 26.04.4**, a named TTY
controller and sibling task containers. It needs the Docker socket and the
existing local image/framework cache; it intentionally runs without networking.
The code/config snapshot and exact input/tool hashes are preserved.

```bash
# Start with a bounded smoke, not a large training sweep.
SNPSTER_MEMORY_GB=2 SNPSTER_TASK_MEMORY_GB=2 \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/experiment.json /absolute/new-run

# The launcher prints the ACTUAL name:
docker logs -f snpster-ancestry-standard-<printed-timestamp>

# Same snapshot, same inputs, same dedicated history; fresh report filenames.
SNPSTER_MEMORY_GB=2 SNPSTER_TASK_MEMORY_GB=2 \
  bash data_analysis/ancestry_training/launch.sh \
  /absolute/experiment.json /absolute/new-run --resume
```

Default task budget: **one CPU / 8 GiB**, plus a **2-GiB controller**.
For a larger run, explicitly set `SNPSTER_CPUS`, `SNPSTER_MEMORY_GB` and suitable
`SNPSTER_TASK_MEMORY_GB`; `SNPSTER_PCA_MEMORY_GB` optionally raises only PCA
base-task memory. The focused comparison above uses 22/160 with 6-GiB tasks
and 16-GiB PCA bases. This is a task
allocation budget, not an aggregate Docker cgroup across sibling containers;
leave headroom for the controller, host, and other experiments. Do not allocate
the same full-machine budget to multiple workflows.

Repository and `/srv/dependencies` inputs are read-only mounted. Other input
directories must be listed, one absolute directory per line, in
`SNPSTER_INPUT_ROOTS_FILE`. Mount paths deliberately reject shell metacharacters
and spaces. Always use **identical absolute paths** in fixture generation too;
do not write `/data/...` into a host-path config unless `/data` truly is its
host path.

Results are copied, not symlinked back into work:

- `frozen/`: source-independent config/plan/input-hash manifests;
- `panels/`, `phased/`: exact feature keys, identities and phase audits;
- `bases/`, `models/`: reference artifacts and fitted serialized heads;
- `predictions/`: per-chromosome native outputs, per-target FLARE models/logs,
  or Gnomix window posteriors;
- `aggregated/`, `interpreters/`: full-scope representations, fitted mappings
  and separate native/learned probabilities;
- `summary/`: validation table, selected-test metrics, raw prediction evidence;
- `verification.json`: separately recomputed headline scores.

The run directory also retains immutable task sources, controller logs and exit
status, each invocation's launcher source/hash, trace, graphical DAG, HTML
report/timeline and plain Nextflow log.
After Nextflow succeeds, the launcher automatically runs the independent durable
audit below and exits unsuccessfully if that audit fails.
Pickles are trusted local artifacts; never load untrusted model files.
Models are panel-specific research artifacts, not a new production serving API.

### Verify durable outputs and clean exactly this successful run

Run [audit.py](audit.py) in the pinned task image with the same read-only
scientific-input mounts and a writable result mount. It independently
reconstructs physical truth, native local integrations, chromosome aggregation
and person-level metrics, then hashes **all** published artifacts. It can run
without Nextflow work:

```bash
python3 audit.py /absolute/new-run/results
# After deleting only this successful run's work/cache:
python3 audit.py /absolute/new-run/results --manifest-only
```

The launcher never globally prunes. Keep failed work/history for recovery.
After successful workflow, independent durable audit and resume validation,
remove only the exact run's `work` and `launch/.nextflow/cache` directories.
Keep scientific inputs, source/hash manifests, final models/predictions and
reports. Removing the cache deliberately ends resumability of those tasks.

## Tests and limitations

Use the pinned task image, `python3 -m unittest -v test_workflow`, from this
directory. [make_smoke.py](make_smoke.py) generates 256 synthetic chr22 SNPs,
24 reference people, eight training and twelve held-out people, both input
routes, all implemented model families, and all decoder heads. Its two-source
order is deliberately reversed to exercise native axis mapping.
This is mechanics-only evidence, **not an ancestry accuracy benchmark**.

Validated locally with the pinned runtimes:

- 22 contract/regression tests passed.
- Both routes completed a real 149-task Nextflow smoke: eight model
  configurations, sixteen route/model aggregates and eighty interpretation
  configurations, including real Beagle, FLARE and upstream Gnomix execution.
- Phase-only processing added zero variants and changed zero diploid calls.
- Two CRF configurations shared their fitted reference base, retained actual
  fitted `c2=0.1` versus `c2=1`, and produced different probabilities.
- Independent reconstruction checked sixteen aggregates and seventy-two local
  predictions; an exact-run resume reused all 149 tasks from cache.
- Only successful smoke work/cache was removed, after verification. All retained
  artifact hashes still verified without mounting original inputs or work.

Excluded intentionally: phenotype/PGS modification, production imputation
tuning, alternate phasing strategies, automatic pedigree simulation, automatic
biological naming of clusters, genetic kinship certification, automatic
reference-data licensing clearance, neural/new Gnomix architectures, and
training a new FLARE2 copying-profile GMM from scratch. The latter can be
researched separately; a compatible explicitly mapped native model file is
already accepted.
