# SNPster v3: 100-model unified pedigree sweep + robust stage-2 recovery

## Current benchmark contract: fixed production preprocessing (2026-10-22)

The [version 2 contract](./config/ancestry_benchmark_contract_v2.json) supersedes
the flexible-preprocessing contract for new benchmarks. Existing experiments,
models and their frozen protocols are unchanged. The new production-matched
workflow is implemented and the all-source donor gate has passed. A real
chromosome-22 smoke exposed insufficient observed-marker coverage in the
historical backbone. Broader original-assay recovery completed, but independently
audited imputation concordance exposed artificially thinned inputs and
reference-overlapping training controls. Dense observed-data recovery and
all-autosome genetic screening have now completed; no new model scores are
available yet.

### Approved dense-split synthetic development sweep

The user-approved [dense split policy](./config/dense_split_protocol.json) reuses
the existing [production-matched Nextflow sweep](./production_matched_benchmark.nf),
not a separate orchestration framework. It starts from the requested 3,033
reference people, 240 pure-training people and 76 descendants from 25 complete
training families. The previously inspected 14 controls are retained, with a
pure SEA control added from the already held-out founder HGDP00720. This produces
15 development controls: eight pure and seven pedigrees.

The completed prediction-independent freeze is
`results_v3_dense_split_20261022_211648/freeze/cohort/frozen/cohort.json`.
All-autosome genetic screening removes ten reference people and one training
founder related to evaluation donors. Removing that founder's whole training
families leaves **3,023 reference people, 239 pure-training people and 56 training
descendants from 20 families**. These changed counts are explicit in the frozen
manifest. The 67 eligible observed SGDP donors remain an external holdout,
not a tuning pool.

This epoch deliberately uses **synthetic genomes defined by frozen dense
reference haplotypes**, as explicitly approved by the user. GT and phasing
confidence do not establish independently observed sequencing-call truth.
Results therefore describe synthetic development performance, not independent
real-person accuracy. This does not weaken the separate observed-data policy
below or relabel reference-imputed calls as observed SGDP evidence.

Original Ped-sim plans and physical breakpoints are reused. Both distinct founder
homologs and native eight-source labels are preserved. Simulated microarray calls
at real rsID-backed mask loci are assembled separately from untyped validation anchors, avoiding
the old common-truth thinning bug. Every mask/chromosome must retain at least
90% of its in-span, in-map keys. Existing per-person genotype concordance,
harmonization, split and QC gates still apply.

Benchmark-only imputation copies exclude simulated training/evaluation founders
and their detected close relatives. Every retained allele key and phased genotype
is independently stream-hashed against the original reference. The production
reference files, production images and phasing/imputation methods are unchanged.
Ancestry fitting references exclude evaluation donors and detected relatives.

Gnomix's base model uses the large fitting reference. Its auxiliary CRF smoother
uses the retained pure and pedigree training people, not the evaluation controls.
Pedigree labels are mapped onto production-inferred homologs using training-only
concordant heterozygotes and nearest informative physical positions; uninformative
mixed homologs fail rather than receiving invented phase. Upstream modal
SNP-window labels are used. The alignment policy and marker support are retained.
Neither FLARE nor the NNLS controls uses this auxiliary CRF-training pool.

Launch with `bash run_dense_split_sweep.sh ABSOLUTE_FROZEN_EXPERIMENT`. The
[launcher](./run_dense_split_sweep.sh) first checks all ten pinned runtime images,
verifies chr22 reference copies and runs a complete chr22 production/model smoke.
Only after independent verification does it extend to all 22 autosomes and both
masks, subject to a measured disk-capacity gate. Graphical Nextflow progress,
DAGs, reports, timelines, plain logs and traces are retained. Successful caches
are cleaned by exact session only after durable verification and then rechecked;
failed work remains resumable.

The same ten configurations and the same TV error apply to every person:
physical-length chromosome aggregation, then mask averaging within each person,
then equal-person averaging. Both native-eight and EAS+SEA-merged-seven errors are
reported. Pedigree truth uses realized physical ancestry, **not equal weighting
of a person's ancestor donor list**. No production promotion is authorized.

### Approved 98-configuration model-only expansion

The [expansion policy](./config/dense_model_expansion.json) freezes **98 distinct
effective settings**, rather than repeated random seeds:

- **64 FLARE:** full reference or fine-population caps 8/16/32, EM false/true,
  and generations 5/7/10/15/20/25/35/50.
- **32 fast Gnomix LR+CRF:** windows 0.25/0.5/1/2/3/4/6/8 cM crossed with CRF
  L2 penalties 0.01/0.1/1/10. A shared reference-only LR base is fitted for each
  chromosome/window; each mask/regularization reloads that untrained-smoother
  base independently. The existing default penalty is effectively 1.
- **Two unchanged NNLS controls:** imputed and raw-observed.

Pre-launch source inspection found that FLARE's array toggle only selects its
minor-allele-count filter. With the pipeline's explicit `min-maf=0,min-mac=0`,
true/false are the same effective setting. Additional generation values replace
that redundant axis. The original ten outputs provide eight unique baselines;
both redundant array variants remain documented as aliases, not counted twice.
Their ancestry predictions must match the canonical baseline exactly before
reuse. Existing source outputs and their original policies are unchanged.

Launch with `bash run_dense_model_expansion.sh ABSOLUTE_SOURCE_EXPERIMENT`.
The [model-only Nextflow extension](./dense_model_expansion.nf) waits for the
existing controller to finish successfully, verifies its durable outputs and
post-cleanup hashes, then reuses preprocessing and eight unique baseline fits.
Only 90 additional settings are fitted. Original completion records are retained
alongside explicit reuse provenance. The existing frozen run is not interrupted.

The cohort, founder exclusions, raw masks, production phasing/imputation, marker
and genotype-quality gates, physical inheritance truth, equal-person TV metric,
and both eight-source and EAS+SEA-merged reporting remain unchanged. The 67 SGDP
people remain external to tuning. Results remain **synthetic development
comparisons**, not independent real-person certification.

All 98 settings must pass chromosome 22 before whole-genome promotion. A
disk-capacity gate uses the actual new panels, models and shared bases; unchanged
scientific directories are linked, not regenerated. Each expanded scope receives
independent posterior/NNLS, inheritance and scoring checks. Exact-session caches
are removed only after durable verification, then retained hashes are rechecked.
The named Docker controller exposes graphical progress with `docker logs -f`;
DAG, report, timeline, trace, plain log and input/source hashes are retained.

The imputation-based expansion is not applicable when the observed-only protocol
is selected. Its per-person untyped concordance gate must not be weakened to
force model fitting. Use `bash run_observed_dense_model_expansion.sh` instead.
That workflow reuses the completed observed-only input and phasing certificates,
runs FLARE with `impute=false`, and evaluates the 64 FLARE plus 32 genuine Gnomix
settings. The two phase-free controls from the 98-setting plan are intentionally
omitted because one consumes imputed genotypes and neither is a native-method
candidate under the approved observed-only protocol.

Resume a retained observed-only expansion with
`bash resume_observed_dense_model_expansion.sh ABSOLUTE_EXPERIMENT ABSOLUTE_SMOKE_RECOVERY`.
The recovery directory must be the retained failed smoke task's `result`
directory. Complete FLARE candidates are accepted only after their full output
sets and exact frozen commands pass validation; incomplete Gnomix work is refit.

The observed-only expansion uses one cacheable preparation task per
chromosome/mask and one model task per candidate. The smoke therefore contains
96 model tasks. The full sweep contains 4,224 model tasks
(22 autosomes x 2 masks x 96 candidates), rather than 44 monolithic tasks that
each hide a sequential 96-model loop. Model tasks request one CPU and 8 GB;
the local executor permits 22 CPUs and 180 GB, allowing up to 22 model tasks
concurrently while leaving four GB below the memory ceiling. FLARE uses one
native thread and a six-GB Java heap. Preparation remains separately cacheable
and is limited to four concurrent tasks.

### Dense observed-support repair

The [data-only repair workflow](./dense_public_support.nf), launched with
`bash run_dense_public_support.sh --run`, reuses 13 verified public SGDP people
and recalls 55 additional FullyPublic/X people from sequencing PL/AD evidence.
The existing Chane person is also recalled as the original-call control.
The metadata census contains 68 candidates: EUR 16, EAS 2, SAS 15, SEA 8,
GME 8, NAM 11, OCE 2 and SSA 6. This is a candidate inventory, not an
independent genetically screened or frozen training/evaluation cohort.

The [repair policy](./config/dense_public_support.json) binds the original
GRCh37/GRCh38 FASTAs, indexes and chain. A real chromosome-22 Nextflow smoke
must pass independent PL/AD rederivation, strict duplicate handling, original
raw-call concordance, normalization verification and >=90% genuinely observed
coverage of **both actual input masks**, before autosomes 1-21 run.
Coverage eligibility is intersected across all 22 chromosomes: different
passing people on different chromosomes cannot satisfy the requirement.
No release-imputed GT/DS/GP is accepted as observed truth.

The named Docker controller provides the live graphical Nextflow overview
through the printed `docker logs -f` command and retains the pipeline DAG,
HTML report, timeline, plain log, trace and source hashes. Only downloads have
network access. Each successful scope is independently hash-verified before
exact-session cache cleanup and rechecked afterward; failed work is retained.

This stage authorizes **no ancestry models**. All-autosome genetic alias/family
screening, independent training/evaluation roles, a new prediction-independent
cohort freeze and decoupling per-person raw observations from common validation
truth remain mandatory before production preprocessing and ten-model testing.
The production phasing/imputation strategy and existing genotype-quality gates
remain unchanged. In particular, only two public OCE candidates are available;
the old two-training-plus-one-evaluation design is not yet feasible.

The completed repair is `results_v3_dense_public_support_20261022_152900`:

```bash
docker logs -f snpster-ancestry-v2-dense-support-20261022_152900
```

The earlier `20261022_152650` launch stopped before scientific execution
because Nextflow automatically mounted the repository and the launcher also
requested that exact destination. Its launcher, logs and failed cache are
retained. The corrected launcher uses scoped read-only code/data mounts,
avoiding the automatic repository-mount collision.

All 68 candidates passed both actual masks across all 22 autosomes; the lowest
individual/chromosome/mask coverage was 92.24%. The 149 full-scope Nextflow tasks
completed, durable published evidence was independently verified, exact-session
caches were cleaned and retained hashes passed the post-cleanup check.

The next [dense donor genetic screen](./dense_donor_qc.nf), launched with
`bash run_dense_donor_qc.sh --run`, screens the 68 candidates plus the three
previously excluded public Dinka people as **historical controls, not eligible
donors**. The comparison universe includes the fixed production reference,
all historical ancestry fitting/development identities and all public-public
pairs. No ancestry predictions select or exclude a donor.

Only exact common identity-marker alleles are merged from the already verified
observed sources. Original missing calls remain missing. A second implementation
compares every merged genotype and the repeated Chane control to its sources;
every fixed-reference, historical and public genotype in the final identity
panels is then independently rechecked. The unchanged all-autosome LD-pruning,
KING, pair-count and callability thresholds are reused.

Candidates in fixed-reference-connected family components are excluded.
Historical non-reference identity matches are reported as identity bridges,
not additional independent people or automatic authorization to reuse old
fitting roles. Such historical identities must be absent from the new fitting
roles. Candidate family components constrain the future donor split, and
conflicting duplicate identities/source labels require review. This screen
does not certify source purity or freeze/authorize ancestry-model execution.
It runs in a named graphical Docker/Nextflow controller with a real chr22 smoke,
reuses that verified chromosome in the all-autosome screen and hash-checks
durable copies before and after exact-session cache cleanup.

The current genetic controller is:

```bash
docker logs -f snpster-ancestry-v2-dense-genetic-20261022_192734
```

Its chromosome-22 smoke passed, including 100% recalled/reused Chane concordance
across 1,871 called identity markers and independent checking of every comparison
genotype. All 22 autosomes completed, independently verifying 903,171,846 genotype
cells and 329,866 candidate-involving pairs on 66,373 LD-pruned markers.
One reference-related Biaka person was excluded, leaving 67 eligible candidates.
Durable evidence passed pre/post-cleanup hashing. The earlier `192517`
probe's missing historical-chain mount was corrected by declaring that provenance
resource as a Nextflow input; its failure evidence is retained.

### Eight-source production-matched development comparison

The [frozen ten-configuration protocol](./config/production_matched_protocol.json)
and [cohort-freezing workflow](./production_cohort.nf) now select eight eligible
evaluation founders: one each from EUR, EAS, SAS, SEA, GME, NAM, OCE and SSA.
All 52 potential evaluation donors, including the 13 eligible public donors,
were screened against the same 4,627-person genetic comparison panel:
239,226 candidate-involving pairs. No evaluated candidate had a flagged
relative in the fixed imputation reference. One related fitting-role person,
`Murut15`, was excluded. Selected evaluation founders, the 16 pure training
people and the 3,033 ancestry-reference people are role-disjoint.

The development cohort has eight pure people and eight distinct-source F1
crosses forming an ancestry ring, each run with the Lonni and Frederik masks.
All 16 evaluation people share founders within one connected founder-family
component. This is a small development comparison, not 16 independent genomes,
an independent final test, a family-bootstrap confidence interval or comprehensive
pedigree certification. Pure-source truth uses the frozen operational population
labels, not a claim that these populations have no ancient admixture.

The [production-matched workflow](./production_matched_benchmark.nf) executes
the unchanged production standardizer, combiner, harmonizer, merger, Beagle
imputer, splitter and per-person PGS QC. Ancestry consumes the phased-and-imputed
branch before per-person PGS QC, while the unchanged QC branch is executed and
audited separately. Each imputation batch contains 32 people: 16 evaluation and
16 training people for one mask. Repeated masks of one person are never pooled
in the same batch. Raw uploads use real rsIDs from the fixed production PVAR,
not invented coordinate-based identifiers.

Evaluation simulation truth contains only genuinely observed founder genotypes. Alleles
must occur exactly in every required observed source and the fixed reference;
unavailable or uncalled loci are explicitly reported and excluded, never filled
from imputed calls. The original >=90% truth-locus retention gate is unchanged.
Production founder phasing must preserve the original observed diploid GT;
pure people inherit complementary homologs, and F1 people inherit one homolog
from each distinct parent.

Before ancestry fitting, every person must pass observed-GT concordance >=99.5%,
untyped-GT concordance >=90%, minimum anchor counts, typed-marker retention and
unchanged split/QC haplotypes. Pooled concordance cannot hide a failing person.
A global gate requires the complete declared chromosome-by-mask preprocessing
scope before any model executes.

The ten frozen candidates are six FLARE configurations (full/fine-capped reference,
array-informed/non-array-informed settings and fixed/EM settings), genuine
Gnomix LR/CRF at 0.5 and 2 cM windows, and two uncalibrated constrained
allele-frequency controls using respectively production-imputed and genuinely
observed calls. These controls are not ADMIXTURE. All predict eight native
sources; no model invents an EAS/SEA split from a seven-source output.
Gnomix remains academic-research-only and no automatic production promotion
is performed.

Launch with `bash run_production_matched_benchmark.sh`. The controller first
executes all ten configurations on chromosome 22 / Lonni, independently checks
Mendelian inheritance, local ancestry aggregation or convex optimality, binds
durable copies, cleans only the verified successful session and rechecks hashes.
Only then does it launch all 22 chromosomes and both masks. Smoke scores are
not whole-genome ancestry errors.

The initial experiment is `results_v3_production_matched_20261022_112337`.
Its first smoke stopped on 37 unavailable exact input loci, before model fitting.
The observed-source intersection was corrected without changing the production
pipeline or relaxing genotype/retention gates; original failed sources, hashes
and controller logs are retained. The resumed controller reused the completed
rsID preparation task and completed real production founder phasing. The next
gate correctly stopped: the 1,840 retained truth loci on chromosome 22 overlap
only 75 Lonni markers (314 Frederik markers), below the required 100 typed loci.
The original wider four-dataset backbone still supplies only 76 Lonni markers,
so simply undoing canonical filtering cannot fix this. No ancestry model ran.

The [observed-support recovery workflow](./production_truth_support.nf), launched
with `bash run_production_truth_support.sh`, instead uses the broader original
GSE73996 unphased calls and original GSE77508 assay matrix already present locally.
GSE73996 alone contains 407 exact Lonni chromosome-22 markers, versus 75 in
the retained historical truth panel; complete shared, callable coverage must
still pass on every chromosome and both masks. This is a prediction-independent
input repair, not a relaxed gate or a denser-mask-only leaderboard.

GSE77508 raw A/B calls require score >=0.8. Allele orientation is inferred using
at least four informative **other** people and checked against published
genotypes; the evaluation donor cannot determine its own orientation. Original
missing/low-quality assay calls are excluded rather than copied from phased
or imputed output. Unresolvable homozygous A/B orientation fails eligibility;
an original heterozygote needs no arbitrary orientation. Training-only genotypes
can come from their declared fixed dense reference role; they are not independent
evaluation truth.

The first broader-support attempt stopped on a QC-accounting error, not an
absolute shortage of SNPs. Chromosome 4 had 9,693 source-compatible loci and
8,201 complete observed loci, including 1,546 Lonni and 7,099 Frederik markers.
The raw-assay checks excluded 1,280 loci before genotype-completeness assessment.
The previous gate incorrectly counted these exclusions against its 90%
genotype-completeness requirement: 84.61% before assay QC versus 97.48% among
the 8,413 assay-eligible loci.

The corrected check applies the unchanged score >=0.8 and allele-orientation
rules first, then requires >=90% complete genotypes among those eligible loci,
plus >=100 complete typed markers for **each** mask. It does not restore any
excluded or missing calls, lower the score cutoff, or change which complete
genotypes are eligible as truth. Every chromosome now writes an
`eligibility.json` audit, including both denominators, per-person called counts,
exclusion reasons and mask counts, before a gate can fail.

The recovery completed all 22 chromosomes and wrote a new derived cohort with
additional observed-input hashes. Durable outputs were verified, the exact
successful support cache was cleaned, and published/source hashes were checked
again. The original frozen protocol, source labels and donor choices remain
unchanged. The automatic child benchmark then stopped before its original
standardizer executed: that production image lacks `ps`, required by Nextflow
resource monitoring. This was an instrumentation failure, not a genotype or
ancestry-model result.

The [monitoring preparation helper](./prepare_nextflow_metrics.sh) now extracts
genuine procps and its loader/libraries from the pinned existing preparation
image. The [wrapper](./pedigree_benchmark/scripts/nextflow_ps.sh) uses those
private libraries only for `ps`; read-only mounts supply it without rebuilding
production images or changing application libraries, phasing or imputation.
All assets are hash-bound. The [runtime smoke](./production_runtime_smoke.nf)
passed actual Nextflow CPU, memory and I/O collection in all ten unchanged task
images, including genuine Gnomix and FLARE. Its exact successful cache was
cleaned and retained evidence rechecked.

The next attempt completed all 32 standardizers but stopped in the harmonizer:
its original root-owned PLINK executables are mode 0664. The application tries
to add execute permission itself, which the benchmark's unprivileged user cannot
do. The earlier monitoring probe did not exercise these subprocess binaries.
The [executable preparation helper](./prepare_production_executables.sh) now
copies both original binaries, verifies their bytes against hashes computed
inside the pinned harmonizer image, and sets only the copies' execute bits.
Read-only mounts apply **only** to the harmonizer; no binary version, application
code, scientific setting or image changes. These copies and their original-image
evidence are included in the source manifest.

The expanded runtime smoke now invokes the required PLINK, PLINK2, bcftools,
Java/Beagle, bgzip and tabix tools in their corresponding original images,
checks the splitter plugin, opens the actual indexed reference FASTA and checks
model runtime imports. All ten checks passed before the new production recovery.
The failed sources, original permissions/hashes, task logs and controller
evidence remain preserved.

The graphical recovery reused the six completed production preparation stages,
all 32 standardizers and 31 completed combiners; real harmonizer tasks then
completed successfully past the former permission error. It
retains the preprocessing gates, automatic smoke/full continuation and
independent verification before successful-session cleanup. Live ANSI progress
uses a controller TTY; each run also saves a whole-pipeline HTML DAG, report,
timeline, plain Nextflow log and task trace under its attempt directory.
The corrected run completed harmonization, merging, Beagle imputation, splitting
and all 32 QC executions. It then stopped at the **unchanged genotype-quality
gate**, before any ancestry model. Independent exact-allele/sample-order
extraction reproduced the failure on 891 untyped chromosome-22 truth loci:
the 16 evaluation people have 62.63--75.87% genotype concordance (mean 68.20%),
below the frozen 90% per-person minimum. Their observed-call concordance is
99.64--100%.

The shared truth support supplies only 279 of the frozen Lonni mask's 2,245
chromosome-22 loci (12.43%), so the synthetic input is substantially thinner than
that mask. All 16 training donors also occur in the unchanged imputation
reference; their 97.08--100% untyped concordance is therefore **not** an
independent imputation-quality control. None of the evaluation donors occurs
there. These observations identify input-density/training-representativeness
concerns, not ancestry-model performance or proof of a particular causal
mechanism. The [per-person audit](./results_v3_observed_truth_support_20261022_114223/benchmark/failures/imputation_concordance/production_imputation_gate_audit.json)
and reproducer are retained with the failure evidence.

The run is **stopped**, not active. No quality threshold has been weakened and
no model scores are claimed. Scientific caches remain resumable while the
input-support problem is addressed; the successful runtime-probe cache was
already verified and cleaned. Inspect the latest controller's preserved log:

```bash
docker logs -f snpster-ancestry-v2-production-recovery-20261022_125800
```

The [resume launcher](./resume_production_matched_benchmark.sh) refuses concurrent
controllers for the same experiment, preserves source-hash epochs and uses the
exact Nextflow session. For the nested expanded-support benchmark, provide its
original derived cohort using `SNPSTER_MATCHED_COHORT`; it must match the frozen
experiment copy byte-for-byte. The [fresh launcher](./run_production_matched_benchmark.sh)
defaults to that completed expanded cohort and shares the same validated
graphical/monitoring runner. Existing experiment directories cannot be silently
overwritten. Failed task evidence is retained rather than broadly
pruned. Focused tests cover contracts, real VCF/PSAM handling, explicit
input-locus exclusions, original assay quality/orientation, per-person gates,
post-assay completeness accounting, convex optimality and durable-copy tamper
detection before and after cleanup. The retention recovery uses the same exact
Nextflow session and reuses the completed assay extraction; its source/hash
epochs and original failure evidence are retained.
The completed donor-cohort cache was independently checked against all 239,226
raw genetic pairs and six durable cohort files, cleaned to 4 KiB, and its
published/source hashes rechecked afterward. Failed model-smoke evidence remains
resumable and is not broadly pruned.

Use one error for every person, pure or mixed:
`100 * 0.5 * sum(abs(predicted_fraction - true_fraction))`.
Thus 100% EUR predicted as 90% EUR + 10% EAS is 10% off. Equal thirds
EAS/SEA/EUR predicted as 40% EAS + 40% SEA + 20% SSA is 33 1/3% off,
not 27%. Inputs must be normalized finite, nonnegative fractions; do not
silently normalize invalid inputs. Sparse vectors in the contract's examples
omit zero components only; actual prediction tables require every source.

Predict eight sources and report the same TV error at both eight- and
seven-source resolution, deriving the latter by adding EAS+SEA in prediction
and truth before scoring. The primary objective is eight-source mean person
error on a fixed, predeclared cohort spanning pure and mixed individuals.
Average the complete predeclared mask set within each person first, then
average unique people equally. Do not use different pure/pedigree objectives,
legacy scenario/family-balanced scores or post-hoc source weights. Show
subgroups and the error distribution as diagnostics. Seven-source-only
historical models remain seven-source baselines; do not fabricate an
EAS/SEA split or claim a smaller merged error is an eight-source improvement.
Mean error <=10% is acceptable and <=5% desirable; neither promises that
every person is within the target.

Keep the existing [production preprocessing](../../../data%20pipeline/)
unchanged, including its resolved Beagle command, reference, maps, batching,
merge/split behavior and QC. FLARE/Gnomix must use actual production
phased-and-imputed output; compatible phase-independent controls can use
its harmonized output. Historical singleton phase-only and inherited-phase
tracks are diagnostics, not selection inputs for this new contract.
Only the ancestry consumer's placement/settings may change.

Before scientific execution, freeze the resolved runtime/input identities,
audit evaluation founders against the fixed production imputation reference
and ancestry/training donors, and freeze source labels, case mix, masks,
truth and family-separated splits. Select eligible founders rather than
changing production imputation to create a disjoint benchmark. Source hashes
in the contract identify inspected entry points, not completed runtime
preflight or evidence that enough independent founders exist.

The initial training-only comparison is limited to ten configurations,
prioritizing FLARE with genuine upstream Gnomix and the phase-independent
control. Candidate settings, panel policies and fitting budgets need a
separate freeze. Complete all 22 autosomes, aggregate fractions by physical
length before calculating person-level TV, and independently verify durable
results before accepting a global score.

Future Nextflow **orchestrators themselves** must run in a named Docker
container, not only spawn short-lived task containers. Report its actual
name and `docker logs -f <container_name>` at launch, with overall progress,
failures, verification and cleanup streamed to stdout/stderr. Pin controller
and task runtime identities, test real task-container execution and persist
launch/history/work/output outside the controller. Retain controller logs,
exit status and failure evidence. Only remove exact-run disposable work/cache
after independent verification, and recheck published results after cleanup.

[Contract acceptance tests](./tests/test_benchmark_contract_v2.py) check the
numerical examples against the existing eight-source TV implementation,
named EAS/SEA merging, equal-person versus legacy scenario weighting,
production source identities and the required pre-launch gates. They do not
certify a new scorer, container launcher or scientific dataset.

### Historical production preflight: missing NAM, now resolved above

The [containerized preflight](./production_preflight.nf), launched with
`bash run_production_preflight.sh --run`, completed the runtime snapshot,
all 22 production-reference audits and cohort/source-label checks on
2026-10-21. The scientific-readiness gate then deliberately exited nonzero.
This is a verified data-eligibility blocker, not an ancestry-model result.

The unchanged production imputer is image
`sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1`,
with Beagle 5.5 (`beagle.17Dec24.224.jar`), Java 21.0.12, heap 6 GB and
one Beagle thread. Its embedded Python sources match the inspected
production sources. All 22 reference headers contain the same ordered
4,091 identities; reference, index and map byte hashes are retained.
This does not complete cohort-specific harmonization/batching/QC freezing.

All 19 evaluation founders in the recent dense pilot occur in this
production **imputation** reference, despite being disjoint from the
previous ancestry inference reference. Reusing them would not provide
an independent production-matched benchmark.

The original array datasets retain these directly ID-disjoint people:

| Existing role | EUR | EAS | SAS | SEA | GME | NAM | OCE | SSA | Total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pure development | 4 | 5 | 4 | 18 | 6 | 0 | 1 | 1 | 39 |
| Pure training | 8 | 10 | 6 | 29 | 10 | 0 | 2 | 2 | 67 |
| Ancestry reference | 52 | 58 | 33 | 180 | 67 | 0 | 10 | 14 | 414 |

These are direct-ID exclusions, not completed alias/relatedness audits.
There is no ID-disjoint NAM donor in any of these three pools, and
OCE/SSA development coverage is very limited. Reassigning these existing
roles cannot resolve the missing NAM source. Do not silently change
production imputation, use overlapping founders, or relabel a limited-
coverage study as full seven/eight-source validation.

The retained controller is
`snpster-ancestry-v2-preflight-20261021_193442`; its logs remain available:

```bash
docker logs -f snpster-ancestry-v2-preflight-20261021_193442
```

Nextflow 26.04.4 ran inside the named production-image container and
launched pinned sibling task containers through the Docker daemon.
The controller has exited with status 1 from the explicit readiness gate;
no ancestry sweep is running. Launch/history, trace, controller logs,
published audit and approximately 4.6 MiB of failed-run task work are
retained; the normal success-only cleanup rule does not apply.
Nextflow warned that execution HTML/timeline rendering failed during
shutdown; these are not evidence of a successful scientific run.

See the [published audit](./results_v3_production_preflight_20261021_193442/summary/preflight.json)
and [independent verification](./results_v3_production_preflight_20261021_193442/preflight_verification.json).
The verifier separately re-queried all 22 reference headers and reconstructed
the exclusion counts from original metadata; it did not independently
rehash all reference bytes. Both EAS/SEA-separated and merged TV reporting
remain required once a suitable cohort/design has been approved.

### Public independent-donor search and SGDP genotype smoke

The user selected finding suitable public independent data rather than
changing production imputation or silently accepting missing NAM coverage.
The official [SGDP distribution](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/)
provides metadata containing original sample aliases and explicit access
categories. Initial metadata checks identified 11 NAM candidates from
Mixe/Chane/Mixtec/Zapotec/Quechua, two Australian OCE candidates and three
Dinka SSA candidates absent under all listed aliases from the fixed
production panel and existing ancestry fitting/development roles.
These are candidates, not certified unrelated or biologically pure donors.
Do not label every geographically Oceanian population as OCE or every
African population as SSA.

The publisher's [2021 release notes](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/Readme)
state that GT/DS/GP are reference-imputed. Those fields cannot serve as
observed benchmark genotypes. Its PL/AD fields instead retain sequencing
likelihoods/read depths, confirmed in the BCF header and a real indexed
query. The [import protocol](./config/sgdp_public_import.json) re-calls
unphased diploid GT from PL/AD only, preserving missing calls. It uses
predeclared depth/likelihood-gap and callability thresholds and checks
concordance against the original high-quality hs37d5 Chane variant calls.
No absent variant is silently assumed homozygous reference.

The [Nextflow smoke](./public_sgdp_import.nf), launched using
`bash run_public_sgdp_import.sh --smoke`, imports chr22 likelihoods for these
16 fully public candidates, verifies publisher MD5s and retrieves only the
required original public tar member using checked byte ranges and its
publisher MD5. Signed-letter/DO_NOT_USE samples are excluded.

The smoke completed successfully (exit 0). The independent verifier
reconstructed all 16,887,264 genotype cells across 1,055,454 chr22 SNPs.
Concordance with original high-quality Chane calls was 99.92117% over
16,492 shared calls, exceeding the predeclared 99.5% threshold.
Its exact-session disposable cache was cleaned from approximately 1.2 GiB
to 4 KiB, then all published genotype evidence was independently verified
again. Reusable calls, metadata, raw-control evidence and provenance remain.
See the [verification report](./results_v3_public_sgdp_smoke_20261021_203540/published_verification.json).

Retained smoke controller:

```bash
docker logs -f snpster-ancestry-v2-sgdp-20261021_203540
```

The [all-autosome import](./public_sgdp_full_v2.nf), launched with
`bash run_public_sgdp_full_v2.sh --run`, reuses verified chr22 and
downloads/re-calls/QCs chromosomes 1-21 with the same frozen protocol.
The full import completed successfully (exit 0) on 2026-10-21 at 21:30
local time. The initial full-import attempt failed before scientific tasks
because chromosome-dependent publish paths were evaluated too early;
its controller, sources and failure evidence remain unchanged.
The separate recovery uses task-scoped publish-path closures and retains
per-chromosome download provenance and collision-safe QC staging.

Retained successful recovery controller:

```bash
docker logs -f snpster-ancestry-v2-sgdp-full-v2-20261021_204615
```

Its final manifest represents producer QC only. Independent all-autosome
genotype verification must precede its own exact-run cache cleanup.
Producer summaries report 77,818,346 biallelic SNP records for the 16
candidates (1,245,093,536 genotype cells including missing calls), with
1,885,342 high-quality Chane raw-call comparison anchors. All chromosome
concordances exceed the predeclared 99.5% threshold; the lowest is 99.8826%.
These are genotype-import checks, not ancestry-accuracy measurements.

The [independent all-autosome verification workflow](./public_sgdp_verify.nf)
is launched with `bash run_public_sgdp_verify.sh --run`. It reconstructs
every GT/GQ/DP from saved PL/AD, checks the locus/sample axes and functional
indexes, and binds all 22 reports to published genotype, metadata and
frozen producer hashes. Its [acceptance tests](./tests/test_public_sgdp_verification.py)
exercise non-chr22 input, separate shared metadata, wrong-chromosome and
altered-genotype rejection, index disagreement, complete coverage and
post-cleanup hash failures. All 31 focused public-import/preflight/contract
tests pass.

Retained failed verification controller:

```bash
docker logs -f snpster-ancestry-v2-sgdp-verify-20261022_080336
```

The verifier exited 1 on 2026-10-22 at 08:10 local time, correctly rejecting
duplicate exact allele `8:89329148:G:A`. Inspection of both published calls
and the original publisher BCF found one sequencing-supported record and
one record missing PL/AD for all 16 selected people. The initial importer
preserved both instead of resolving the redundant source record.
Chromosomes 1, 2, 3, 4, 5 and 7 had independently passed; all original
sources, reports and failure evidence remain unchanged. No cache cleanup
occurred. HTML/timeline rendering warnings during failed shutdown are not
evidence of successful execution reports.

The separate [recovery workflow](./public_sgdp_recovery.nf), launched with
`bash run_public_sgdp_recovery.sh --run`, rechecks hashes before reusing those
six chromosome reports. It audits every unverified chromosome for duplicate
exact keys and creates versioned corrected inputs only where needed.
It keeps one unchanged sequencing-supported record when other copies are
identical or wholly missing, or one copy of identical missing records.
Conflicting raw evidence fails explicitly; records are never merged and
no genotype is filled or recalculated under a different threshold.
Removed and retained records are saved as duplicate evidence. Corrected
inputs receive fresh raw-control QC, indexes, callability counts and the
same strict independent genotype/uniqueness checks.

Retained successful recovery controller:

```bash
docker logs -f snpster-ancestry-v2-sgdp-recovery-20261022_082600
```

All 36 focused tests pass, including real compressed-VCF/index repair,
retained genotype equality, reversed duplicate order, conflicting raw
evidence rejection and durable hash binding through Nextflow-style staging.
The recovery manifest exposes per-chromosome data roots so later workflows
cannot inadvertently reuse the uncorrected duplicate input.

The recovery completed successfully (exit 0) on 2026-10-22 at 08:35 local
time. All 22 autosomes passed independent verification, reusing the six
previously verified chromosome reports after hash checks. Five redundant
records were removed: one on chr8, two on chr12, one on chr14 and one on
chr17. No conflicting duplicate evidence was accepted.

The corrected dataset contains 77,818,341 SNP records and 1,245,093,456
independently reconstructed genotype cells across 16 people.
All 1,885,342 high-quality Chane comparison anchors remain; the minimum
chromosome concordance is 99.882557%. This is sequencing-call concordance,
not ancestry prediction accuracy.

The exact completed import session
`2d2f3b41-7812-4261-aa63-91fe306fe82f` was cleaned using a checked Nextflow
dry-run plan, reducing its approximately 74 GiB task work to 4 KiB. All
222 durable evidence hashes passed again afterward. Following a fresh
pre-cleanup hash check, the recovery session
`9eb94212-edeb-4ef6-bea3-5265e043d221` was also cleaned from approximately
546 MiB to 4 KiB; all 222 evidence hashes passed once more.
Original published calls, corrected calls, duplicate evidence, sources,
execution logs and failed-run evidence remain retained.
See the [all-autosome verification](./results_v3_public_sgdp_recovery_20261022_082600/summary/published_verification.json)
and [final post-cleanup check](./results_v3_public_sgdp_recovery_20261022_082600/summary/post_recovery_cleanup_verification.json).
Coordinate/allele harmonization and genotype identity/relatedness remain
explicitly incomplete regardless of the genotype-verification outcome.

This is data-quality preparation, not the ancestry sweep. Full-autosome
data, coordinate/allele harmonization, genotype identity/relatedness,
operational source-label review and the fixed cohort/candidate manifests
still gate scientific launch. Both separate and merged EAS/SEA ancestry
validation remain mandatory. The [independent published-data verifier](./pedigree_benchmark/scripts/verify_public_sgdp.py)
rederives every genotype from saved sequencing evidence before exact-run
cache cleanup is authorized.

### Autonomous public-data normalization and genetic eligibility

Routine public-data normalization is authorized without further user
confirmation. The [normalization protocol](./config/sgdp_normalization.json)
and [Nextflow workflow](./public_sgdp_normalization_v3.nf) bind the verified
all-autosome recovery manifest, original hs37d5 FASTA, the production
GRCh38 FASTA, and the existing hg19-to-hg38 chain by content hash.
Production phasing/imputation and all frozen historical protocols remain
unchanged. New sequencing inputs are normalized to verified GRCh38 before
constructing GRCh38 microarray inputs; this does not change the production
standardizer's handling of user uploads.

VCF positions are one-based and UCSC chain positions are zero-based:
map `POS - 1`, then emit `mapped_position + 1`. Accept one mapping on the
original autosome; record unmapped, ambiguous, different-contig and
reference-incompatible exclusions explicitly. Verify REF against the
source FASTA, complement both alleles on reverse mappings, and align REF
to the production GRCh38 FASTA. If REF/ALT must swap, change allele indices
and reverse the PL/AD axes together, preserving nucleotide genotypes,
quality/depth and missingness. Reject all exact target-key collisions.
Do not infer phase, fill missing calls or substitute imputed genotypes.

The independent verifier parses UCSC chain blocks separately from
`pyliftover`, restores normalized records to their original source-row
order, and compares every retained diploid nucleotide genotype and raw
likelihood/depth axis with its original sequencing evidence. It separately
checks mapping exclusions, sample order, sorted unique loci, reference
bases and functional index counts.

Genetic eligibility uses exact shared allele keys on all 22 autosomes
across normalized public calls, the canonical ancestry panel and the
unchanged production reference. The comparison universe is the union of
all production-reference and existing ancestry fitting/development
identities, plus the 16 candidates. No absent genotype is filled.
Predeclared filters are MAF >=0.05, marker missingness <=0.05 and
500 kb / r-squared 0.2 LD pruning. Require at least 20,000 pruned markers,
10,000 jointly called markers in every candidate-involving pair, and
90% sample callability. KING kinship >=0.0442 or genotype concordance
>=99.5% flags candidate/reference relationships; candidate/candidate
relationships become split constraints. These are conservative operational
screens, not source-purity or absolute-unrelatedness certificates.

The task image is pinned to
`sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b`,
built from the existing scientific image with pinned `pyliftover==0.4`
and `zstandard==0.23.0`. The latter was absent in the original image.
All 44 focused tests pass, including wrapped FASTA boundaries, forward
and reverse chain endpoints, one-based coordinates, REF/ALT likelihood
reordering, real VCF/index normalization, independent genotype validation
and actual PLINK KING output parsing. The new integration checks run the
all-22 synthetic concat/filter/LD-pruning/KING path, detect a candidate
duplicate, verify complete candidate-pair coverage, and reject mutated
durable collector evidence.

The real chr22 smoke passed: 1,053,778 of 1,055,454 source sites retained,
16,860,448 genotype cells independently verified, and 1,876 exact shared
identity markers across 4,627 people. Explicit exclusions were 899
different-contig mappings, 773 unmapped sites, two target-reference
mismatches and two target-collision records. A staged-input/output folder
collision initially blocked the identity-panel task; a separate v2
recovery staged its input as `normalized_input`, reusing all three
completed tasks through Nextflow `-resume`. Its verified task cache was
cleaned to 4 KiB and published evidence rechecked.

The first full attempt exposed a separate scope bug: CLI `--smoke false`
was truthy in a conditional and ran only chr22. It correctly failed the
all-autosome completion gate, without reporting successful full results.
The v3 workflow validates the parameter as literal `true` or `false` and
uses the resulting boolean consistently. A real Nextflow preview and
the completed controller explicitly report chromosomes 1 through 22.
Its [recovery launcher](./run_public_sgdp_normalization_v3.sh) resumes the
same session, reusing assembly preparation and all three chr22 tasks;
original failure logs and producer hashes remain preserved.

Completed all-autosome controller (exit 0; logs retained):

```bash
docker logs -f snpster-sgdp-normalization-v3-full-20261022_100116
```

The complete dependency chain finished without another confirmation:
normalization, independent genotype verification, shared-marker panels,
all-candidate KING comparisons, durable-manifest verification, exact
successful-session cache cleanup and post-cleanup verification.
There is no ancestry-model sweep running.

The [completed normalization manifest](./results_v3_sgdp_normalization_full_v2_20261022_095828/summary/normalization_complete.json)
contains 77,768,904 retained SNP records and 1,244,302,464 independently
verified genotype cells across all 22 autosomes and 16 people.
The [genetic audit](./results_v3_sgdp_normalization_full_v2_20261022_095828/kinship/genetic_eligibility.json)
used 192,903 exact shared markers, 66,420 LD-pruned markers, 4,627 people
and all 73,896 candidate-involving pairs; the minimum pair had 65,632
jointly called markers.

Thirteen candidates passed this screen: 11 NAM and two OCE. All three
new Dinka candidates were excluded for strong genetic matches to existing
comparison samples under different identifiers:

| Public candidate | Existing identity | KING kinship |
| --- | --- | --- |
| B_Dinka-3.DG | GRC10046841 | 0.438630 |
| S_Dinka-1.DG | GRC10046839 | 0.437544 |
| S_Dinka-2.DG | GRC10046845 | 0.434369 |

These high-relatedness matches are not evidence for independent SSA
donors, even though metadata aliases did not overlap. Do not use those
three candidates as independent benchmark founders. No candidate-to-
candidate pairs crossed the relationship thresholds. Passing this
screen does not prove biological source purity or absolute unrelatedness.
The dedicated work directory was reduced to 4 KiB; the
[post-cleanup verification](./results_v3_sgdp_normalization_full_v2_20261022_095828/summary/post_cleanup_verification.json)
rechecked all 175 retained evidence-file hashes.

The launcher verifies the real published outputs before exact-session
cache cleanup and checks them again afterward. It retains normalized
inputs, assembly resources, exclusion/provenance reports and genetic
identity evidence. Failures retain their work, sources and logs.
Completing normalization and relatedness screening does not authorize
ancestry scoring: source-label review, family/role separation, frozen
cohorts/candidates and the production-matched preprocessing smoke remain.

## Historical benchmark contract: flexible pipeline placement (2026-10-18)

The user clarified that ancestry can run at whichever pipeline stage
performs best: raw harmonized array GTs, phase-only output, or phased
and imputed output. The [shared benchmark contract](./config/ancestry_benchmark_contract.json)
supersedes a fixed pre-imputation placement requirement. No production
module is moved until matched accuracy/integrity evidence supports it.

FLARE and Gnomix require phased haplotypes; unphased raw GTs are only a
direct input for compatible phase-independent methods. Arbitrary GT
allele order must not be treated as phase. Original inherited phase
remains a mechanistic diagnostic, not a deployable pipeline step.
Processing tracks must use matched people/truth/labels/donor exclusions,
with stage-native panel coverage and shared-panel ablations reported.
Methods within each stage share exact allele panels. Native and learned
calibration scores remain separate, and selection uses donor-purged
training evidence rather than a real individual's expected ancestry.

The completed observed-only pilot below remains valid for that specific
track. Its frozen inputs and no-imputation protocol are unchanged.
Prior partial imputed jobs remain stopped; permitting an imputed track
does not make their different/incomplete cohorts a fair stage comparison.
Any resumed/new comparison requires an explicit compatible version.

### Training-only FLARE mechanism diagnostic

Before a broad stage or hyperparameter sweep, the fixed
[native FLARE diagnostic](./config/flare_native_diagnostic.json) isolates
copying/mismatch initialization on seven hash-selected pure-training
controls, one per source, disjoint from ancestry reference and evaluation
people. It uses chr22, both masks, and paired inherited/statistical phase
with the same observed GTs. This is a mechanism check, not a seven-person
whole-genome leaderboard or new independent certification.

Five contrasts hold the exact panels, reference, seed, physical
aggregation and rho parameters constant:

1. Native reference-estimated copying, fixed parameters.
2. The same initial model with singleton EM.
3. 90% shrinkage of the copying matrix toward named source-matching panels.
4. 100% source-matching copying.
5. 90% source anchoring with mismatch probability 0.001.

Explicit fixed-model files retain T=10 and uniform seven-source priors.
The source-to-panel identity matrix is constructed by names, not array
position. No target ancestry truth is supplied to inference. These are
controlled hypotheses, not an assumption that cross-source copying is a
bug or that diagonal copying is scientifically preferable.

The [producer](./pedigree_benchmark/scripts/diagnose_flare_native.py)
audits all 3,033 reference labels against authoritative source/fine
metadata, verifies reference/training/evaluation donor separation and
compares stock FLARE with the precision-export build. The jars differ
only in the global-export class; runtime checks compare SAS training
control local-posterior fractions under both phases. Every diagnostic
prediction is also independently reaggregated from its local posteriors,
and phase-only GT preservation is checked with a separate VCF reader.
Three focused tests cover named/permuted copying axes, disjoint
hash-selection and commands without target truth priors.

The [runner](./run_flare_native_diagnostic.sh) uses a fresh pinned,
network-disabled container with 12 CPUs/32GiB and two mask workers.
Use `--prepare`, then `--run`; default output is
`results_v3_flare_native_01`. Sources and inputs are SHA-bound.
No frozen benchmark/model is changed and no production/PGS promotion
is part of this experiment. A subsequent stage comparison will depend
on what this FLARE diagnostic actually establishes.

The diagnostic completed successfully on 2026-10-18 in approximately
five minutes. It produced 140 method/phase/mask predictions from seven
unique training people, not 140 independent controls. All 140 local
posterior aggregations passed independent recomputation, all phase-only
diploid GTs were preserved, and all four stock-versus-precision SAS
comparisons agreed exactly in physically aggregated local fractions.
The full reference label audit found no label mismatch or donor overlap.

Statistically phased chr22 training-control mean TV:

| FLARE contrast | Lonni mask | Frederik mask |
|---|---:|---:|
| Native fixed | 8.37% | 5.75% |
| Native singleton EM | 5.96% | 6.35% |
| 90% source anchoring | 16.38% | 12.55% |
| 100% source anchoring | 24.04% | 19.37% |
| 90% anchoring, theta=0.001 | 15.89% | 12.41% |

Source anchoring worsened these controls under inherited phase as well.
The results do not support treating native cross-source copying as a
simple bug or replacing it with a diagonal identity model. They do not
isolate a native-model mismatch-probability effect, since the theta
perturbation was tested only together with 90% anchoring.

The GME control HGDP00741 remains the worst native case: native fixed
statistical-phase TV is 41.96% / 32.96% on Lonni / Frederik, versus
34.05% / 31.76% with inherited phase. The selected SAS training control,
HG04080, instead has native fixed statistical-phase TV of 4.58% / 3.44%.
Thus the previous difficult SAS individual cannot be generalized to
every SAS donor; ancestry-specific means can also hide a large GME tail.
These are operational source-label errors, not evidence that an
individual is biologically ancestry-pure.

Results and verification manifests are saved in
[the diagnostic output](./results_v3_flare_native_01/), including
[comparison](./results_v3_flare_native_01/comparison.tsv),
[predictions](./results_v3_flare_native_01/predictions.tsv) and per-mask
verification. No candidate is promoted from this seven-person, one-
chromosome diagnostic; these means must not be substituted for the
previous all-autosome pedigree scores. The next substantive comparison
needs broader donor/population coverage and matched marker-density /
reference-model contrasts, with pedigree controls and locked
all-autosome follow-up. A pipeline stage has not yet been selected.

### Broader training-only reference/density screen (2026-10-18)

The next [locked protocol](./config/flare_reference_density.json) expands
from seven controls to all 240 approved pure-training people and all 76
descendants from 25 complete training families. It uses no evaluation
person, evaluation founder or new external data. Source references and
all training founder identities remain disjoint. Related descendants,
reused donors, chromosome blocks, panels and phase conditions are not
independent extra people.

This is a three-chromosome screen (22, 12, 1), **not a whole-genome
accuracy result or independent certification**. Chromosomes are reported
separately. No final candidate or pipeline stage is chosen automatically.

The native FLARE contrasts are:

| Inference reference | People | Purpose |
|---|---:|---|
| Full existing reference | 3,033 | Baseline |
| At most 16 per fine population | 905 | Fine-population representation contrast |
| Broad-source-count-matched hash sample | 905 | Control for reference size and source counts |

Both 905-person panels have exactly the same counts in every broad
source, including all 24 OCE and 49 NAM people. Capping uses authoritative
fine labels; selection is deterministic and does not use prediction
errors. The single subset seed is exploratory; promising differences
will need replication, not a causal claim from one reference draw.

Each reference is tested with fixed parameters and singleton EM,
`update-p=false`. Native reference-initialized copying/rho/theta are
estimated afresh for each exact panel/reference/chromosome. Explicit
fixed models retain T=10 and uniform seven-source priors; no ancestry
truth, diagonal anchoring or learned calibrator is supplied. A direct
implicit-native versus serialized-model check verifies that model
injection does not change the initialization control.

Each chromosome has three exact marker panels:

- The retained Lonni mask.
- The retained Frederik mask.
- A spatially thinned subset of Frederik, retaining exactly the Lonni
  marker count and chromosome endpoints.

The thinned panel is a density ablation, not another real array or
imputed data. It separates density from differences between the two
actual arrays' SNP choices, though thinning also changes spacing/content.
Parent full-reference marker filters and physical truth spans are held
fixed, rather than re-filtering per inference-reference candidate.

Statistical phasing is singleton and uses the full reference on that
exact marker panel. The same phased GTs are then supplied to all inference
reference contrasts, so the reference comparison does not also change
the phasing algorithm's reference. Original inherited phase is a paired
mechanistic control. Raw-GT compatible methods and phased+imputed tracks
remain permissible future pipeline placements; this screen does not
test imputation or decide its placement.

The [producer](./pedigree_benchmark/scripts/screen_flare_reference_density.py)
validates every phase-only GT with a separate VCF reader and independently
reaggregates every native local posterior. Chromosome-specific pedigree
truth is reconstructed from complete half-open ancestry segments over
both copies, not global pedigree fractions or retained-SNP counts.
Per-call verified output hashes preserve useful completed evidence if
a JVM call fails; an incomplete task is not silently overwritten/retried.

Output includes per-person predictions, mean/p95/max TV, maximum
component error, all-components-within-5/10pp coverage, zero-truth
off-target mass, family-balanced pedigree TV, source-balanced pure TV,
fine-population diagnostics and inference time. Family/source balancing
is descriptive, not donor-purged supervised CV or an independent
confidence interval. Existing donor-purged calibration protocols remain
unchanged.

Preflight passed eight focused tests and a real 72-call smoke on two
hash-selected training people, all three panels, six configurations and
both phases. The [separate verifier](./pedigree_benchmark/scripts/verify_flare_reference_density.py)
also recomputed all 228 training descendant/chromosome truth vectors
directly from original Ped-sim BP files and plans, checked 33,306 supplied
diploid GTs against the original people/haplotype phase, independently
reaggregated all 72 smoke local posteriors and recomputed persisted
metrics/coverage and hashes. A Python 3.10 compatibility failure in the
verifier's hash helper was fixed and regression-tested before launch;
the frozen scientific producer was not changed.

The [runner](./run_flare_reference_density.sh) supports `--prepare`,
`--smoke`, then `--run` and refuses an existing run container. Default
output is `results_v3_flare_reference_density_01`. The full screen
contains 34,128 method predictions from 316 unique training people,
across nine panel/chromosome blocks. It runs in a pinned network-disabled
container with 12 CPUs, 32GiB and at most two six-thread workers. Prior
models, results and pipeline/PGS code remain untouched.

After the full run, the separate verifier can verify its completion,
original GT concordance, BP truth, summary tables and additional
local-posterior spot checks without `--smoke`. It writes an independent
verification manifest; producer completion alone is not that manifest.
Genuine Gnomix baselines are retained. Any subsequent locked
all-autosome finalist comparison must refit Gnomix on the same exact
new panels rather than reusing an incompatible earlier model.

#### Audited results and incomplete-run status (2026-10-21)

The run exited **1**, not successfully, on 2026-10-19 at 09:07 UTC,
after approximately 21.5 hours. Eight of nine blocks completed. The
dense Frederik chr1 block completed full-reference fixed/EM but stopped
at the first fine-capped fixed initialization control. Four configurations
there are missing: fine-capped fixed/EM and size-control fixed/EM.
No top-level completion manifest was produced; the original failed
state and outputs remain unchanged.

The failure was the implicit-native versus explicit-model initialization
assertion, not the independent local-posterior aggregation check.
Maximum fraction difference was 5.26278349e-8 (approximately 0.0000053
percentage points), above the prescribed 1e-10 tolerance. The compared
model files have identical copying, rho, theta and generation values;
the native exported uniform prior is 0.14285715 versus the writer's
double-precision 1/7. This identifies a representation difference but
does not establish the exact numerical cause of the posterior mismatch.
Do not silently relax the frozen check or report the sweep as complete.

A separate [partial-results audit](./results_v3_flare_reference_density_01/report_20261021/verification.json)
verified all 30,336 predictions in eight completed blocks and the 1,264
checkpointed full-reference predictions in the incomplete block.
It rechecked frozen inputs/producers and completed-table hashes,
reconstructed all 228 descendant/chromosome truth vectors directly from
raw Ped-sim BP/plans, verified 50,121,392 typed diploid GTs against the
original person/phase and phase-only output, and independently aggregated
200 local-posterior spot checks across every completed configuration.
Every saved TV, maximum component error and off-target fraction was
recomputed. The [audited comparison](./results_v3_flare_reference_density_01/report_20261021/comparison.tsv)
explicitly flags incomplete blocks; no partial-autosome global score
was substituted for the planned result.

For singleton EM on actual statistically phased arrays, mean TV is:

| Chromosome | Mask | Full reference: pure | Fine-capped: pure | Full reference: pedigree | Fine-capped: pedigree |
|---|---|---:|---:|---:|---:|
| 12 | Lonni | 17.33% | 6.51% | 7.53% | 7.45% |
| 12 | Frederik | 19.53% | 7.86% | 7.36% | 7.11% |
| 22 | Lonni | 14.14% | 7.25% | 11.64% | 11.34% |
| 22 | Frederik | 18.75% | 9.25% | 10.44% | 10.12% |

These are within-screen paired comparisons on the same training people,
not improvements relative to the old 14-person evaluation pilot.
The smaller reference substantially improves pure-source mean TV, with
much smaller pedigree changes. It is not uniformly superior: on chr1
Lonni, fine-capped EM improves pure TV from 13.65% to 6.50% but worsens
pedigree TV from 8.15% to 9.35%.

The source-count-matched 905-person control also helps. On chr12
Frederik, it has pure/pedigree TV 9.13%/6.88%, compared with
7.86%/7.11% for the fine cap. Thus the evidence does not attribute all
improvement to fine-population balancing. Reference size/composition and
its newly estimated native parameters change together, including the
default theta; a single deterministic subset draw is not a causal
isolation of any one mechanism.

Source tails still fail important guardrails. On chr12 Frederik,
fine-capped EM GME mean TV is 40.24% over 12 controls (none within 10pp),
versus 29.41% with full-reference EM. SAS is 11.13%. Conversely,
fine-capping lowers SSA mean TV from 25.17% to 4.66%, EUR from 22.41%
to 5.34%, and EAS/SEA from 12.34% to 3.64%. Full-reference EM assigns
the SSA controls an average 24.5% NAM, versus 4.2% for the fine-capped
reference. These are errors against operational source labels, not proof
that each control has biologically pure ancient ancestry.

Even the fine-capped chr12 Frederik result has only 71.1% of descendants
within 10pp on every component and 55.3% within 5pp; pedigree p95 TV is
25.74%. Pure-source mean TV is 7.86% but source-balanced TV is 10.57%,
highlighting the GME tail and unequal source sample sizes.

Phasing remains relevant: fine-capped chr12 Frederik pedigree TV is
3.87% with inherited phase versus 7.11% with statistical phase, despite
unchanged GTs. More markers are not an unconditional fix. At the same
905-person fine-capped reference, dense versus thinned Frederik chr12
has pure TV 7.86% versus 5.77% and pedigree TV 7.11% versus 7.67%;
on chr22 it has pedigree TV 10.12% versus 12.51%.

No model or pipeline stage is promoted. The next version should use
Nextflow to reuse verified phase/reference inputs and explicitly recover
the four missing dense-chr1 configurations, after characterizing the
initialization discrepancy. Further training-only controls should isolate
reference size/default mismatch parameters from population composition.
Any locked finalists still require all-autosome development validation
and an exact-panel genuine Gnomix comparison. The current screen does
not benchmark imputation or provide independent certification.

### Nextflow recovery and crossed-mismatch follow-up (2026-10-21)

The [follow-up protocol](./config/flare_nextflow_followup.json) preserves
the failed screen and its exact scientific parameters. Its new
[Nextflow DSL2 workflow](./flare_followup.nf) has separate cacheable
preflight, initialization characterization, singleton inference and
aggregation/verification processes. It uses the existing pinned image,
up to 12 CPUs/30GiB of scheduled task memory and at most two six-thread,
12GiB inference tasks. Python implements tasks, not parallel scheduling.

Two deliverables are kept separate:

1. Recover only the 2,528 missing dense Frederik chr1 predictions. Reuse
   the 31,600 audited predictions, all references and all phase-only
   inputs. The resulting recovered screen must have exactly 34,128
   distinct panel/chromosome/method/phase/person combinations.
2. Run 3,792 new singleton-EM crossed-theta controls on chr12, both
   actual masks, both phase conditions and all 316 training people.
   The full-reference model takes the native 905-person theta; each
   905-person model takes the native full-reference theta. Preserve its
   own copying, rho, T=10 and uniform priors, mapping matrix axes by
   names. Reuse 3,792 matched native-theta baselines.

The second experiment isolates a theta intervention within each native
reference model. It does not independently isolate every aspect of
reference composition/initialization or constitute a sweep over many
theta values. No evaluation families or real-person expected answers
enter selection. No whole-genome certification, automatic winner,
production move or PGS change is part of this follow-up.

Before inference, the actual failed initialization is replayed four
ways: implicit native initialization, two identical explicit fixed-model
runs, and the native exported model as explicit input. In the real
smoke, the explicit repeats agree exactly; implicit repeat and exported
rounded-prior injection both reproduce the original implicit fractions
exactly. Double-precision uniform-prior injection reproduces the
5.26278349e-8 fraction difference. This ties the replayed discrepancy
to parameter representation, not changed genotype data. The original
failed producer/check remains unchanged.

The new characterization gate explicitly requires repeated explicit
predictions within 1e-10 and all aggregate initialization differences
within 1e-6. These gates are declared in the new protocol rather than
retroactively changing the old failure. Independent per-call
local-posterior aggregation remains at 1e-10.

Eleven focused tests and a real 40-call smoke pass. Preflight rechecks
frozen input/source hashes, training-only identities, original person
GT/phase concordance, independent phase-only GT invariants and raw BP
truth. All new local posteriors are independently aggregated at inference
and again at aggregation. A separate
[published-output verifier](./pedigree_benchmark/scripts/finalize_flare_nextflow.py)
then checks durable per-call models/local VCF hashes, metadata, metrics,
summary tables and exact full-run Cartesian coverage before cleanup.

Use the [version-2 launcher](./run_flare_followup_nextflow_v2.sh):

```bash
bash run_flare_followup_nextflow_v2.sh --prepare
bash run_flare_followup_nextflow_v2.sh --smoke
bash run_flare_followup_nextflow_v2.sh --run
# For an interrupted full run:
bash run_flare_followup_nextflow_v2.sh --resume
```

Default output is
[results_v3_flare_nextflow_02](./results_v3_flare_nextflow_02/).
The initial launcher is retained as a SHA-bound implementation dependency;
do not use it directly for new launches.

#### Verified cache cleanup

The user requested removal of analysis cache after completion.
The launcher waits for successful Nextflow completion **and** verified
durable outputs before running `nextflow clean` for the exact session
UUID inside that mode's dedicated launch directory. It checks the dry-run
paths against the dedicated work directory and confirms task files are
gone before marking cleanup complete. Failure before verification keeps
the work/cache for explicit resume.

The initial cleanup smoke exposed that a new display name on `-resume`
did not identify the original cached session: a name-based clean returned
without releasing task files. The version-2 wrapper was prepared in a
new output version rather than editing the frozen launcher. Its real
resume test reused all 43 tasks, then UUID-scoped cleanup reduced the
work directory from approximately 40MiB to 156KiB of empty directory
entries, while retaining all 16MiB of published smoke call evidence.
The first smoke's own cache was also removed by its exact session UUID.

Published predictions, compact models/local posteriors, manifests,
characterization, execution traces, timeline and resource reports are
retained. Shared scientific inputs, original failed-run evidence,
global Nextflow assets and other runs are never pruned. Cache cleanup
means a verified completed run no longer has its disposable task cache
available for `-resume`; a new scientific experiment uses a new output
version and can reuse the persisted verified inputs.

#### Completed Nextflow follow-up (2026-10-21)

The follow-up completed successfully at approximately 18:48 local time.
All 6,320 inference tasks passed, restoring the exact 34,128-row
three-chromosome screen and producing 3,792 new crossed-theta
predictions alongside 3,792 matched native baselines.
The [durable-output verification](./results_v3_flare_nextflow_02/full/published_verification.json)
passed before and again after disposable cache removal. It verifies
published per-call local-VCF/model hashes, every prediction key, source
metadata, saved error metrics, summary tables and Cartesian coverage.
The full preflight independently rechecked 35,591,396 reused typed GTs
and all 228 descendant/chromosome raw-BP truth vectors. Each new local
posterior was independently reaggregated both at inference and at
aggregation.

There was one explicit same-parameter recovery. At approximately
13:56 local time, Temurin 21.0.12.1 crashed in compiled
`java.io.BufferedWriter.write` during
`recovery_frederik_1_finecap16_em_array_phase_only_NA19435`.
The native SIGSEGV surfaced as SIGABRT. This was not a scientific
assertion or established OOM event. Crash reports, commands, partial
native outputs and workflow/trace logs were copied outside the work
cache before `--resume`. The successful resume reused 1,778 inference
tasks and both preflight tasks, reran the failed/aborted work with the
same image, JVM command, model, input, seed and threads, and completed
the remaining tasks. Original failure evidence remains under
[failure_evidence](./results_v3_flare_nextflow_02/full/failure_evidence/).

**Changing the native mismatch probability is not the main explanation
for the reference-size effect in these controls.** For chr12 Frederik,
singleton-EM statistical-phase mean TV is:

| Reference | Pure, native theta | Pure, crossed theta | Pedigree, native theta | Pedigree, crossed theta |
|---|---:|---:|---:|---:|
| Full 3,033 | 19.53% | 19.57% | 7.36% | 7.29% |
| Fine-capped 905 | 7.86% | 7.89% | 7.11% | 7.08% |
| Size-matched 905 | 9.13% | 9.16% | 6.88% | 6.92% |

Native theta is approximately 8.949088e-6 for the full reference and
3.452333e-5 for the smaller panels. Swapping these while holding the
rest of each model fixed changes the primary pure/pedigree group mean
TV by less than 0.1 percentage point for both masks. It does not
reproduce or remove the large pure-source reference effect.
This does not exclude effects of other theta values or other
initialization parameters.

Individual sensitivity is not zero: under full-reference EM on
statistically phased Frederik data, two of 316 people move by more
than 5pp in at least one ancestry component, with maximum component
shift 13.14pp. Stable group means must not be interpreted as identical
per-person predictions or universal numerical robustness.

GME remains problematic. On chr12 Frederik, pure GME TV is
29.41% / 30.55% for full-reference native/crossed theta and
40.24% / 40.33% for the fine-capped reference. The theta intervention
does not solve the source-tail problem. The inherited/statistical-phase
pedigree gap also persists: fine-capped Frederik native theta has
3.87% / 7.11% TV, versus 3.86% / 7.08% after crossing theta.

The recovered dense chr1 results confirm the reference trade-off:
full-reference versus fine-capped singleton EM gives pure TV
17.25% versus 8.61%, but pedigree TV 8.46% versus 9.20%.
There is no uniform reference winner across pure/pedigree tracks.

Results are in the
[recovered-screen comparison](./results_v3_flare_nextflow_02/full/summary/summary/recovered_screen_comparison.tsv),
[crossed-theta comparison](./results_v3_flare_nextflow_02/full/summary/summary/crossed_theta_comparison.tsv)
and their per-source/per-person companion tables.
No automatic candidate or pipeline placement is promoted. These remain
training-only, three-chromosome diagnostics, not all-autosome or
independent certification. Subsequent experiments should focus on
reference structure/native copying-rho estimation, GME separation and
phase robustness rather than another theta-only sweep.

Exact-session Nextflow cleanup completed only after durable-output
verification. Remaining empty work directories were removed separately;
the full work root is now approximately 4KiB. Approximately 2.4GiB of
published call evidence and 3.9MiB of crash evidence are retained as
results, not disposable cache. No other run, shared reference or
global Nextflow dependency was pruned.

## Previous observed-only benchmark contract (2026-10-17)

**Historical scope of the completed observed-only pilot, not a current
pipeline-placement restriction.** The user initially requested ancestry
before imputation. Phase-only
processing of typed SNPs is permitted; adding untyped variants, filling
missing GTs or changing supplied diploid GTs is not. Imputed-array results
remain historical diagnostics, not deployment acceptance gates.

The owned shared imputation scheduler and three waiting evaluators were
stopped with user approval. Completed models, verified CV/development
results, chromosome outputs and immutable old protocols are preserved.
Old state files may still say `waiting_for_shared_array_verification`;
they describe the superseded experiment, not an active requirement.

The initial [common-panel protocol](./config/observed_ancestry_protocol.json)
and [compatibility runner](./run_observed_ancestry_benchmark.sh) established
native FLARE and genuine pinned Gnomix LR+CRF interfaces. Its four variants
are FLARE with `array=true` and fixed native parameters, FLARE with
`array=true` and singleton EM, an `array=false` fixed-parameter ablation,
and Gnomix with a predeclared2-cM target window for sparse input.
No inherited GMM decoder or old dense-panel Gnomix model is substituted.

Standardization requires:

- Identical exact CHROM/POS/REF/ALT panels, reference identities, source
  labels, evaluation people and founder exclusions across methods.
- Singleton phase-only processing using `impute=false` and a reference
  restricted to those same typed sites. Every typed diploid GT and every
  marker key must survive unchanged; missing GTs fail explicitly.
- Separate primary statistically phased and diagnostic inherited-phase
  conditions, with identical GTs/markers in both.
- Physical-midpoint posterior aggregation over the original full
  half-open chromosome spans, paired copies equally weighted, and the
  same physical-span ancestry truth for every method.
- Separate mask, pure/pedigree and native/calibrated reporting. The
  initial pilot has no learned calibration or winner selection.
  Report TV, tails, all-component5/10-point coverage, per-source errors,
  off-target mass, phase diagnostics and native runtimes; never score
  partial-autosome results.
- Refit Gnomix on each observed-marker panel. Do not pad missing sites
  or reuse a dense-marker model with a different input axis.

Four focused invariant tests pass. A real chr22 smoke completed all four
methods under both masks on one pure and one pedigree control, using
both phase conditions. An [independent bcftools reread](./results_v3_observed_benchmark_02/independent_smoke_verification.json)
verified790 typed diploid GTs and zero added/untyped SNPs. Beagle omits
contig metadata in this phase-only output; explicit metadata-only
normalization precedes exact key/GT validation. The initial failed
metadata smoke is retained under `results_v3_observed_benchmark_01/`.

**Do not treat this as a full-array accuracy benchmark.** The frozen
197,260-site common panel retains only6,205/148,019 Lonni typed SNPs
(4.19%) and27,783/839,144 Frederik typed SNPs (3.31%). The
[coverage report](./results_v3_observed_benchmark_02/coverage.tsv) and
prepared compatibility pilot are persisted, but a full run is held
out of that restricted panel rather than presented as full-array accuracy.

Existing dense data are a viable alternative without downloading people:
chr22 contains every one of the2,245/12,355 observed mask keys. However,
only3,033/3,447 current reference identities and240/307 pure-training
identities occur in that resource. All seven source groups remain
represented, but OCE pure-training support falls from six people to four.
Using that resource therefore requires an explicitly versioned shared
reference/training subset and regenerated compatible training inputs,
not imputation or silent reuse of incompatible old models.

The compatibility runner defaults to `results_v3_observed_benchmark_02`.
For a distinct experiment choose a fresh `OBSERVED_BENCHMARK_OUTDIR`,
then use `--smoke`, `--prepare`, and `--run`.
The common-panel pilot contains14 already-inspected people, not new
independent donors, and cannot certify seven-source array accuracy.
No production model or PGS calibration has changed.

### Approved dense-supported observed-array benchmark

The user approved the broader existing-data subset. The active
[dense observed-only protocol](./config/observed_dense_protocol.json) and
[runner](./run_dense_observed_benchmark.sh) use3,033 reference people and
240 pure-training people, with all seven sources represented. Whole
pedigree-training families are retained only when every founder is in
that240-person pool:25 of108 families,76 descendants. No unavailable
founder is filled or imputed. This is a smaller training cohort than
the historical307+492 and must not be disguised as an identical-budget
rerun. The [selection manifest](./results_v3_observed_dense_01/selection.json)
records identities, counts and source composition before inference.

Every mask uses all exact typed alleles available in the existing dense
resource, inside the original physical spans and genetic-map bounds,
with common reference-only MAF>=0.005/MAC>=50 filters. The old197,260-site
cap is removed. These quality-filtered observed panels are shared by
FLARE, phase-only processing and freshly fitted Gnomix models.
Source definitions and all identity exclusions remain unchanged.

Five native variants are predeclared: the same three FLARE policies plus
Gnomix LR+CRF at0.5-cM and2-cM target windows. Native EM is per individual,
not fitted across a heterogeneous evaluation cohort. No evaluation
ancestry truth enters marker selection, phasing, model fitting or parameter
selection; there is no learned global calibration in this pilot.

The real chr22 smoke completed all five methods on both masks/phases:

| Mask | Typed SNPs available | Retained after common reference QC | Old common-panel count |
| --- | ---: | ---: | ---: |
| Lonni | 2,245 | 2,225 (99.1%) | 76 |
| Frederik | 12,355 | 12,203 (98.8%) | 319 |

Six focused tests passed. Independent bcftools comparisons verified
28,856 supplied diploid GTs and zero added/untyped SNP keys. Dense
reference/training/reconstructed Ped-sim genotypes also reproduced
1,323,645 original common-marker genotype-and-phase anchors.
The [independent smoke record](./results_v3_observed_dense_01/independent_smoke_verification.json)
documents scope: these are compatibility/integrity checks, not chr22
accuracy rankings or full-autosome certification.

Reproduce with a fresh `DENSE_OBSERVED_OUTDIR`:

```bash
bash run_dense_observed_benchmark.sh --smoke
bash run_dense_observed_benchmark.sh --prepare
bash run_dense_observed_benchmark.sh --run
```

The default output is `results_v3_observed_dense_01`; sources and inputs
are SHA-bound after smoke. The pinned network-disabled container allows
12 CPUs/64GiB and two concurrent chromosome/mask tasks. All22 autosomes
are required before any global comparison is emitted. Each completed
task is checkpointed; interrupted partial tasks require inspection, not
silent overwrite. Monitor its `state.json` and
`docker logs snpster-results_v3_observed_dense_01`.

Reports retain per-chromosome native fit/inference runtimes and panel
coverage, plus per-mask/phase/pure-pedigree global TV, tails, source errors,
off-target mass and5/10-point coverage. The14-person development pilot
does not establish independent seven-source or minor-component accuracy,
and neither array masks nor phase variants count as additional people.
No production or PGS change is made.

### Verified observed-only results (2026-10-18)

The full44-task/22-autosome native-method pilot is complete after an
explicitly recorded recovery. The original container exited after a
native JVM SIGSEGV in a BGZIP/arraycopy routine during the Lonni chr4
non-array FLARE ablation. Its executor completed the other43 tasks;
the original failed `state.json`, stopped container, empty failed-call
log and `hs_err_pid122181.log` are preserved.

The separate [recovery producer](./pedigree_benchmark/scripts/recover_dense_observed.py)
verified frozen inputs and all43 completed task checkpoints, reused71
successful FLARE calls from chr4, executed only13 missing FLARE calls
and the two missing chr4 Gnomix fits, then generated the full comparison.
The failed call succeeded with the same JVM/scientific parameters and
seed under a distinct output prefix. No frozen producer/configuration,
input, phase, variant or reference identity was modified. Recovery is
recorded separately in `recovery.json` and `recovery_state.json`;
the original `state.json` intentionally remains the failed-run evidence.

Independent [verification](./results_v3_observed_dense_01/verification.json)
checks all6,160 chromosome/method predictions, recomputes all280 global
predictions and summary errors from CSV, verifies frozen hashes/donor
exclusions, rereads all13,623,582 supplied typed diploid GTs and exact
alleles across inherited/unphased/phased VCFs, and independently
reaggregates eight FLARE local-posterior files. Every supplied GT is
preserved; no untyped SNP is added.

Coverage is146,080/148,019 typed SNPs (98.69%) for Lonni and
827,033/839,144 (98.56%) for Frederik after shared reference-only QC
and physical-span/map restrictions. This is not the old3-4%-coverage panel.

Primary phase-only observed-array mean TV, lower is better:

| Method | Lonni pure (7) | Lonni pedigree (7) | Frederik pure (7) | Frederik pedigree (7) |
| --- | ---: | ---: | ---: | ---: |
| FLARE array, singleton EM | 15.62% | 12.76% | 18.55% | 15.22% |
| FLARE array, fixed | 17.69% | 18.05% | 18.80% | 17.54% |
| Gnomix LR+CRF,0.5 cM | 2.11% | 20.30% | 1.56% | 13.59% |
| Gnomix LR+CRF,2 cM | 2.34% | 20.68% | 1.51% | 14.85% |

The FLARE `array=false` fixed ablation gives identical saved fractions
to `array=true` fixed in this protocol; the switch alone does not repair
these results. Neither method meets the10% mean pedigree target on both
observed-only masks. All-component10-point pedigree coverage is4/7
(Lonni) and2/7 (Frederik) for FLARE EM, versus3/7 on each mask for
Gnomix0.5 cM. The two phase conditions and masks are repeated measurements
of the same seven pedigrees, not new independent people.

The inherited-phase diagnostic sharply separates the failure modes:

| Method | Lonni pedigree TV, inherited / statistical phase | Frederik pedigree TV, inherited / statistical phase |
| --- | ---: | ---: |
| FLARE array, EM | 9.17% / 12.76% | 11.42% / 15.22% |
| Gnomix LR+CRF,0.5 cM | 0.63% / 20.30% | 0.43% / 13.59% |

Gnomix recognizes these inherited haplotypes very accurately but loses
mixture accuracy after phase-only processing with unchanged GTs; this
is not an imputation-induced genotype error. FLARE's errors already
persist with original phase, especially for the single SAS and GME
pure controls, so phasing alone cannot explain its performance.
For example, SAS donor HG03829 receives only45.6-48.9% SAS and about36.3%
NAM under statistical phase. Even inherited phase leaves35.0-44.5% TV
for that donor. This is a specific diagnostic lead, not proof of its
root cause or of universal FLARE performance.

The next FLARE-focused investigation should isolate native copying,
source-panel structure and parameter initialization on training-only
positive controls before a large hyperparameter sweep. Parameter
learning/calibration must remain separate from native scores and use
donor-purged training evidence. Gnomix remains an informative comparator,
not a production recommendation. No additional training run has been
launched from these scores.

See the [comparison](./results_v3_observed_dense_01/comparison.tsv),
[per-person diagnostic](./results_v3_observed_dense_01/per_person_diagnostic.tsv),
[coverage](./results_v3_observed_dense_01/coverage.tsv) and
[independent verification](./results_v3_observed_dense_01/verification.json).
The14 previously inspected controls, one pure person per source and
seven50:50 mixtures, do not certify all sources, low-level ancestry,
pedigree scenarios or new individuals. Models remain development-only.

**Status (2026-10-13, corrected after founder controls):** The historical sweep scored 80 decoders across 24 shortlisted
representations on 1,148 non-founder individuals from 252 pedigrees and all 22
autosomes. Final unified selection stopped because **no candidate passed the
configured gates**; there are no approved unified finalists. The reference-only
leaderboard is diagnostic, not a replacement for unified selection. Do not
replace an active run's working files. Its pedigree scores are **confounded by
reference-overlapping founders** and must not be treated as the corrected
models' performance. The new default output is
`results_v3_sweep_100_disjoint_reduced/`; the historical output is preserved.

## Phase-independent development candidate found (2026-10-16)

### User-approved accuracy targets (2026-10-17)

The user selected **both** of the following supplemental reporting goals:

1. Mean pedigree TV<=10% is acceptable; <=5% is desirable.
2. Report the fraction of people whose **every** reported source is
   within10 percentage points of truth, and the fraction within5 points.
   These are absolute percentage-point bands, not relative percentage
   errors. No minimum cohort coverage percentage has been specified,
   so coverage is reported rather than assigned a pass/fail gate.

The frozen candidate's currently verified **clean development** results:

| Cohort | Mean TV | Every source within10 points | Every source within5 points |
| --- | ---: | ---: | ---: |
| Pure,159 people | 2.62% | 146/159 (91.8%) | 142/159 (89.3%) |
| Pedigree,123 people | 7.58% | 96/123 (78.0%) | 65/123 (52.8%) |

Thus the mean pedigree error meets the user's acceptable target, but not
the desirable target; it does not imply every pedigree individual is
within10 points. Improving per-person coverage and reducing tail errors
remain important alongside the average.

The earlier1% figure is the **minor-component detection cutoff**, not
a requirement to estimate ancestry within one percentage point. For
example, predicting0% for a true3% component has a3-point error: that
component lies within the desired5-point band even though the separate
detection diagnostic counts it as missed. Small-component recall should
not be confused with these overall accuracy goals.

The targets are persisted in
[phasefree_accuracy_targets.json](./config/phasefree_accuracy_targets.json).
They do not modify any running scientific configuration, frozen
selection gate, model, production QC or PGS calibration.
`bash run_phasefree_accuracy_report.sh` reports both metrics from verified
predictions. It explicitly leaves array results pending until their
independent verification is complete, and then reports each mask/track/
stage separately rather than counting repeated conditions as new people.
For a refreshed report, set `PHASEFREE_ACCURACY_OUTDIR` to a fresh
`results_v3_phasefree_accuracy_targets_*` name.

The initial [accuracy report](./results_v3_phasefree_accuracy_targets_01/accuracy_report.json)
and [per-person errors](./results_v3_phasefree_accuracy_targets_01/per_person.tsv)
cover the existing clean cohorts only. Three focused tests verify
inclusive5/10-point boundaries, the distinction between TV and worst
component error, and rejection of invalid inputs. Separate scalar CSV
calculations reproduced both cohort means and all four coverage counts.
The broader array validation continues without interruption.

The follow-up search completed **25 fixed candidates in four batches**,
retaining all earlier baselines. One candidate,
**`phasefree_stack_boost_clean`**, passes every currently evaluable
unchanged development gate and the prospectively fixed array-control
diagnostic (meanTV<=0.10). Search stopped at this point. It is frozen for
independent validation, **not promoted to production**.

### Results and the important trade-off

Mean TV, lower is better:

| Model | Pure (159) | Original pedigree (123) | Array-imputed (5) |
| --- | ---: | ---: | ---: |
| Frozen LR+CRF | 0.032818 | 0.027295 | 0.243200 |
| Previous clean haplotype-output Ridge | 0.027077 | 0.025121 | 0.232300 |
| Retained fine-frequency shrink50 | 0.091600 | 0.076015 | 0.063573 |
| Reference-augmented PCA/ExtraTrees, leaf1 | 0.021934 | 0.147173 | 0.077036 |
| **Disjoint phase-independent boosted stack** | **0.026171** | **0.075774** | **0.017791** |

The new candidate retains strong pure-source accuracy and removes most
of the observed array-transfer error in these controls. Its clean-pedigree
TV is, however, worse than the best phase-aware decoder's2.51%.
It passes the **unchanged frozen-v1 no-regression thresholds**
(pure0.045825; pedigree0.098466), not a newly invented claim of matching
the best historical result on every track.

The five array-control errors range from1.18% to3.17% TV. They are all
EUR/Asian50:50 mixtures and share some ancestors: their1.78% mean is not
a general seven-source array-performance estimate.

Real-case predictions are EUR44.88%, EAS/SEA48.63%, GME5.75%,
NAM0.19%, SSA0.55%, SAS/OCE0%. Combined off-target mass is6.50%, below
the12% constraint. These satisfy regression bounds, not independently
known ancestry truth. EAS/SEA reporting remains merged; PGS calibration
and production ancestry adjustment remain unchanged.

### Why this differs from the failed global models

The winning model combines two complementary, phase-independent views:

This is **not a Gnomix or FLARE predictor**. Shared container/helper names
are infrastructure, not model provenance. Disjoint-reference discipline
and processing-matched global calibration could also be tested with
FLARE, but this14-feature decoder is not a drop-in decoder for FLARE
outputs, nor evidence that calibration repairs upstream phase distortion.

1. A **reference-only classifier**: 400 ExtraTrees, leaf1, over the exact
   frozen64-PC reference-only representation. Its3,447 training people
   are source/fine-population balanced.
2. The frozen **100-population shrink50 frequency projection**. Its
   populations are aggregated to seven global-source fractions.
3. A **global boosted decoder** over the resulting14 features, fitted on
   the307 clean pure-training people and492 training descendants.
   Pure/source-balanced and pedigree/family-balanced weights are half
   each. Seven independent scalar histogram-boosting regressors use
   300 iterations, learning rate0.05, at most15 leaves, at least10 people
   per leaf, zero L2 penalty and no validation-based early stopping.
   The existing simplex projection produces the final seven fractions.

The classifier never sees any decoder-training person or training
founder: the latter come from the disjoint pure-training pool. Thus the
decoder's classifier features are out of sample without an in-sample
stacking shortcut. Neither held-out genotypes nor truth fit any component.
The winning decoder uses **no imputed training people**. A separately
fixed Ridge/matched-training contrast was tested, but did not pass.

Both views use the same26,066 physically thinned hard-GT markers.
No phase-specific supervision, new imputation, new marker selection or
local-ancestry retraining is involved. Physical thinning is not formal
LD pruning, and complete0/1/2 allele counts are not posterior FORMAT/DS.

### Search record

| Batch | Fixed candidates | Outcome |
| --- | ---: | --- |
| All-marker weighted Ridge | 9: alpha0.1/1/10 x clean/reference/matched policies | Reference augmentation helps; no pure-gate pass |
| Global nonlinear kernels | 9: quadratic/RBF0.5/RBF2 x alpha1e-5/1e-4/1e-3 | None passes all gates |
| Reference-augmented PCA trees/boosting | 4: ExtraTrees leaf1/5, histogram boosting150/400 | Pure accuracy recovers; pedigree gate still fails |
| Disjoint reference-only classifier/frequency stack | 3: clean/matched Ridge, clean boosting | Clean boosting passes all evaluable gates and array diagnostic |

The first two batches retain every standardized genotype marker, rather
than fitting64-PC or population-centroid coefficients. Reference-augmented
policies weight reference50%, pure25% and pedigrees25%; the matched
contrast weights75% of that policy plus25% over the seven actually-imputed
training descendants. The final stacking experiment instead gives its
matched decoder37.5% pure,37.5% pedigree,25% imputed, with **no reference
people in the decoder**.

Kernel Ridge uses an explicit unpenalized-intercept solve and checks
normal equations. Before configuration, focused tests caught and fixed
numerical cancellation in a centered dual solve and rejected accidental
reference inclusion in the stacking matched-policy construction.
A missing historical `/fine` read-only provenance mount also failed
preparation explicitly; its evidence is retained under
`results_v3_supervised_genotypes_01_failed_prepare_01/`. No candidate
fitting occurred in that failed preflight.

Each batch was fixed before its evaluation. The sequence nevertheless
uses the same inspected development cohorts to choose the next
hypothesis: it is adaptive development, not an independent final test.
`results_v3_phasefree_stack_01/overall_comparison.tsv` contains all25 new
candidates and10 retained baselines, including negative results.

### Verified, but not certified

All13 focused tests pass in the pinned scientific container. Independent
verification checked326 distinct hashes for the first22 models, then347
for the final stack and its inputs. It independently reconstructed
training identities, global truth, weights, Ridge/kernel normal equations,
saved predictions and validation errors. The best reference-augmented
tree was independently refitted.

For the final stack, the reference-only classifier and all three decoders
were independently refitted. All1,094 feature rows reproduced within
2.3e-16; held-out decoder refit predictions reproduced **exactly**.
All25 candidates had zero saved-model and phase-invariance prediction
discrepancy across the7,507,008 previously phase-erased/reoriented held-out
genotype cells. Baselines remain unchanged.

The winner's balanced meanTV is0.063830, balanced p95TV0.203084 and
balanced off-target mass0.029496. **The2-5% minor-ancestry gate is not
evaluable**: support is18, below the unchanged minimum20. Observed recall
is7/18=38.89%, below the40% target if it were applied; do not describe
minor-ancestry sensitivity as certified. Pure p95TV remains15.08%;
mean pure GME TV is7.68%, and OCE/NAM validation counts are only3/4.

Next is prospective validation, not more tuning on these cohorts:
additional pure and mixed array controls across sources (especially GME),
adequate minor-ancestry support, other array masks, and faithfully
reproduced production GP/DS-QC missingness. Current inputs are complete
hard calls and the model rejects missing genotypes. Identity disjointness
is not a kinship-disjointness guarantee. No independent final, clinical,
PGS-calibration or general microarray certification is claimed.

### Prospective-validation readiness audit (2026-10-16)

Before selecting another cohort, the audit checked available donor
identities and compatibility with the current production QC module.
Two constraints prevent treating another local run as independent
production certification:

- The union of3,447 fitted reference people,307 pure-training/founder-pool
  people and159 inspected pure-validation people contains3,913 identities.
  The previously eligible founder pool and five historically inspected
  additional donors are also excluded. Of the remaining973 annotated
  identities,698 occur in the dense4,091-person imputation reference.
  **All973 are affinity-only/admixed populations; none has established
  pure seven-source truth.** All dense identities have metadata. These
  people cannot be assigned one-hot source labels merely to create a
  larger test set. New recombinations using existing evaluation donors
  remain seen-donor stress tests, not a donor-independent final test.
- Executing the **unchanged production QC method**, with its current
  GP_MIN0.90/DS_TOL0.05 and biallelic-SNP defaults, on the five saved
  array controls retains21,726-23,752 of26,066 required markers per
  person (83.35-91.12%). It removes2,314-4,340 required markers
  (8.88-16.65%). Consequently, **none of these post-QC inputs satisfies
  the frozen model's complete-genotype contract**.

This does not invalidate the1.78% TV measured on complete imputed GTs.
It means those results do not establish performance on production-filtered
VCFs. No predictions or accuracy estimates were produced from unsupported
incomplete inputs. There is no silent mean filling, marker dropping,
DS-for-GT substitution or production-QC relaxation.

The next implementation choice must be explicit: validate a separate
ancestry branch using complete imputed GTs before strict PGS QC, leaving
PGS QC unchanged, or develop and validate a new missingness-aware model.
The former preserves the frozen predictor but still needs broader
processing/ancestry stress tests; the latter requires new training.
Genuinely independent accuracy validation additionally requires new
source-labelled donor data. No new cohort, training sweep or production
promotion was started by this audit.

Four focused tests cover donor eligibility, reporting-label consistency,
exact sample/allele/marker axes and the actual production-QC boundaries,
including zero-retention behavior. Separate verification checked146 input
and five output hashes, independently parsed all110 preserved VCFs and
reproduced QC for130,330 genotype cells exactly. The candidate digest
remains unchanged. This checks the current QC module's defaults in the
isolated scientific runtime, not a deployed container's environment.

Results: [readiness report](./results_v3_phasefree_readiness_01/readiness.json),
[QC coverage](./results_v3_phasefree_readiness_01/qc_coverage.tsv) and
[independent verification](./results_v3_phasefree_readiness_01/verification.json).
Reproduce with `bash run_phasefree_readiness_audit.sh`; an existing result
directory is rejected. Set `PHASEFREE_READINESS_OUTDIR` to a fresh
`results_v3_phasefree_readiness_*` name for a separate rerun. The runner
uses the pinned scientific image, no network, and read-only scientific
and production-QC sources. It does not import or run upstream Gnomix.

### Approved frozen-model array stress test

The user selected broader seen-donor tests on complete imputed GTs,
not a new missingness-aware training sweep. The fixed design is under
`results_v3_phasefree_array_validation_01/`:

- Seven pure controls, one per source, selected by a fixed hash from
  dense-panel-available evaluation donors.
- Seven50:50 descendants from distinct families and14 distinct donors,
  with every source appearing twice. The five previous array-control
  individuals are not selected again.
- Each of the14 people is tested under both Lonni's148,019-site and
  Frederik's839,144-site observed masks:28 array conditions,616 singleton
  imputation fits across22 autosomes, using19 previously inspected donors.
  The same person under two masks is not two independent individuals.
- The3,731-person imputation reference follows the previous exclusion
  policy. Every selected donor is disjoint from the3,447 fitted ancestry
  reference people and307 training/founder-pool people and is excluded
  from the imputation reference.

No predictions guide selection. The model,26,066-marker axis, reference
policy, seed and meanTV<=0.10 array diagnostic are fixed before scoring.
Dense reconstruction must reproduce both original canonical model-panel
haplotypes exactly. Unphased input sorts allele pairs, preserving dosage
without leaking the original orientation. Every supplied observed
genotype, not just observed model markers, must survive imputation
unchanged. Complete hard GTs are used before strict QC, without
substituting posterior DS or padding missing markers.

This bounded stage tests seven-source pure and50:50 transfer under a
second mask. **It has zero2-5% minor-ancestry support** and cannot certify
that gate or improve its previous18-contribution support by counting masks as
new people. It also cannot establish donor/kinship-independent accuracy,
PGS calibration, or broad pure-population performance from one person per
source.

Twelve focused tests pass, including frozen-stack tests, input/sample axes,
phase erasure, production-QC boundaries, the independent verification VCF
reader and rejection of changed checkpoints. The Java21 runner is
`bash run_phasefree_array_validation.sh --prepare`, then `--smoke` for
chromosome22 and `--run` for all autosomes. It runs two chromosomes at a
time, each capped at6 CPUs/30GiB, with network disabled and source inputs
read-only. Each completed sample/chromosome is hash-verified and reusable.
`--resume` verifies completed work and preserves any incomplete case in
an `.incomplete_*` directory before explicitly retrying it. Incomplete
chromosome preparation still requires inspection.

The final stage will write `predictions.tsv`, `metrics.tsv` and
`completion.json`, comparing each imputed result with its paired clean
prediction. No training or production/PGS-QC change is part of this run.
The chromosome22 preflight completed all28 conditions. Independent
re-reading reproduced4,312 original and8,624 imputed model-marker genotype
cells and confirmed204,400 supplied observed genotypes were preserved.
`preflight_verification.json` records these checks; no ancestry accuracy
was scored from this single chromosome. The full two-chromosome scheduler
was then started, with a separate post-completion verification job queued.
The separate `pedigree_benchmark/scripts/verify_phasefree_arrays.py`
must then run in the same pinned container/mount layout. It independently
rereads original and retained imputed GTs, checks supplied observed calls,
reproduces frozen-model predictions, TV/off-target metrics and diagnostic
decisions, and writes `verification.json`. Its full-result verification
cannot run until all chromosomes complete.

#### Durable recovery (2026-10-17)

The first shell-driven full run did not complete overnight. Chromosomes2
and22 completed all28 conditions each; chromosome1 stopped during
reference preparation. Its compressed cache fails `gzip -t` with an
unexpected end of file. The old scheduler/verification PIDs and worker
containers are absent, and the old shell logs do not establish the cause.
These facts do not support an ancestry-accuracy conclusion.

`bash run_phasefree_array_durable.sh --recover` archives the stale lock
only after checking recorded host PIDs and active validation containers.
Inside one detached, daemon-owned Docker container, the new scheduler
hash-verifies and reuses completed chromosomes, preserves incomplete
chromosome preparation in a named archive, and runs at most two pending
chromosomes concurrently. The container is capped at12 CPUs/64GiB and
retained after exit for inspectable exit status and logs. It neither
modifies the frozen scientific producer/configuration nor changes the
model, masks, donors, reference exclusions or acceptance diagnostic.

Worker output goes directly to `durable_chrN.log` files rather than
through a shell/terminal pipe. Persistent `durable_scheduler.json` records
active/completed chromosomes, archived attempts and any failed child
exit codes. Scoring and independent final verification execute inside
the same container, not in a separate host process waiting on a PID.
Successful completion requires a matching final `verification.json`;
completed imputation alone is not treated as verified accuracy.

Inspect with
`docker logs snpster-results_v3_phasefree_array_validation_01-durable-01`
and the output's `durable_scheduler.json`. The launch script refuses an
existing same-name container instead of overwriting its failure evidence.
Fifteen focused tests pass, including safe archival, preservation of
completed work and persistent nonzero child-failure reporting.

### Targeted frozen-candidate error diagnosis (2026-10-17)

While the recovered array run proceeds, a separate read-only diagnostic
reproduces all282 pure/pedigree predictions from the frozen decoder and
exact saved14-feature values. It changes no models or thresholds.
Independent CSV-based verification checks1,974 person/source cells,
training support and all14 input/three output hashes.

The existing2-5% support is **18 ancestry contributions in14 people from
10 families**, not18 independent people. At the unchanged1% detection
threshold the final decoder detects7/18, the raw classifier8/18 and the
raw frequency view11/18. The minimum support remains20 contributions,
so the minor gate is still unevaluable.

Small SAS contributions are the clearest diagnostic failure: the decoder
detects0/7 across six families, while each raw view detects3/7. The clean
training cohort has only11 people with2-5% SAS, and only six each for GME
and SSA. These are support counts, not independent-founder counts or
proof of a causal training defect.

Raw-view recall alone is not a fix: among1,465 contributions whose truth
is below1%, the final decoder flags171 at>=1%, compared with280 for the
classifier and308 for frequency. The two raw views are uncalibrated
features, not validated replacement models. Lowering the threshold or
forcing frequency estimates into the answer would risk more spurious
small-component calls; neither change was made.

The audit also caught a **feature-export numerical contract issue**:
pandas' default TSV parser changes tiny floating-point values at
histogram-tree split boundaries. Compared with the saved predictions,
default parsing changes seven pure and50 pedigree answers, with a largest
single-component shift of5.06 percentage points. Exact round-trip parsing
reproduces all282 answers within1.2e-16. The diagnostic now uses
`float_precision="round_trip"`, covered by a regression test.
Original training and the active array evaluator construct features
directly in memory, so this does not invalidate their predictions.
Future feature caching must preserve exact values; rounding the14
features is not a supported deployment shortcut. No quantization or
retraining was applied to the frozen model.

The next bounded training experiment tests calibration while preserving
pure-source, pedigree, off-target and array performance. Training-only
cross-validation can proceed alongside the array run; final comparison
must wait for independently verified two-mask results. Better balanced
low-level mixture training remains a subsequent hypothesis if calibration
alone is insufficient. Neither hypothesis proves the cause or that
accuracy has already improved. Genuinely independent certification still
requires additional donor data.

Reproduce with `bash run_phasefree_error_audit.sh` and a fresh
`PHASEFREE_ERROR_AUDIT_OUTDIR` if rerunning. Results are
[error audit](./results_v3_phasefree_error_audit_01/error_audit.json),
[minor contributions](./results_v3_phasefree_error_audit_01/minor_contributions.tsv),
[interpretation](./results_v3_phasefree_error_audit_01/interpretation.json)
and [independent verification](./results_v3_phasefree_error_audit_01/verification.json).
All20 focused phase-independent tests pass in the pinned scientific
container. The editor diagnostic/test-discovery tools remain unavailable;
runtime tests are the validation evidence, not a claimed editor type-check.

### Standardized donor-purged calibration sweep (2026-10-17)

The fixed [sweep protocol](./config/phasefree_sweep_protocol.json) tests
**nine new calibration variants plus an architecture-matched refit
control**, without rerunning FLARE, Gnomix or imputation. The frozen
reference-only classifier and fine-frequency views supply the same14
features. Only the decoder is retrained on the existing307 pure people
and492 training descendants; no imputed or development people enter
decoder training.

The candidates compare squared/absolute-loss histogram boosting,
smoother/flexible trees, linear-residual correction, quadratic Ridge and
linear/boosted blends. All new candidates round features to12 decimals
at both fit and prediction, with one-ULP perturbation and serialization
checks. This is a newly trained input transform, **not a modification
or supported quantization shortcut for the frozen baseline**.

Validation and selection are standardized as follows:

1. **114 training-only folds:**108 leave-family-out pedigree folds plus
   six source-stratified pure-donor folds. Every one of the799 training
   people receives exactly one out-of-fold prediction.
2. **Founder purging:** remove held-out donors from pure training and
   remove every other training family sharing any held-out founder.
   Ordinary coarse family folds were not used because reused founders
   leave almost no pedigree training after strict purging. These folds
   are donor-purged, not independently sampled families or a
   kinship-independent certification cohort.
3. Recompute50% source-balanced pure/50% family-balanced pedigree
   training weights within each fold.
4. Rank on the existing scenario/combination/family-balanced TV; use
   pedigree mean TV and5/10-point coverage as fixed tie breakers.
   **Lock two new finalists using training CV only**, before evaluating
   them and the control on the existing development cohorts.
5. Retain historical pure/pedigree, balanced-tail and off-target
   guardrails. Report the user's mean-TV10% acceptable/5% desirable
   targets, every-source5/10-point coverage, source/scenario errors and
   minor-component diagnostics. No minimum coverage gate is invented.
6. Wait for independent shared-array verification, then reuse the
   existing genotype matrices for both masks. Report clean/imputed,
   pure/mixed and mask conditions separately;28 conditions are14 people,
   not28 independent people. No new imputation is launched.

An independent preflight audited all114 train/test donor splits and
exactly reproduced all1,094 original saved-feature predictions from a
full-data control refit (maximum difference0). All ten models also
passed one complete real-fold smoke fit. Preparation and smoke checks
are not evidence of improved ancestry accuracy.

**Verified CV/development outcome:** all114 folds completed, all799
people received exactly one OOF prediction, and independent donor,
metric, model-refit and saved-prediction checks passed with zero
reproduction error. Training CV locked `stable_squared` and
`stable_blend25` before their development scores were read.

| Model | Training-CV balanced TV | Clean pure TV | Clean pedigree TV | Pedigree within10 points | Pedigree within5 points |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen baseline / exact refit control | 8.56% | 2.62% | 7.58% | 78.0% | 52.8% |
| `stable_squared` | 8.21% | 2.48% | 7.57% | 78.0% | 52.0% |
| `stable_blend25` | 8.65% | 3.12% | 7.69% | 81.3% | 52.0% |

CV and development columns describe different cohorts and weighting;
they are not interchangeable estimates. The new squared model's
donor-purged pedigree mean TV is9.28%, versus9.63% for the refit control.
It improves the training-CV ranking modestly but barely changes clean
development pedigree TV. The blend increases10-point coverage by3.3
percentage points while worsening mean TV, pure-source error and the
original five-array diagnostic. Both finalists pass the unchanged
development gates and original five-array mean-TV diagnostic, but
**neither establishes the5% pedigree target or an overall accuracy
breakthrough**.

The [CV comparison](./results_v3_phasefree_calibration_01/cv_comparison.tsv),
[development comparison](./results_v3_phasefree_calibration_01/development_comparison.tsv)
and [independent verification](./results_v3_phasefree_calibration_01/verification.json)
are persisted. The durable sweep is now waiting for verified shared
arrays; it will perform the locked comparison automatically. No model
has been promoted. Rather than choosing further variants by these
already-inspected development scores, a subsequent substantial
improvement should test richer phase-independent representations or
better balanced, independently founded training mixtures under this
same fixed protocol.

The [runner](./run_phasefree_calibration_sweep.sh) uses the pinned,
network-disabled scientific container, four CPUs and16GiB RAM, with two
CV workers. Individual fold predictions and completion records support
checkpointed recovery. Full-data model bundles, locked selection,
development scores, independent verification and the eventual shared
array comparison are retained under
[results_v3_phasefree_calibration_01/](./results_v3_phasefree_calibration_01/).
The final recommendation requires training-CV improvement, unchanged
development gates, the original five-control array diagnostic, all
expanded mask-by-track mean-TV diagnostics<=0.10 and the user's pedigree
mean-TV<=0.10. It can recommend further validation, never production
promotion.

For a fresh experiment, use a new `PHASEFREE_CALIBRATION_OUTDIR` and run:

```bash
bash run_phasefree_calibration_sweep.sh --prepare
bash run_phasefree_calibration_sweep.sh --smoke
bash run_phasefree_calibration_sweep.sh --run
```

The default prepared experiment must not be prepared again or overwritten.
Inspect the durable run with
`docker logs snpster-results_v3_phasefree_calibration_01` and its
`sweep_state.json`. The runner refuses an existing same-name container;
review retained exit status/logs before explicitly restarting it.
The shared array scheduler continues independently. Final reporting
waits at most48 hours for its verified completion and fails explicitly
if the source fails.

Six focused sweep tests cover donor exclusion/OOF coverage, fold weights,
invalid metadata, numerical stability, canonical JSON ordering of every
candidate kind, serialization and balanced family scoring. No production
QC, PGS calibration, frozen baseline or running array configuration is
changed. Current local data still lack unused source-labelled donors,
post-QC missingness support and independent final-test certification.

### Richer phase-independent features, current dataset only (2026-10-17)

The next fixed [richer-feature protocol](./config/phasefree_rich_protocol.json)
uses **only the existing dataset**. No external samples are downloaded,
no new people or pedigrees are simulated, and no additional imputation
is launched. The original307 pure and492 pedigree training people,
114 donor-purged folds, source/family-balanced weights, accuracy targets
and historical guardrails are retained.

Earlier direct PCA, fine-frequency, full-marker Ridge and kernel models
did not pass all gates. This experiment instead augments the successful
stack's14 inputs with **64 frozen reference-only PCA coordinates and
100 unaggregated fine-population mixture weights**, giving178 features.
The PCA, classifier and reference-frequency profiles are not refitted.
Feature extraction may transform evaluation genotypes, but evaluation
truth or scores never enter training, preprocessing or finalist ranking.

Six fixed heads comprise a stable14-feature control, three boosted
ablations (base+PCA, base+fine weights, and all features) and two
conservative blends of the control with a richer Ridge head (25%/50%
linear contribution). Ridge centering/scaling is learned only from
each fold's training rows and weights. New PCA inputs are rounded to
eight decimals; base/fine fractions to12, at fit and prediction.
Signed PCA coordinates are retained. Missing/invalid genotype or
feature values fail explicitly; there is no padding or DS substitution.

The stable control must reproduce the preceding sweep's
`stable_squared` OOF predictions in **every fold**, providing a matched
comparison rather than changing the baseline during representation
testing. Each model is also checked for serialization and one-ULP
feature stability. Seven focused tests pass in the pinned container;
a complete six-model real-fold smoke fit passed in8.65 seconds.
These are software/preflight checks, not accuracy improvements.

Training CV locks two new finalists before development evaluation.
An independent verifier re-audits donors, OOF coverage, balanced TV,
accuracy bands, saved predictions and separate pedigree/pure-fold
refits. Afterward the durable runner waits for the existing shared
array verification, extracts the richer features from the retained
genotypes and checks both masks without repeating imputation.
Final eligibility requires training-CV improvement over the stable
control and the unchanged development/user/array guardrails.
No automatic production or PGS promotion is permitted.

**Verified outcome:** all six models completed114 folds and all799
training people received one OOF prediction. Independent donor,
metric, saved-model and refit checks passed; spot-check differences
were at most2.3e-16. Every stable-control OOF prediction reproduced
the preceding standardized sweep.

| Model | Training-CV balanced TV | Training-CV pedigree TV | Clean development pedigree TV | Development within10 points | Development within5 points |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stable14-feature control | **8.21%** | **9.28%** | 7.57% | 78.0% | **52.0%** |
| Fine-population richer boost | 8.56% | 9.93% | **7.35%** | **82.1%** | 50.4% |
| Combined25% linear blend | 8.82% | 9.66% | 7.52% | **82.1%** | **52.0%** |

The two new finalists pass the unchanged development gates, but
**none of the five new models improves the predeclared training-CV
ranking over the stable control**. The fine-weight model's attractive
development mean/10-point coverage is therefore not sufficient to
replace the baseline;5-point coverage also decreases. Do not select
it by the already-inspected development cohort or claim independent
improvement. See the [CV table](./results_v3_phasefree_rich_01/cv_comparison.tsv),
[development table](./results_v3_phasefree_rich_01/development_comparison.tsv)
and [independent verification](./results_v3_phasefree_rich_01/verification.json).

During this sweep, the shared array run stopped in chromosome4
preparation, and both dependent sweeps correctly stopped rather than
scoring incomplete inputs. All1,855 chr4 model markers are present.
The cause is one coordinate with two **different valid exact alleles**:
the model requires`4:164203107:G:A`, while the Frederik mask observes
`4:164203107:G:T`. Both exist in the dense reference. The old global
unique-position assertion and position-keyed extraction cannot
represent that union; substituting either allele or dropping a model
marker would be incorrect.

An explicit [exact-allele recovery](./pedigree_benchmark/scripts/recover_phasefree_alleles.py)
keeps both variants, checks the complete canonical haplotypes and
every supplied observed genotype, and extracts retained calls by
`(CHROM,POS,REF,ALT)`. It preserves the original frozen producer,
manifest, masks, reference exclusions and completed chromosomes.
Its separate`allele_recovery.json` fingerprints the recovery code and
is included in new preparation/case input provenance. Twelve focused
array/accuracy/allele tests pass. A real Frederik chr4 imputation
passed in66.7 seconds. A separate`bcftools` reread verified all52,773
supplied observed genotypes, all1,855 model genotypes and the presence
of both collision alleles, without using the recovery loader.
The [independent smoke verification](./results_v3_phasefree_array_validation_01/allele_smoke_verification.json)
also confirms the original frozen model hash.
The [recovery runner](./run_phasefree_allele_recovery.sh) requires this
real-case smoke stage before durable scheduling:

```bash
bash run_phasefree_allele_recovery.sh --prepare
bash run_phasefree_allele_recovery.sh --smoke
bash run_phasefree_allele_recovery.sh --run
```

The recovered run reuses completed chromosomes1,2,3,5,22 and the smoke
case; it does not restart their imputation. Only after the array
scheduler is running again should the stopped calibration and richer
containers be explicitly restarted. Their verified CV/full-fit
checkpoints are reused; they still require independent shared-array
verification before reporting transfer results.
Recovery is now running in
`snpster-results_v3_phasefree_array_validation_01-alleles-02`
with chromosomes4/6 active. Both dependent sweep containers have
been restarted from their completed CV/development checkpoints.

Run a fresh experiment with:

```bash
bash run_phasefree_rich_sweep.sh --prepare
bash run_phasefree_rich_sweep.sh --smoke
bash run_phasefree_rich_sweep.sh --run
```

The [runner](./run_phasefree_rich_sweep.sh) defaults to
[results_v3_phasefree_rich_01/](./results_v3_phasefree_rich_01/),
four CPUs,16GiB RAM and two CV workers in a pinned, network-disabled
container. For a distinct experiment use a fresh
`PHASEFREE_RICH_OUTDIR`; never overwrite prepared/frozen producers.
Inspect progress with
`docker logs snpster-results_v3_phasefree_rich_01`
and the output's `sweep_state.json`. Completed folds are checkpointed.

Additional simulations from these same founders would not create new
independent donor evidence. The limited OCE/founder diversity,
already-inspected development cohorts, complete pre-QC GT contract and
lack of unused source-labelled local donors still constrain what can
be claimed, regardless of whether richer features improve performance.

### Balanced training mixtures, unchanged donor-purged validation

The fixed [balanced-mixture protocol](./config/phasefree_balanced_protocol.json)
changes training support rather than adding more decoder features. It uses
only the existing307 pure-training founders and the same stabilized14
classifier/frequency fractions. Reference-only preprocessing is frozen;
no samples are downloaded and no imputed data enter training.

The independently verified bank contains2,626 training-only mosaics in371
scenario/source-role/fraction strata:50:50 F1,75:25 asymmetric,2.5/3/4/4.5%
minor components,50:25:25 three-way and25:25:25:25 four-way mixtures.
Donor sampling balances fine populations within each source. An initial
1,428 plans receive donor-safe coverage repairs using only training-fold
identities, with a fixed4,096-plan ceiling. Every stratum retains at least
two safe mosaics in every fold.

These are ancestry-training mosaics, **not new independent people or a
Ped-sim/meiosis benchmark**. Exact physical-base-pair ancestry truth is
recomputed from saved half-open segments across both chromosome copies,
not inferred from SNP counts. Complementary founder haplotypes preserve
the original diploid genotype wherever both copies use the same donor.
A regression test caught artificial homozygosity in the initial generator
before data preparation; the corrected generator passed all seven focused
tests. The independent [bank verification](./results_v3_phasefree_balanced_01/augmentation_verification.json)
reconstructed all68,449,316 mosaic genotypes and their physical truth.
All22 founder chromosome extracts also matched the frozen dosage caches.

Five models retain identical boosted-head settings:

| Model | Original pure mass | Original pedigree mass | Augmentation mass |
| --- | --- | --- | --- |
| Stable control | 50% | 50% | 0% |
| Balanced20 | 40% | 40% | 20% |
| Balanced40 | 30% | 30% | 40% |
| Balanced60 | 20% | 20% | 60% |
| Minor-half40 | 30% | 30% | 40% |

Augmentation weights balance scenarios, then strata, rather than row
counts. Minor-half40 puts half its augmentation mass on minor mixtures;
the other policies give each of the five scenarios equal mass.

The original799-person OOF cohort and114 donor-purged folds are unchanged.
Every augmented sample sharing any held-out donor is removed. Ranking
uses only original-person OOF predictions under the established balanced
TV/coverage ordering; synthetic training accuracy cannot select a winner.
Two new finalists are locked before development evaluation. Historical
pure-source and other guardrails remain unchanged. The five-model real-fold
smoke passed, including exact stable-control reproduction, serialization,
simplex output and one-ULP stability.

Run a fresh experiment with:

```bash
bash run_phasefree_balanced_sweep.sh --prepare
bash run_phasefree_balanced_sweep.sh --smoke
bash run_phasefree_balanced_sweep.sh --run
```

The [runner](./run_phasefree_balanced_sweep.sh) defaults to
[results_v3_phasefree_balanced_01/](./results_v3_phasefree_balanced_01/),
four CPUs,16GiB RAM and two CV workers in the pinned, network-disabled
container. Use a fresh `PHASEFREE_BALANCED_OUTDIR` for a distinct experiment;
do not modify frozen producers or overwrite completed results. Monitor
`docker logs snpster-results_v3_phasefree_balanced_01` and the output's
`sweep_state.json`. Verified fold checkpoints are reusable.

After CV/development verification, the sweep waits for the existing shared
two-mask array validation and reuses its complete genotype outputs without
additional imputation. No production model, GP/DS-QC policy or PGS
calibration is changed. More mosaics do not increase independent founder
diversity or establish a new final certification cohort.

#### Verified CV and development results

All114 folds and799 original-person OOF predictions are complete.
Independent verification confirms zero original/augmented held-out donor
overlaps, exact stable-control reproduction, recomputed metrics and exact
refits of all three evaluated models on one family and one pure-donor fold.
The two finalists were locked using training CV only.

Minor-half40 ranks first: balanced CV TV improves from8.205% to7.115%;
its pedigree CV TV improves from9.285% to8.098%. Balanced40 ranks second,
with balanced CV TV7.148% and pedigree CV TV8.036%.

| Model | Clean pure TV | Clean pedigree TV | Pedigree within10 points | Pedigree within5 points | Original five array controls TV |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stable control | 2.480% | 7.569% | 78.0% | 52.0% | 1.538% |
| Minor-half40 | 3.079% | 6.088% | 84.6% | 61.0% | 0.997% |
| Balanced40 | 3.221% | 6.096% | 85.4% | 60.2% | 1.380% |

This is a meaningful pedigree/mixture improvement, with a pure-source
trade-off rather than a universal gain. All three pass the unchanged
evaluable guardrails. The five old array controls are already-inspected,
related EUR/Asian mixtures, not independent seven-source certification.
The expanded two-mask array evaluation remains pending, and no production
promotion has occurred.

See the [CV comparison](./results_v3_phasefree_balanced_01/cv_comparison.tsv),
[development comparison](./results_v3_phasefree_balanced_01/development_comparison.tsv)
and [independent verification](./results_v3_phasefree_balanced_01/verification.json).

### Persisted candidate and reproducibility

`results_v3_phasefree_stack_01/selected_candidate.json` freezes the model,
image, experiment and verification digests. The self-contained bundle is
`candidates/phasefree_stack_boost_clean/model.pkl`; its SHA256 is
`95f3e3c9abaf76fb357d4ef51b2e29bf212a033002c02a04b825ec01b0231686`.
All component models, marker keys and source ordering are saved.
`interpretation.json` documents the decision and its limitations.

Run each relevant shell runner with `--prepare`, then `--smoke` where
supported, then `--run`:

- `run_supervised_genotype_screen.sh` -> `results_v3_supervised_genotypes_01/`
- `run_kernel_genotype_screen.sh` -> `results_v3_kernel_genotypes_01/`
- `run_reference_tree_screen.sh` -> `results_v3_reference_trees_01/`
- `run_phasefree_stack.sh` -> `results_v3_phasefree_stack_01/`

All stages use fresh network-disabled containers with scientific sources
and baselines read-only. Existing completed outputs must not be
overwritten. Full fits are short bounded batches; incomplete outputs
require explicit inspection rather than silent reuse. These original
screens are complete; subsequent training sweeps and the separately
approved frozen-model array validation are described above.

## Fine-population frequency profiles (2026-10-16)

The four completed genotype-screen models are retained as immutable
development baselines. A separate completed screen tests whether finer
reference profiles fix the broad-frequency model's pure-source errors.
It reuses the verified dosage matrices without further imputation.

All **100 fine populations** in the 3,447-person reference supply one
allele-frequency profile each. For a population with `n` reference people,
the fixed profiles are:

```text
raw:       fine_AF
shrink10:  n/(n+10) * fine_AF + 10/(n+10) * broad_AF
shrink50:  n/(n+50) * fine_AF + 50/(n+50) * broad_AF
```

These three settings were frozen before evaluation. Counts and broad
frequencies use reference individuals only; neither training nor held-out
genotypes or labels influence profile construction. All three use the
same 26,066 physically thinned marker panel and the exact broad-source
AF-range>=0.03 informative mask from the earlier frequency baseline.
This is not formal LD pruning.

For each person, nonnegative population weights sum to one and minimize
squared allele-frequency error. Broad-source fractions are sums over
their populations; no new intercept or supervised calibration is fitted.
Initialization assigns equal mass to each broad source, avoiding a larger
initial prior for sources with more populations. SLSQP uses
`ftol=1e-12`, at most1,000 iterations, and an explicit convex-simplex
dual-gap check (maximum1e-6). Subtracting a shared per-marker center
improves numerical conditioning without changing this constrained
objective. The older helper is seven-profile-specific, so this separate
generalized solver does not modify that immutable baseline implementation.

### Completed comparison

Mean TV, lower is better:

| Model | Pure | Original pedigree | Matched array-imputed |
| --- | ---: | ---: | ---: |
| Previous clean haplotype-output Ridge | 0.027077 | 0.025121 | 0.232300 |
| Retained broad-frequency baseline | 0.116679 | 0.088692 | 0.073522 |
| Fine profiles, raw | 0.099333 | 0.088586 | 0.089931 |
| Fine profiles, shrink10 | 0.090414 | 0.078049 | 0.077591 |
| Fine profiles, shrink50 | 0.091600 | 0.076015 | 0.063573 |

Finer profiles with shrinkage improve the broad-frequency model's
pure/pedigree trade-off. Shrink50 also improves array-control TV to
**6.36%**, compared with7.35% for the broad-frequency baseline and24.32%
for the frozen CRF under the same imputation-reference policy.
Nevertheless, **all three fail pure no-baseline-regression and the real
off-target gate**; none is promoted. Their real combined GME/NAM/OCE/SSA
fractions are16.20%,15.85%,15.51%, versus the unchanged12% ceiling.
Reasonable EUR/Asian balance alone is not sufficient.

The hoped-for GME fix did not materialize: raw fine profiles increase
pure GME TV from27.55% to36.57%; shrink50 returns it to27.20%.
Shrink10 improves SAS from20.59% to12.86% and EUR from13.21% to10.70%,
but this remains far above the required overall pure accuracy.
These results support retaining phase-independent baselines, not
replacing the prior models or launching a broad profile/shrinkage sweep.

All six focused tests pass, including exact toy pure/mixed recovery,
count-based shrinkage, reference-only construction, source-balanced
initialization, persistence, phase invariance and explicit optimizer
failure. The real-data compatibility check used one clean training person
per source, with no held-out evaluation. The full fit/prediction stage took
about5 minutes.

Independent verification checked **159 input/output hashes**, recomputed
all100 profiles from reference genotypes, reconstructed every broad
prediction from population weights, independently checked convex
optimality gaps, and recomputed all validation TVs and real gates.
Twelve additional tighter-tolerance, uncentered solver checks agree in
objective within1e-7 and in broad fractions within7.4e-5. Population-weight
aggregation agrees with saved predictions within1.6e-15.
All7,507,008 held-out phase-erased/reoriented genotype cells were reused;
every candidate's phase and serialization prediction discrepancy is
exactly zero.

```bash
bash run_fine_profile_screen.sh --prepare
bash run_fine_profile_screen.sh --smoke
bash run_fine_profile_screen.sh --run
```

Results are in `results_v3_fine_profiles_01/`; the comparison includes all
four retained genotype baselines. Inputs and baselines are mounted
read-only in fresh network-disabled containers using the original pinned
image. Candidate checkpoints are hash-verified for interrupted-run reuse;
completed results cannot be overwritten. See `interpretation.json` and
`verification/final_verification.json` for conclusions and verification.
Population weights are fitting coefficients, not certified fine ancestry.
The same development-cohort, shared-donor five-control, single-array and
hard-GT/QC limitations from the genotype screen still apply.

## Phase-independent genotype-model screen (2026-10-16)

This completed four-candidate screen reuses the verified matched-imputation
inputs; no further imputation or local ancestry fitting is performed.
It predicts global ancestry from **hard-GT alternate-allele counts**,
accepting phased or unphased calls. These are not posterior-mean FORMAT/DS
dosages. Missing, partial or non-biallelic calls fail explicitly.

The fixed common panel is physically thinned to **26,066 SNPs across all
22 autosomes**: the first available marker in each 100kb bin is retained
without consulting evaluation genotypes or labels. This reduces marker
density imbalance but is **not formal r-squared LD pruning**. Exact
CHROM/POS/REF/ALT identities and sample ordering are verified.

| Candidate | Fixed method |
| --- | --- |
| `frequency_simplex_100kb` | Seven reference source allele-frequency profiles; nonnegative least-squares projection constrained to sum to one, reusing the transfer-audit fitter |
| `dosage_pca64_ridge_clean` | Reference-only allele scaling and 64 whitened PCs; clean-training Ridge |
| `dosage_pca64_ridge_matched` | Same preprocessing and Ridge, with actual imputed training |
| `dosage_pca64_trees_matched` | Same preprocessing; 200 ExtraTrees, minimum leaf size5, all PCs available per split |

Allele centering, binomial allele-frequency scaling and randomized PCA are
fitted on the 3,447 reference individuals only. The supervised heads use
the existing 307 pure-training people and 492 training descendants; the
matched heads also use the seven verified, balanced imputed training
descendants. Clean-only source/family weight is 50%/50%; matched weights
are 37.5%/37.5%/25%. No held-out individual or truth enters preprocessing,
head fitting or hyperparameter selection.

Both Ridge heads use fixed `alpha=0.01` and an intercept. Unlike the
previous probability-feature heads, the PCs are centered coordinates:
an intercept is needed to represent ancestry proportions around that
origin. Trees use the fixed seed and one worker for reproducible inference.
All models produce the existing seven-source output; simplex projection
uses the existing helper. No local ancestry, production PGS adjustment,
reporting labels or acceptance gates change.

### Completed results

Mean TV, lower is better:

| Model | Pure | Original pedigree | Matched array-imputed |
| --- | ---: | ---: | ---: |
| Previous clean haplotype-output Ridge | 0.027077 | 0.025121 | 0.232300 |
| Source-frequency simplex | 0.116679 | 0.088692 | 0.073522 |
| Dosage PCA/Ridge, clean | 0.115549 | 0.095062 | 0.126774 |
| Dosage PCA/Ridge, matched | 0.119122 | 0.101432 | 0.140404 |
| Dosage PCA/ExtraTrees, matched | 0.065077 | 0.172847 | 0.119679 |

The frequency baseline reduces array-control mean TV from the frozen
CRF's `0.243200` to `0.073522`, but its pure-source error is too high.
Its real EUR/Asian estimates are 40.91%/46.53%; it narrowly fails the
combined off-target constraint (12.22% versus maximum12%). Both dosage
Ridge models pass all real-case constraints. The clean Ridge estimates
45.18% EUR and45.89% Asian, but fails pure no-baseline-regression.
The matched Ridge also fails the pedigree regression gate. Trees have the
lowest pure TV among the new models but worse pedigree performance and
fail off-target and real-case gates. **No candidate passes all gates;
none is promoted.**

The pure-frequency errors are largest in GME (27.55%), SAS (20.59%) and
EUR (13.21%). This suggests testing whether finer within-source reference
profiles improve pure-source fit while retaining phase invariance; it is
a hypothesis, not a verified fix. Do not tune around the one real sample
or launch a broad sweep based only on the five array-imputed controls.
The representation and marker panel also changed, so the accuracy
difference cannot be attributed solely to removing phase dependence.

Independent verification checked **725 distinct input/output hashes**,
refitted PCA on reference individuals only, and refitted all three
supervised heads using training-only features/targets. PCA components
reproduced exactly; head predictions reproduced within `4.4e-15`.
All **7,507,008 held-out genotypes** were phase-erased and alternately
reoriented: every model's prediction difference was exactly zero.
Serialized predictions and all three TV means were independently
recomputed. Five focused tests cover actual phased/unphased VCF reading,
sample/allele ordering, missing-call rejection, persistence and held-out
feature exclusion from training.

```bash
bash run_dosage_model_screen.sh --prepare
bash run_dosage_model_screen.sh --smoke
bash run_dosage_model_screen.sh --run
```

The default isolated output is `results_v3_dosage_screen_01/`.
Chromosome extraction runs two jobs concurrently in fresh network-disabled
containers; prior inputs and models are read-only. Completed extraction
checkpoints are hash-verified and reused. The final fit is one fresh
container and must not overwrite a completed screen. See `completion.json`,
`comparison.tsv`, `interpretation.json` and `verification/final_verification.json`
there for evidence. This remains development validation, with sparse
OCE/NAM support, seven imputed training genomes, five shared-donor controls
and one observed-array mask; no independent final certification is claimed.

## Matched-imputed global-decoder training (2026-10-16)

After rejecting the two synthetic phase-switch recipes, the user approved
a small **global-ancestry decoder**, with the local LR+CRF models frozen.
This distinction avoids an invalid training shortcut: after imputation
changes phase, the original ancestry labels do not necessarily belong to
the same numbered haplotype. Individual-level global ancestry truth remains
valid and is the only supervised target used in this experiment.

Dense reference genotypes are available for all ancestors of 194 existing
training descendants. A deterministic integer-program selection chooses
seven two-source descendants from seven families, using 14 distinct
training founders. Each individual is 50/50 between two merged sources;
each of the seven sources appears twice, giving exactly balanced ancestry
exposure. Selection uses donor availability, family/donor identity and known
training truth, never validation predictions. All founders come from the
existing pure-training pool and are disjoint from the ancestry reference
and held-out individuals.

The original saved Ped-sim breakpoints extend these same inherited genomes
to the observed array markers. They must exactly reproduce the original
canonical phased model genotypes before proceeding. All seven training
individuals and five held-out controls are imputed **separately** using the
existing Java21/Beagle runtime, exact marker mask and maps.

**Consistent reference policy:** exclude all 307 pure-training identities
plus the original 159 held-out/exclusion identities from the imputation
reference. Thus no selected training ancestor can be copied directly.
The five validation controls are re-imputed under this expanded policy,
rather than comparing training and validation generated with different
reference exclusions. Their old results remain preserved and are not used
as the primary matched-reference comparison.

Frozen LR base and CRF posterior means give 14 features per individual,
aggregated with the existing primary marker weights across all 22 autosomes.
The 307 clean pure-training individuals and 492 separate clean training
descendants are retained. Their original cached training base profiles are
checksum-verified and reused, avoiding base or local-model retraining.

Two fixed, inexpensive global heads isolate the effect of matched data:

| Head | Training policy |
| --- | --- |
| `global_ridge_clean` | Clean-only control: 50% pure / 50% pedigree |
| `global_ridge_matched` | 37.5% clean pure / 37.5% clean pedigree / 25% genuinely imputed training |

Both use the previous Ridge settings, `alpha=0.001` and no intercept.
Pure weights are seven-source balanced, pedigree weights family balanced,
and the seven imputed training individuals equally weighted. Prediction
uses the existing simplex projection helper. No hyperparameter selection,
per-person correction or fitting to validation cases is performed.

The unchanged pure, original-pedigree and real regression gates are
reported alongside the five newly imputed validation controls. The frozen
CRF itself is scored on these **same** inputs as a matched baseline.
Original clean frozen-model predictions must reproduce within `1e-10`.
Reloaded global heads must reproduce every held-out prediction before
completion is accepted. The local ancestry outputs, production models
and PGS calibration remain untouched; no automatic promotion occurs.

```bash
bash run_matched_global_training.sh --prepare
bash run_matched_global_training.sh --smoke
bash run_matched_global_training.sh --run
```

Preparation requires a new output directory, default
`results_v3_matched_global_01/`. The smoke stage checks chromosome22;
the full runner reuses completed checkpoints and runs two chromosomes
concurrently in fresh network-disabled containers. Original data/models
and imputation references are read-only; outputs are isolated. Completed
chromosomes are hash-verified on restart. Incomplete chromosomes require
inspection rather than silently overwriting partial results.

This is 264 singleton imputations (12 individuals x 22 autosomes), not a
broad model sweep. The expensive step is generating realistic training and
matched validation inputs; both global fits run only after full-genome
features are complete. This first small cohort cannot establish broad
population robustness. It uses one observed-array mask and complete Beagle
calls, not production GP/DS-QC missingness. Previously inspected controls
are development validation, not independent final certification.

The initial full run stopped on chromosome8 after a native JVM SIGSEGV
in compiled `java.text.DigitList.set`, during the first validation control's
Beagle execution. Chromosomes1-7 and22 were already complete and are reused.
The partial chromosome8 directory and crash log are preserved under
`failed_chr8_java_sigsegv_01/`; its seven completed training imputations are
repeated because the original implementation commits whole chromosomes,
not individual samples. Recovery retries the same inputs, image, seed and
thread count in a fresh container; it does not change model parameters or
treat the runtime crash as a performance result. See `runtime_failure_01.json`
and `full_run_resume_01.log` in the isolated results directory.

### Completed result: modest array gain, unacceptable matched-head trade-off

All 264 unique singleton imputations and all 22 feature checkpoints are
complete. Both global heads were fitted and independently verified:

| Model | Pure mean TV | Original pedigree mean TV | Matched array-imputed mean TV |
| --- | ---: | ---: | ---: |
| Frozen LR+CRF | 0.032818 | 0.027295 | 0.243200 |
| Clean-only global Ridge | 0.027077 | 0.025121 | 0.232300 |
| Matched-data global Ridge | 0.046018 | 0.037974 | 0.229262 |

These array results use the new shared exclusion policy; compare against
`0.243200`, not the historical `0.251716` obtained with a different reference.
Adding imputed training improves array mean TV by only **0.003038** beyond
the clean-only head, while worsening pure and original-pedigree mean TV.
In pure validation, original EAS/SEA group errors increase from
`0.004439/0.012348` in the clean head to `0.047841/0.064651`; GME increases
from `0.073390` to `0.130399`. The clean-only head improves aggregate pure
and pedigree scores but still fails the real regression gate.

The matched head predicts only **14.33% EUR** and **11.07% OCE** for the real
case and fails both its EUR and combined off-target constraints, plus the
pure no-baseline-regression gate. None of the three rows passes all gates;
no promotion or production changes occur.

Independent verification checked **4,165 input/output hashes**, all
264 singleton outputs, all supplied observed genotypes, and all model-marker
genotype/switch metrics. Both global heads were independently refitted
from training-only features and targets. Marker-weighted global features
and validation TV were independently recomputed; serialized predictions
agree within `1.12e-16`. All seven successfully imputed training genotypes
from the archived chromosome8 attempt exactly match their repeated outputs.
The original frozen clean-validation predictions reproduce within
`8.89e-16`. The JVM recovery did not change scientific parameters.

See `results_v3_matched_global_01/completion.json`, `comparison.tsv`,
`interpretation.json` and `verification/final_verification.json`.
Retain the clean global head as a development baseline, but reject this
matched-head recipe as a deployment fix. These small cohorts do not rule
out matched training in general. Before a broad hyperparameter sweep,
the next suggested model-family test is a phase-insensitive genotype/dosage
global model against the frozen haplotype models, reusing the verified
imputed training and evaluation inputs. No further experiment has started.

## Small phase-robustness smoother experiment (2026-10-16)

The matched array-imputation controls exposed a phase-related transfer gap:
mean TV increased from `0.0189823` to `0.2517156`; restoring identifiable
true heterozygote phase reduced it to `0.0497411`, without changing inferred
diploid genotype calls. The next experiment is **two fixed smoother-only
training arms**, not another model-family sweep.

Both reuse the complete frozen `gnomix_fast_050` LR base classifiers. Each
trains a fresh upstream CRF on the original 307 pure-training individuals
and 492 training descendants, plus one equally weighted phase-perturbed
view of each individual. The 234 training pedigree founders are checked
against the ancestry reference and evaluation-founder pools. No evaluation
controls, evaluation donors, or real sample enter training.

| Arm | Mean genetic phase-block length |
| --- | --- |
| `phase_blocks_050` | 0.5 cM |
| `phase_blocks_200` | 2.0 cM |

These are predeclared screening scales, **not error-rate estimates fitted
to the five validation controls**. A deterministic continuous-distance
two-state phase process exchanges paired haplotype alleles and their
marker-level ancestry labels together. It exactly preserves each diploid
genotype and both ancestry copies. Window labels are reduced from the
switched marker truth using the existing modal-label rule. Original
training data and upstream CRF hyperparameters are retained; doubling
the sequence count does change the effective data/regularization balance.
There is no additional base fitting or calibration.

This first experiment addresses phase robustness only. It does not use
matched-imputed training genotypes, reproduce genotype-call errors, or
claim that the synthetic switch process models actual real-person errors.
The original pure, original pedigree, and real regression tracks remain
unchanged. The five previously inspected singleton array-imputed controls
are scored separately as a fourth development-validation track. They are
not five independent founder replicates. Clean and array-imputed mean TV
are compared against the original frozen CRF, in addition to all unchanged
frozen-v1 acceptance gates; **no automatic production promotion** occurs.

Preparation and chromosome tasks run in fresh network-disabled containers,
with the reference, original models, and completed imputation pilot mounted
read-only. Two chromosome jobs run concurrently; each fits the two arms in
separate containers. Every serialized model must reproduce the full held-out
posterior tensors. Chromosome checkpoints verify input and output hashes
before reuse; incomplete jobs fail explicitly rather than silently restarting.
An initial chromosome-22 evaluation failed because historical pilot hashes
used `/pilot` paths. Its partial fit and log are preserved separately; mounting
the immutable pilot at its historical path resolved the preflight without
changing any model or input data.

```bash
bash run_phase_robustness_screen.sh --prepare
bash run_phase_robustness_screen.sh --smoke
bash run_phase_robustness_screen.sh --run
```

After a normal interruption, reuse `--run` to verify completed chromosome
checkpoints. An interrupted, incomplete chromosome requires inspection
before retry; it is not a resumable internal CRF optimizer. Do not launch
a second scheduler while one is active. Results are isolated under
`results_v3_phase_robustness_01/`. `completion.json` is written only after
both arms cover all 22 autosomes and all 197,260 shared model markers.
The existing SVM jobs remain excluded and are never restarted by this runner.

### Completed result: reject both phase-augmentation arms

Both arms completed all 22 autosomes. Neither improved the array-imputed
controls, and both substantially worsened clean pure and pedigree validation:

| Model | Pure mean TV | Original pedigree mean TV | Array-imputed mean TV |
| --- | ---: | ---: | ---: |
| Original frozen LR+CRF/0.5 | 0.032818 | 0.027295 | 0.251716 |
| Phase blocks, 0.5 cM | 0.228694 | 0.159253 | 0.302584 |
| Phase blocks, 2.0 cM | 0.153767 | 0.105509 | 0.259588 |

The real regression sample's EUR estimate increases to 24.49%/25.27%,
but both still fail its EUR and combined off-target constraints. Moving
that one prediction closer to its expected range is not an accuracy win:
all three known-truth means worsen. Neither model is promoted.

An independent audit checked **469 distinct input/output hashes**, all 66
profile/model checkpoints, all 44 fitted chromosome models, and every
marker-weighted full-genome prediction. Pure, pedigree and array-control
TV were independently recomputed. Serialized evaluation posterior error
was exactly zero. As a training-path sanity check, clean-only CRFs were
refit on chromosomes 1 and 22 using the new cached training profiles:
both reproduced the original CRF's full training and held-out posterior
tensors exactly. This supports interpreting the deterioration as an
outcome of these augmentation recipes, not a changed clean-training path.

See `results_v3_phase_robustness_01/completion.json`, `comparison.tsv`,
`interpretation.json` and `verification/` for the preserved evidence.
This rejects these two synthetic phase-switch recipes, not phase-robust
training in general. The next recommended bounded experiment is to use
actual array-masked, independently imputed views of a small balanced
subset of the **separate training pedigrees**, excluding their ancestor
donors from imputation references. Retain clean training, all validation
tracks and unchanged gates. Do not start another synthetic-rate sweep
or tune toward the inspected five controls or single real case.

## Twelve-candidate FLARE/Gnomix screen (2026-10-14)

The approved next screen contains **six source-anchored FLARE candidates and
six Gnomix candidates**, with ADMIXTURE explicitly excluded. It reports seven
broad ancestries: `EUR,EAS_SEA,SAS,GME,NAM,OCE,SSA`. Original EAS/SEA labels and
fine populations remain in cohort metadata and pure-validation diagnostics.
This is a coarser reporting task, not evidence of equivalent PGS calibration;
no production labels, `pgsc_calc` ancestry adjustment, or score calibration
are changed.

| Family | Configurations |
| --- | --- |
| FLARE direct | Seven labelled reference groups, fixed versus training-only EM refinement |
| FLARE fine | Seven labelled copying profiles with fine reference panels, fixed versus refined |
| FLARE subclusters | Two within-source kmeans components per ancestry, fixed versus refined |
| Gnomix fast | Logistic-regression base + CRF, nominal 0.2/0.5 cM windows |
| Gnomix default | Logistic-regression base + XGBoost, nominal 0.2/0.5 cM windows |
| Gnomix best | Random-string-kernel SVM base + XGBoost, nominal 0.2/0.5 cM windows |

FLARE copying profiles use only labelled pure-training donors. Its subcluster
initialization never clusters across source labels. Refinement uses the
307 pure-training genomes and 492 separate training descendants with known
ancestry priors: pure donors cannot acquire another source label during
training. Evaluation receives no truth priors. The refined copying/mismatch
matrices and haplotype-switch rates are retained, but both policies use
`T=10` and equal broad ancestry priors for evaluation, rather than importing
the training cohort's unrepresentative global ancestry frequencies.

Gnomix is pinned to upstream commit
`cd15f65eadfca8e2bbf209eb69a82978a200c5fd`. The actual upstream models are used;
the implementation is not an independent approximation of their architectures.
Reference donors train the base classifier; disjoint pure-training donors and
existing Ped-sim haplotypes train the smoother. No new mosaics, calibration,
or post-smoother base retraining are used. The base and smoother use upstream
unweighted training; these fits do **not** inherit the previous decoder
screen's 75/25 pure/pedigree weighting. Evaluation remains source/scenario/
combination/family balanced. A CRF output-axis adapter maps its first-seen
class order to the fixed source order without changing fitted probabilities.
Reloaded chromosome models must reproduce complete evaluation posterior
tensors before their output is accepted.

The upstream implementation translates nominal cM window sizes into
**fixed-SNP windows**, so their local genetic spans are not exactly constant.
The XGBoost smoother sees 75 windows at nominal 0.2 cM and 31 at nominal
0.5 cM, giving approximately comparable genetic context. The CRF uses the
whole chromosome and does not use this finite context setting. All other
settings are held fixed. These settings are screening choices, not claimed
optimal published hyperparameters.

Gnomix's licence is academic/non-commercial only; the user confirmed this
experiment is academic non-commercial research. The upstream notice is
preserved in `third_party/gnomix`. This does not authorize commercial SNPster
deployment or redistribution under a different licence.

### Common-marker and genotype-integrity controls

Preflight reused the existing 108-family/492-descendant training genomes and
27-family/123-descendant development genomes, with 307 pure training and
159 pure validation donors. All 3,447 actual reference donors have authoritative
pure-source labels. Affinity-only records in the larger role audit are not
silently assigned to new reference categories.

Two compatibility problems were caught before full experiments:

* Some legacy development VCFs contain different alternate records at the
  same position. For example, chr22:17082284 has both C>G and C>T records.
  FLARE skipped the matching C>T record in the compatibility test when the
  other alternate came first. Exact REF/ALT canonicalization retains the
  correct phased genotype, instead of dropping the position or changing its
  alleles. A regression test checks this second-record case.
* Eighteen reference SNPs lie outside the training pedigree map spans.
  The real-user regression VCF also lacks approximately 38,000 reference
  SNPs. Gnomix's upstream missing value is the numeric code 2, not genetic
  imputation; LR treats it as an ordinary feature and kernels treat it as an
  allele mismatch. Padding was therefore not used as an unverified shortcut.

The user explicitly approved a shared **197,260-SNP** panel for this screen.
Reference-only MAF/MAC filters, training map spans, and the real sample's SNP
identities define this panel; evaluation ancestry labels and allele frequencies
do not select markers. Every canonical training/evaluation VCF must have the
exact same ordered CHROM/POS/REF/ALT records and unchanged sample order.
Haplotype segment truth is reconstructed from existing `--bp_all` records,
and its physical-span totals must reproduce every saved pedigree truth vector.

The cached seven-source frozen-v1 baseline (pure TV 0.03969, pedigree TV
0.09762) remains diagnostic. An additional **matched frozen-v1 inference**
uses the shared panel and historical three-decimal decoder inputs, without
retraining, so input normalization is not mistaken for a model improvement.
Candidates are compared with that matched baseline. Marker-weighted posterior
means are primary, matching historical FLARE's global estimand. Physical-span
midpoint-weighted probabilities are separately exported, not mixed into the
cached comparison.

### Execution and acceptance

The new implementation is `pedigree_benchmark/scripts/screen_model_families.py`,
with the Gnomix adapter in `model_family_gnomix.py`. The current isolated
output is `results_v3_model_family_screen_05/`. Runs 01-04 preserve failed
preflight/compatibility evidence; none is a completed leaderboard. Existing results/models
are not overwritten.

```bash
export MODEL_FAMILY_SCREEN_OUTDIR=results_v3_model_family_screen_05
# --prepare refuses an existing output; these setup stages are one-shot.
bash run_model_family_screen.sh --prepare
bash run_model_family_screen.sh --initialize
bash run_model_family_screen.sh --matched-baseline
bash run_model_family_experiments.sh
```

The experiment scheduler first runs chromosome-22 compatibility checks,
then the complete 22-autosome candidates. Smoke outputs are not an accuracy
leaderboard. Every candidate/stage uses a fresh, network-disabled container
with read-only project/dependency inputs and isolated writable results. At
most two experiments run concurrently (12 CPU/80 GiB combined limits).
Failures stop the scheduler explicitly and its known active containers are
stopped; incomplete runs cannot produce the final leaderboard.

### Ten completed candidates; SVM stopped for runtime (2026-10-15)

The user approved moving forward with the ten completed configurations and
stopping `gnomix_best_020` and `gnomix_best_050`. Both containers have stopped.
Their saved chromosome-1 checkpoints were verified and retained, alongside
all completed predictions/models. These are explicit runtime exclusions,
not failed accuracy tests or completed full-genome SVM evaluations.
`results_v3_model_family_screen_05/screen_disposition.json` records the
decision. The twelve-candidate `completion.json` is intentionally absent.

| Configuration | Pure mean TV | Pedigree mean TV |
| --- | ---: | ---: |
| Matched frozen v1 | 0.045825 | 0.098466 |
| FLARE direct fixed | 0.183109 | 0.146779 |
| FLARE direct refined | 0.247144 | 0.229727 |
| FLARE fine fixed | 0.184830 | 0.177576 |
| FLARE fine refined | 0.254697 | 0.277748 |
| FLARE subclusters fixed | 0.190501 | 0.166087 |
| FLARE subclusters refined | 0.240733 | 0.216705 |
| Gnomix LR+CRF, 0.2 cM | 0.036111 | 0.028910 |
| Gnomix LR+CRF, 0.5 cM | **0.032818** | **0.027295** |
| Gnomix LR+XGBoost, 0.2 cM | 0.248103 | 0.165286 |
| Gnomix LR+XGBoost, 0.5 cM | 0.174281 | 0.115469 |

TV is half the L1 ancestry-proportion error, not classification accuracy.
All configurations use the same seven-source reporting task and canonical
197,260-marker panel. Both LR+CRF models pass the evaluable validation gates
but fail the unchanged real-case gate. At 0.5 cM the real sample is predicted
as 12.34% EUR and 70.41% EAS_SEA, outside the respective 35-65% bounds.
The baseline also fails the real case (30.97% EUR and 23.76% SAS).
No candidate is promoted. Minor 2-5% component recall remains unevaluable:
only 18 supported components, below the unchanged minimum of 20.

Pure-source averages also need scrutiny: LR+CRF/0.5 improves the overall
pure mean but worsens GME mean TV from 0.0581 to 0.1068 (14 individuals).
OCE and NAM have only three and four pure validation individuals.
The next priority is a targeted real-sample transfer/QC and difficult
population audit, not optimizing to the single real sample or launching
another broad hyperparameter sweep. Previously inspected development
cohorts and one real case do not establish independent production validity.

For a future implementation, consider Nextflow tasks keyed by candidate
and chromosome, each in a fresh pinned container, followed by a strict
complete-chromosome aggregation task. Bound concurrency by available CPU
and memory; preserve input/configuration hashes, deterministic seeds,
marker-weighted primary and physical-weighted secondary estimands, and
pure/pedigree/real evaluation separation. Nextflow caching alone is not a
substitute for these checks. No Nextflow migration is implemented here.

### Real-sample transfer audit (2026-10-15)

`results_v3_real_transfer_audit_03/` contains the completed, all-22-autosome
audit of frozen Gnomix LR+CRF/0.5. Its `audit.json` records input/model/output
hashes and `interpretation.json` separates findings from unresolved causes.
No models were trained, production files changed, or acceptance gates relaxed.
The original screen was mounted read-only in a fresh network-disabled container.

**Input handling passed:** all 197,260 canonical real-sample phased genotypes
match both the original concatenated VCF and the 22 surviving upstream QC
VCFs. The evaluation cohorts have no missing calls on this panel. Every
saved ancestry prediction reproduces to maximum absolute error below
`5e-16`; swapping complete real-sample haplotypes leaves the paired result
unchanged. This does not establish concordance with the unavailable raw
array calls or measure true phasing accuracy.

**A phase-independent diagnostic contradicts the very low EUR prediction:**

| Real-sample diagnostic | Sites | EUR | EAS+SEA |
| --- | ---: | ---: | ---: |
| Reference-frequency projection, all sites | 197,260 | 47.3% | 45.2% |
| Projection, observed sites only | 6,205 | 43.6% | 51.4% |
| Projection, observed sites thinned to 100 kb bins | 5,359 | 42.3% | 51.3% |
| Projection, imputed sites only | 191,055 | 47.5% | 44.9% |
| Saved LR+CRF prediction | 197,260 | 12.3% | 70.4% |

The frequency projection uses reference-only ancestry allele frequencies
and constrained least squares on diploid allele counts, not haplotypes.
It is a diagnostic, **not** a new certified model or corrected real-sample
truth: its pure/pedigree mean TVs on all markers are 0.1154/0.0876, worse
than LR+CRF on these validation cohorts. The observed/imputed comparison
does not suggest a simple wholesale loss of European diploid signal during
imputation, but cannot exclude genotype errors.

**Haplotype transfer/smoothing is the main remaining lead.** The real
sample's unsmoothed LR means are 28.0% EUR, 36.0% EAS+SEA and 22.1% SAS;
the CRF shifts these to 12.3%, 70.4% and 7.7%. Five existing 50/50
EUR/Asian pedigree controls have mean TV 0.0190 with the original phase.
Randomly swapping alleles at markers, while exactly preserving their
diploid genotypes, increases their mean TV to 0.7211. This establishes
severe phase sensitivity, **not** actual switch errors in the real sample.

Simply shortening CRF sequences without retraining did not fix transfer:

| Frozen-weight CRF control | Pure mean TV | Pedigree mean TV | Real EUR |
| --- | ---: | ---: | ---: |
| Original whole chromosome | 0.0328 | 0.0273 | 12.3% |
| 100-window sequences | 0.0793 | 0.0545 | 13.1% |
| 20-window sequences | 0.2448 | 0.1645 | 17.5% |

These are inference-only mechanism controls, not accepted additional
candidates. Raw LR probabilities are not a validated replacement either.
Relevant Filipino groups are present (30 reference, five pure-training,
three pure-validation donors), but coverage is sparse and not an
independent population-disjoint test.

The singleton imputation export's DR2 is almost a genotype indicator:
mean DR2 is 0.00169 for homozygous sites and 0.98990 for heterozygous sites.
Low DR2 alone therefore does not prove poor individual calls; do not add
a DR2 threshold as an ancestry fix. The audit also retains the inherited
two-decimal GP rounding (some triplets sum to 0.99/1.01); GP was not used
to change genotypes or fit the ancestry models.

**Next justified experiment:** pass a small set of existing held-out
pedigrees through the real array-mask/imputation/phasing workflow, with an
audited donor-disjoint imputation reference. Compare recovered diploid
genotypes and switch errors with the known original haplotypes, then
score the frozen LR+CRF model. This distinguishes upstream phase/genotype
degradation from estimator/domain mismatch before another training sweep.
No such re-imputation/retraining has been launched by this audit.

To reproduce the audit in a new isolated directory:

```bash
REAL_TRANSFER_AUDIT_OUTDIR=results_v3_real_transfer_audit_04 \
    bash run_real_transfer_audit.sh
```

`results_v3_real_transfer_audit_01/` preserves the initial GP-rounding
preflight failure. Run 02 preserves the completed first audit and exact
script snapshot; run 03 adds the bounded CRF sequence controls and correct
independent thinning of observed markers.

### Matched array/imputation pilot (2026-10-16)

**Completed and independently verified at approximately 13:34 CEST.** All
110 singleton imputation tasks finished across 22 autosomes. The pilot
used 148,014 of 148,019 observed array sites; the five unavailable sites
were counted explicitly. All supplied observed genotypes were preserved.
The original haplotypes reproduced exactly and their frozen genome-wide
predictions reproduced within `5.6e-16`. Final verification checked 2,299
input/reference/model/output hashes; five targeted tests passed.
No experiment containers remain running.

| Five-control stage | Mean TV | Mean EUR | Mean EAS+SEA |
| --- | ---: | ---: | ---: |
| Original known haplotypes | 0.01898 | 48.13% | 50.15% |
| Array masked, then statistically phased/imputed | **0.25172** | **24.83%** | **61.23%** |
| Restore known phase at correctly imputed true heterozygotes | 0.04974 | 45.11% | 50.16% |

The full-panel genotype error rate was 4.88% pooled across controls
(per-control range 2.35-6.23%). The pooled switch-error rate was 4.31%
between adjacent correctly imputed truth heterozygotes, with per-control
rates from 0.15% to 6.56%. Sparse observed-only switch rates use much wider
marker spacing and must not be compared directly with full-panel rates.

Restoring only identifiable phase, without changing imputed diploid
genotypes, reversed approximately 87% of the increase in mean TV. The
worst controls reproduce the very low EUR estimates seen in the real
sample. This is **controlled evidence of a phase-related transfer failure
in the frozen model**, not proof of the real individual's phase errors
or exact ancestry. Perfect-phase pedigree scores substantially overstate
this model's robustness to the matched statistically phased array domain.

The next justified work is a small phase-robustness experiment using only
separate training donors/pedigrees, with original pure and pedigree
validation retained and a separate array-imputed/phase-error validation
track. Do not train on these five controls, tune to the single real case,
relax existing gates, or promote the current model. `summary.json`
contains machine-readable aggregate results, `predictions.tsv` contains
each control/stage, and `interpretation.json` records scope and limitations
under `results_v3_matched_imputation_02/`.

The approved follow-up reuses the five existing 50/50 EUR/Asian controls.
`pilot_matched_imputation.py` traces the original Ped-sim breakpoints into
the same ten founder donors' dense reference haplotypes to recover their
real-array SNP genotypes. Every chromosome must reproduce the original
197,260-model-marker phased genotypes exactly before proceeding. This
extends marker coverage, not the family design or original inheritance.

Each control is stripped to **unphased diploid genotypes** at exact
CHROM/POS/REF/ALT matches to the real sample's observed sites. Each is then
imputed **separately**, not jointly with related controls. The full dense
imputation reference excludes all development founders and pure-validation
individuals present in it: 120 of 4,091 reference individuals on chr22,
leaving 3,971. This is identity disjointness, not proof of kinship
disjointness. Missing array sites outside original breakpoint spans or
absent from the reference are counted explicitly, never imputed into truth.

The runtime uses the Beagle 5.5 jar extracted from the pinned installed
`imputer` image, a fixed seed of 20261013, six threads and the production
GRCh38 maps. The exact historical Beagle binary/seed for the real sample
cannot be verified from the surviving QC VCFs; this is a matched workflow
control, not a byte-identical reconstruction of that historical job.
The installed jar requires Java 21, whereas the previous ancestry image
has Java 17. `Dockerfile.matched-imputation` adds a digest-pinned Java 21
runtime without changing the frozen Gnomix Python dependencies. Run 01
preserves the caught class-version preflight failure; run 02 uses the
corrected runtime.

```bash
docker build -f Dockerfile.matched-imputation \
    -t snpster-ancestry-benchmark:matched-imputation .
export MATCHED_IMPUTATION_OUTDIR=results_v3_matched_imputation_02
bash run_matched_imputation_pilot.sh --prepare
bash run_matched_imputation_pilot.sh --smoke
bash run_matched_imputation_pilot.sh --run
```

Preparation requires a new result directory. The five-control chromosome22
preflight must complete before the all-autosome run. The runner uses fresh,
network-disabled containers, at most two chromosomes simultaneously
(12 CPUs/80 GiB limits), with read-only source genomes/references/models.
Completed chromosomes are checksum-verified and reused by `--run`;
incomplete ones fail explicitly and require inspection, not automatic
overwriting. This bounded runner is not a Nextflow migration; its
chromosome task boundaries are suitable for a future Nextflow workflow.

The pilot records:

* Exact preservation of observed diploid genotypes, including typed sites
  outside the frozen model panel.
* Per-sample genotype errors and dosage squared error, separately on all,
  observed and imputed model-marker sites.
* Switch errors between consecutive correctly imputed true heterozygotes
  within a chromosome, with explicit denominators; genotype-error sites
  are excluded. Hamming error allows a global chromosome haplotype swap.
* Frozen LR+CRF predictions before and after the array/imputation workflow.
* An oracle diagnostic restoring known phase at correctly imputed true
  heterozygotes while preserving inferred diploid genotypes. False
  heterozygotes and genotype-error sites do not have identified phase;
  this is not a complete separation of genotype and phase error.

The strict genome-wide summary requires all 22 chromosomes and must
reproduce the original frozen five-control predictions. No ancestry model
is retrained or promoted. These controls share some donors, are not five
independent family replicates, and are not necessarily EUR/Filipino;
results will diagnose transfer mechanisms, not certify real-user ancestry.
Dense temporary imputation outputs and reference subsets are removed only
after successful scoring; exact truth/array/model-panel genotypes, commands,
logs, hashes, and original sources are retained for reproducibility.

### Historical pause and resume instructions

**After the 2026-10-15 runtime exclusion, do not run the resume command below
unless intentionally restarting the two SVM candidates.** The original
scheduler still expects all twelve; the ten-model assessment does not
silently redefine its completion requirements.

Run 05 was stopped at the user's request on 2026-10-14. All six full FLARE
fits, Gnomix fast/0.5, all twelve compatibility checks, and chromosomes 1-10
of Gnomix fast/0.2 are saved. Chromosome 11 had begun preparation during
container shutdown but has no completed model checkpoint. No saved fit is
discarded. The partial screen is not an accepted final leaderboard.

From this project directory, resume without repeating completed fits:

```bash
export MODEL_FAMILY_SCREEN_OUTDIR=results_v3_model_family_screen_05
bash run_model_family_experiments.sh --resume
```

Optionally check readiness without training or changing saved predictions:

```bash
bash run_model_family_experiments.sh --verify-resume
```

Resume verifies the pinned Gnomix image, input hashes, configuration,
dependencies, original implementation, saved model checksums, marker/window
identities, and completed predictions. The original adapter and scheduler
are archived under `results_v3_model_family_screen_05/revisions/`. Legacy
checkpoints are accepted only when their original adapter checksum matches
and its fitting/evaluation code is unchanged. Saved chromosomes are reloaded
and re-evaluated to reconstruct genome-wide totals, not retrained. Missing
chromosomes and pending configurations run in fresh containers.

For subsequent runs, request a cooperative pause from another terminal:

```bash
export MODEL_FAMILY_SCREEN_OUTDIR=results_v3_model_family_screen_05
bash run_model_family_experiments.sh --pause
```

Gnomix finishes, saves, and evaluates each active chromosome before exiting;
the scheduler starts no new candidate. Exit status 75 means an intentional
checkpoint pause, not a completed screen. FLARE candidates finish their
current candidate before pausing. `--resume` clears the pause request.
Once the scheduler exits, no experiment container or its RAM needs to remain
running; saved work survives a host restart. The resume command runs in the
foreground: keep that terminal/session open, or request a separately managed
background launch. Do not rerun the one-shot preparation commands or remove
the result directory.

Pure-validation and pedigree mean TV are both primary: candidates must not
regress either relative to matched frozen v1. Existing numeric tail/off-target/
minor-support and real-user constraints are retained, but seven-source results
are explicitly not the original eight-source selection objective. Insufficient
minor-component support stays unevaluable. No automatic production promotion
is performed. This remains previously inspected development data, not a new
population-disjoint final test, phenotype-based PGS validation, or production
microarray/imputation certification.

Shared-panel preparation and FLARE initialization completed. The matched
frozen-v1 control has pure mean TV **0.045825** and pedigree mean TV
**0.098466**; it still fails the unchanged real-user constraint. These are
the new comparison values, not replacements for the cached historical
seven-source metrics on the old marker sets.

Validation passed **26 targeted tests** (9 common-screen and 17 Gnomix
adapter tests), plus 6 existing precision regression tests. A real
Gnomix fast/0.2 smoke fit on chr22 completed in 174 seconds, using all 3,447
reference donors, 1,913 shared chr22 markers and 283 evaluation samples.
Every evaluation posterior reproduced exactly after model reload; an
independent fresh-container audit checked all 76 input hashes and both output
checksums. This is compatibility evidence, not full-genome performance.
Window-count extrapolation suggests roughly 2-4 hours for that full candidate;
the XGBoost/kernel configurations need their own timing measurements.

Run04's four direct/fine FLARE compatibility tests passed, but its subcluster
models were rejected by FLARE because sklearn's subtraction-based kmeans
centres contained tiny negative probabilities (approximately `-7e-18`).
The initializer now computes empirical means of assigned nonnegative copying
profiles, preserving exact zero coordinates, and model readers/writers enforce
FLARE's strict `[0,1]` bounds. This is a numeric parameter-construction fix,
not a change to cluster assignments or an accuracy fallback.

Run05 verified all previous input hashes against the archived original
implementation, reused the canonical genomes/truth, matched frozen baseline,
raw copying profiles, and five completed compatibility tests, and rebuilt
only the cluster parameters. Previously valid direct/fine initial models
remain byte-identical. Both repaired subcluster compatibility tests and
the second fast Gnomix test passed, giving **eight verified configurations**.
The scheduler is running the remaining XGBoost/kernel compatibility tests;
full-autosome fits follow only after all twelve pass. Recovery evidence is
in `results_v3_model_family_screen_05/recovery_audit.json`.

The actual run05 recovery/continuation commands were:

```bash
MODEL_FAMILY_SCREEN_OUTDIR=results_v3_model_family_screen_05 \
  bash run_model_family_screen.sh --recover results_v3_model_family_screen_04
set -o pipefail
MODEL_FAMILY_SCREEN_OUTDIR=results_v3_model_family_screen_05 \
  bash run_model_family_experiments.sh --resume-verified-smokes \
  2>&1 | tee results_v3_model_family_screen_05/experiment_scheduler.log
```

## Ten-model frozen pilot and work cleanup (2026-10-13)

The old `work/` directory was removed **before pilot validation/execution**,
freeing approximately 729 GiB. Published historical results remain intact.
Before deletion, all 80 saved candidate pairs were checked for complete
artifacts and absence of symlinks into work. The exact 6,001-entry stage-2
manifest was repointed to published `stage2_candidates/` directories after
byte-for-byte checks of metrics, candidate configuration and predictions.
The durable manifest is now `selection/recovery_stage2_result_dirs.txt`.
All 572 simulation tasks' founder assignments, logs and command lines were
preserved under `pedigree/simulation_provenance/`.

The initial corrected pilot evaluates **10 frozen model/decoder pairs on
27 reduced-design families (three per class), across all 22 autosomes**.
`config/frozen_pilot_10.tsv` includes frozen v1, the historical pedigree winner,
one representative of each other shortlisted strategy family chosen by
reference-validation TV, and a separate spherical-GMM/HistGradientBoosting
control. There is one pair per representation. This is a small diagnostic
comparison, not an exhaustive model search or an independent final test.

The approved pilot completed all 66 chromosome/batch simulations, 30 FLARE
inference tasks and ten evaluations of 123 descendants. Output is
`results_v3_pedigree_pilot_10/`, with isolated work directory
`work_pedigree_pilot_10/`. Frozen v1 has the lowest pedigree mean TV (0.16318;
p95 0.48544). Unified selection exited nonzero because no candidate passed the
unchanged acceptance gates, not because simulation or inference failed. Its
diagnostic unified leaderboard remains in
`work_pedigree_pilot_10/5f/0d93ff2e1aff1e4f2b8590a91e37bc/`.
The smaller, redesigned cohort is not directly comparable to the historical
252-family benchmark.

Frozen mode imports the existing pedigree processes from `main.nf` using
`pedigree_pilot.nf`. It skips representation training, decoder/calibrator
fitting, Stage 1/2 selection and reference/real-user inference. Saved reference
predictions and real-user gate results are reused without recomputation.

Reusable maps and reference marker sites were archived in
`results_v3_sweep_100/pedigree/reusable_inputs/`. Reuse requires exact SHA-256
matches for the original source genetic map, reference VCF and every archived
map/sites file. **Old offspring genotypes, ancestry truth and pedigree features
are not reused**: the corrected founder pool changes donor assignments and the
cohort must be regenerated. Its new pedigree outputs can later be used for
additional frozen candidates on the same cohort.

```bash
FROZEN_CANDIDATES_MANIFEST=/workspace/config/frozen_pilot_10.tsv \
PEDIGREE_REUSE_DIR=/workspace/results_v3_sweep_100/pedigree/reusable_inputs \
RUN_OUTDIR=/workspace/results_v3_pedigree_pilot_10 \
RUN_WORKDIR=/workspace/work_pedigree_pilot_10 \
PEDIGREE_DESIGN=reduced PEDIGREE_FAMILIES=27 PEDIGREE_BATCH_SIZE=10 \
./run_in_container_v3.sh
```

Do not launch a second copy while the pilot is active. Frozen mode cannot be
combined with stage-2 recovery. Its input audit records the model/configuration
checksums and its trace/report/timeline are written inside the new output
directory, not over the historical run's reports. Pilot preparation passed
26 unit tests, shell checks, real artifact validation and founder preflight.

## Small frozen-v1 decoder screen (2026-10-13)

Run `bash run_decoder_screen.sh` on the Docker host. Each audit, decoder
experiment and summary starts a **fresh named container**, with two CPUs,
8 GiB memory, no network, read-only input mounts and a separate writable result
mount. No FLARE inference or genome simulation is repeated. An existing output
directory is never overwritten; set `DECODER_SCREEN_OUTDIR` to a new
`results_v3_decoder_screen_*` directory name for another screen.

The screen reuses the corrected pilot's v1 features and realized physical-span
truth, 307 pure-source training feature vectors and 159 pure-validation
vectors. Frozen artifact hashes are checked against the pilot input audit.
Feature schemas, authoritative labels, actual reference/mosaic donor exclusions,
sample IDs and reproduction of saved v1 predictions are validated again.
Every input and experiment script is SHA-256 recorded.

Five explicitly fixed experiments (no hyperparameter search) compare:

| Experiment | Training | Decoder |
|---|---|---|
| `ridge_pure` | pure sources | no-intercept Ridge, alpha 0.001 |
| `ridge_mixed` | pure + pedigrees | same Ridge |
| `stochastic_pure` | pure sources | nonnegative feature-to-source matrix, every row sums to one |
| `stochastic_mixed` | pure + pedigrees | same constrained matrix |
| `extratrees_mixed` | pure + pedigrees | 200 trees, minimum leaf size 3 |

Ridge/tree predictions use the existing simplex projection. The constrained
matrix preserves mixture additivity without a post-hoc calibrator; optimization
failures are explicit. Its mean-squared-error regularization is 0.00001.
Pure training gives equal total weight to each source. Mixed training gives
equal total weight to the pure and pedigree sets, with equal pedigree-family
weight within the latter. No calibrator is fitted in these five experiments.

Pedigree predictions are out-of-fold: hold out each of the 27 families and
**purge every training family sharing any founder donor**. This leaves 9-24
training families per fold. Splitting only by family or descendant ID would
leak repeated founders. Fold models and donor lists are saved. Pure-training
donors must be disjoint from all pedigree and pure-validation donors.

Some pure-validation donors are also pedigree founders. Reference predictions
therefore average only fold models that did **not** train on that validation
donor. This keeps all eight sources represented without scoring a donor on a
model trained on its descendants. These reference averages are diagnostics,
not predictions from the full-data development model. The known real-user
constraints are checked unchanged using that development model.

Outputs under `results_v3_decoder_screen_01/` include residual/family/donor
audits, the uncalibrated frozen baseline, expected-parent-feature diagnostics,
paired family error deltas, all fold models, development models and
`decoder_screen.tsv`. Donor associations do not establish causation. Expected
parent averages ignore realized recombination and are not production inputs.

**This is development screening, not independent final validation.** The pilot
has already been inspected, related founders recur across folds, and no model
is promoted or gate relaxed. Improvements need a new founder-disjoint cohort
and production microarray/imputation validation before deployment. A
full-data development model must never be evaluated on these same pedigrees
as if they were held out.

### Completed screen 01

The audit and all five experiments completed successfully in separate fresh
containers. All 135 fold models were checked for zero training/test donor
overlap, with 123 out-of-fold descendants and 159 donor-safe reference
predictions per experiment. The main regression suite passed 28 tests,
including seven new decoder-screen tests.

| Decoder | Descendant mean TV | Descendant p95 TV | Pure-validation mean TV |
|---|---:|---:|---:|
| frozen v1 baseline | 0.16318 | 0.48544 | 0.08785 |
| Ridge, pure | 0.16295 | 0.37093 | 0.23159 |
| Ridge, pure + pedigree | 0.17742 | 0.38275 | 0.24553 |
| stochastic matrix, pure | 0.18192 | 0.38857 | 0.26180 |
| stochastic matrix, pure + pedigree | 0.19193 | 0.40991 | 0.26879 |
| ExtraTrees, pure + pedigree | 0.30405 | 0.49888 | 0.14306 |

The pure-trained Ridge reduces the tail but does not materially improve mean
TV and worsens pure-source validation. It improves 11 families and worsens 16.
All five new development models fail the unchanged known real-user regression
constraints. Frozen v1 already failed that regression gate in the pilot.
**No model was promoted and no gates were changed.** The reference column uses
donor-safe fold averages, not the full-data model (see evaluation protocol above).

Residual diagnostics suggest a more specific next experiment:

- Collapsing EAS/SEA **for diagnosis only** reduces v1 mean TV from 0.16318 to
  0.09762: 0.06556 of mean TV disappears, without changing the production
  ontology or the actual leaderboard.
- Actual descendant features differ from expected-parent feature averages by
  mean TV 0.03722. Decoding those ideal feature averages still gives TV 0.15965
  against expected ancestry, so imperfect realized transmission is not enough
  to explain the residual decoder errors.
- In `fam000001`, actual and expected-parent feature vectors differ by only
  0.0005 TV, but both produce 0.5 ancestry TV. Its SEA founder `HGDP00720` is
  predicted as 100% EAS when inferred independently; its EUR founder
  `NA12344` is predicted as 93.8% EUR. The mixed descendant is predicted as
  100% SEA. This exposes both source confusion and non-additive decoding;
  it is not evidence to relabel or exclude the founder.
- Removing the frozen calibrator is not a solution: uncalibrated mean TV is
  0.17442 versus calibrated 0.16318.

These are negative screening results, not a limit on attainable performance.
The new training variants deliberately use pure sources and/or pedigrees, not
the original 1,008 many-donor mosaics, and change decoder/calibration choices.
They therefore do not isolate the causal effect of adding pedigree training.
A follow-up should preserve the original pure-plus-actual-mosaic training data
and add modest, family-balanced pedigree weight with a matched no-pedigree
control; separately inspect EAS/SEA population coverage and problematic pure
founder predictions. Such a follow-up has **not** been launched.

## Separate pedigree training cohort and matched screen (2026-10-13)

Following discussion, the next experiment replaces old mosaic training with
**pure-source individuals plus real Ped-sim pedigrees**. Mosaics are retained
only inside the untouched historical frozen-v1 baseline; no old mosaic
examples are used to fit the new decoders.

The 307 pure-training donors have no overlap with FLARE's reference, the
pure-validation set or the corrected pilot's founders. Their existing phased
`source_train.autosomes.vcf.gz` contains all 307 donors and 235,503 marker
records. It includes extra populations absent from the larger chromosome
reference files, so those files must not be silently substituted. Split the
panel-ready founder VCF by chromosome, then use the existing simulation,
founder-assignment verification, physical-span truth, panel projection and
FLARE inference processes. No consumer-array masking or imputation is done.

The deterministic full-design training plan has **108 families**, 12 per
class, **492 descendants**, and uses **234 unique pedigree donors** from the
307-donor pure-training pool. The pedigree donors include six OCE, nine NAM
and 35 SEA individuals. Reuse across training families is allowed, but within
a family donors remain unique. `trainfam*` node names prevent collisions
with the existing development cohort. A separate training-only audit is used;
the stricter held-out-founder selection and scoring safeguards are **not**
disabled or misrepresented.

Evaluation remains fixed on the previously inspected 27-family reduced-design
pilot (123 descendants), 159 pure-validation individuals and the unchanged
known real-user regression constraints. **Every training donor is disjoint
from every evaluation founder.** This is development evaluation, not an
independent final test. Training may use the full design while evaluation
stays reduced; their scores are never pooled as if they were one cohort.

Three fixed decoder architectures each have a matched pure-only control and
a pure-plus-pedigree arm:

- No-intercept Ridge, alpha 0.001.
- Row-simplex-constrained linear mapping, regularization 0.00001.
- ExtraTrees, 200 trees, minimum leaf size 10, all features available.

Within a pair, architecture, hyperparameters, seed, feature ordering and
calibration are identical. No calibrator is fitted in either arm. The mixed
arm gives pure sources 75% total training weight and pedigrees 25%, with
source-balanced pure weights and family-balanced pedigree weights. This is
six fixed fits, not a broad hyperparameter search. Frozen v1 is an additional
historical baseline, not a matched pure-only control.

Run on the Docker host:

```bash
bash run_pedigree_training_screen.sh
```

The launcher uses a **new named container for cohort generation, each of the
six decoder experiments, the baseline and the summary**. The cohort uses
20 CPUs/180 GiB; decoder containers use two CPUs/8 GiB with read-only inputs.
No network is enabled. Output is `results_v3_pedigree_training_screen_01/`;
work is `work_results_v3_pedigree_training_screen_01/`. Existing output/work
directories are rejected. Set `TRAINING_SCREEN_OUTDIR` to another isolated
`results_v3_pedigree_training_screen_*` name for a new experiment.

The cohort has 242 chromosome/batch simulations and 11 frozen-v1 inference
tasks. Only after successful completion are the decoders trained and scored.
The launcher records the image identity and source hashes and refuses to
score if experiment code changes during generation. SHA-bound static maps
and marker sites are reused. Historical results and previous screens are
not overwritten. Do not edit active experiment scripts or workflows.

Preparation passed **38 unit tests**, Pylance compatibility checks for the
optional family-name prefix, previews of both the new training workflow and
the existing frozen pilot, and a real chromosome-22 Ped-sim smoke test:
47 distinct donors, 52 assigned founder nodes, 94 total nodes, all nine
classes, and normalized realized truth.

The user approved persistent execution after the chat closes. The run
completed successfully: all 242 simulations, 11 genome merges and 11
frozen-v1 inference tasks finished; cohort generation took 10 minutes
32 seconds. The baseline and all six matched decoder fits then completed in
their separate containers, and `matched_decoder_screen.tsv` was generated.
All six saved models reproduced every development prediction on reload.
Input hashes remained unchanged, and coverage was verified for 492 training
descendants, 123 development descendants and 159 pure-validation individuals.

| Architecture | Pure-only pedigree mean TV | Pure + pedigree mean TV |
|---|---:|---:|
| Ridge | 0.16295 | 0.16350 |
| constrained linear mapping | 0.18192 | 0.18205 |
| ExtraTrees | 0.35049 | 0.27225 |

Frozen v1 remains 0.16318 mean pedigree TV, with pure-validation mean TV
0.08785. Ridge with pedigree training has p95 pedigree TV 0.37063, but its
pure-validation mean TV is 0.23364. Adding pedigrees improved the matched
ExtraTrees control by 0.07824 absolute TV, but did not bring it close to v1.
The two linear architectures changed very little with pedigree training.
These findings apply to the tested configurations, not to all possible
pedigree-trained models.

All six new models failed the unchanged known real-user regression constraints;
their balanced p95 TV also exceeded the existing 0.40 threshold. No model
was promoted and no acceptance threshold was relaxed. The results do not
establish a performance ceiling; targeted investigation of EAS/SEA source
coverage and feature/decoder discrimination remains warranted.

## Completed ten-representation diagnostic audit

[run_representation_audit.sh](run_representation_audit.sh) evaluates the ten
saved representations in separate, fresh, network-disabled containers, then
uses separate containers for summary and exact-feature/genotype controls.
It refuses existing output directories. To rerun without overwriting screen 01:

```bash
REPRESENTATION_AUDIT_OUTDIR=results_v3_representation_audit_02 \
  bash run_representation_audit.sh
```

Screen 01 completed in
[results_v3_representation_audit_01](results_v3_representation_audit_01).
It reused existing features for 307 pure-training donors, 159 pure-validation
donors and 123 descendants in the corrected 27-family development cohort.
No new genomes or FLARE inference were required. Every representation's
saved reference and pedigree predictions were reproduced to absolute
tolerance 1e-10. Reloaded source probes also reproduced their predictions;
input hashes and complete sample coverage were verified.

The audit includes source/population errors, founder errors and nearest
training vectors, feature-inheritance and decoder-mixture diagnostics, and
a fixed EAS/SEA logistic probe (C=10, balanced classes). The probe uses only
63 EAS and 39 SEA pure-training donors, with no validation-label fitting or
hyperparameter search. Population-held-out probes additionally exclude the
validation population from label training. The frozen representation itself
was **not** retrained without that population: this is not a fully
unseen-population representation benchmark.

Selected findings from [representation_diagnostics.tsv](results_v3_representation_audit_01/representation_diagnostics.tsv):

| Frozen representation | Pedigree mean TV | Pure-validation mean TV | SEA pure-validation mean TV |
|---|---:|---:|---:|
| v1 K18 | 0.16318 | 0.08785 | 0.33708 |
| unrefined diagonal K30 | 0.21199 | 0.05611 | 0.09976 |
| diagonal K36 | 0.21754 | 0.03329 | 0.08407 |

For v1, combining EAS/SEA **only as a diagnostic** reduces pedigree mean TV
to 0.09762. The production labels and gates remain unchanged. The v1 probe's
validation AUC is 0.93887, falling to 0.69436 with population-held-out label
training. Among the 28 pedigree-eligible EAS/SEA validation donors, including
only four SEA donors, balanced accuracy is 0.50000 and 0.43750 respectively.
The unrefined K30 population-held-out AUC is 0.88715. These results suggest
population-dependent discrimination, not a universal absence of signal;
the small eligible SEA subset limits generalization.

All ten frozen decoders show substantial non-additive mixture behavior.
For v1, decoding expected parental-average features gives mean TV 0.15965
against expected pedigree truth, versus 0.08344 when averaging the decoded
known founders; the mean consistency gap is 0.12704. These are diagnostic
comparisons, **not** deployable predictors, achievable-performance claims or
an additive causal error decomposition. Expected truth is distinct from
the realized recombination-weighted scoring truth.

### Exact-feature collisions and global-export precision

[audit_feature_collisions.py](pedigree_benchmark/scripts/audit_feature_collisions.py)
persists controls under [collisions](results_v3_representation_audit_01/collisions).
Without rounding the input vectors, v1 has two cross-source ambiguous groups
containing 37 of the 466 combined pure donors. All ten representations have
at least one such group, with different affected donors. HGDP00720 shares
its exact recorded v1 vector with both EAS and SEA training donors.
Chromosome-22 controls confirm these are not identical genotype copies:
896 and 813 unphased genotype disagreements across 2,522 aligned variants.
Any deterministic decoder restricted to an identical vector must give
identical predictions despite the differing labels. Collision-derived mean
TV lower bounds are finite-dataset diagnostics, not biological ceilings or
pedigree-performance limits.

Verified [public FLARE source](https://github.com/browning-lab/flare/blob/d6a6c8fd4479ebdb68529824acd11d2b6bbac4a1/src/admix/EstimatedGlobalAncProportions.java)
rounds global proportions to three decimal places. Global values already
average posterior probabilities across markers and haplotypes, including
with `probs=false`; `probs=true` controls local posterior output, not
global-file precision. Installed jar constants corroborate the formatter,
scale and accumulator; the jar hash and public-source identifiers are
recorded in [global_precision_evidence.json](results_v3_representation_audit_01/collisions/global_precision_evidence.json).
Identical recorded vectors do **not** prove identical unrounded posteriors.
The relative contribution of export rounding versus genuine posterior-summary
overlap remains unmeasured.

All **44 regression tests** passed after adding the audit controls. No model
was promoted, labels were not changed and acceptance gates were not relaxed.
These are development diagnostics, not independent final validation.
The matched high-precision global-export control below has since completed.
It does not establish that precision-aware retraining will improve prediction.

## Matched global-export precision control (2026-10-14)

[run_global_precision_screen.sh](run_global_precision_screen.sh) completed
three matched v1 inference arms: installed stock jar, an isolated jar with
the original exporter recompiled from verified public source, and an isolated
jar with only that exporter's output precision increased. Preparation, each
inference arm and summary ran in fresh, network-disabled containers. The
public Java source was downloaded before execution and its Git blob hash
verified. The existing host JDK was mounted read-only for compilation; no
dependencies were installed and the installed jar was not modified.

The screen uses 11 reference-disjoint pure targets: nine validation donors
(the Cambodian and three KHV pedigree founders, two EAS controls, one EUR
control and both eligible OCE founders), plus the two previously identified
feature-identical EAS/SEA training donors. These are selected diagnostic
controls, not a representative validation cohort. The frozen v1 model,
seed 20261002, reference, markers and inference parameters were held fixed.
No decoder was trained and no model was promoted.

Results are persisted in
[results_v3_global_precision_01](results_v3_global_precision_01), including
[summary.json](results_v3_global_precision_01/summary.json),
[sample_predictions.tsv](results_v3_global_precision_01/sample_predictions.tsv),
the three inference outputs, isolated jars, source/build inputs, command logs
and checksum provenance. Each jar changed exactly one class entry. The
precision variant preserves posterior accumulation and normalizes the
accumulated values in double precision when exporting, rather than selecting
a three-decimal lookup-table value.

### Verified outcome

- All three arms processed 235,453 markers and all 11 targets. Their local
  ancestry VCF data and sample-order hashes were identical.
- Stock and recompiled-source global outputs matched, and their frozen-decoder
  predictions matched to tolerance 1e-12. Stock normalized features reproduced
  the historical features to numerical precision (maximum difference 1.11e-16).
- Rounding the high-precision global output to three decimals reconstructs
  every stock global value. The maximum raw component difference was 0.0004964.
- Both tested exact collisions disappeared. HGDP00720's feature TV distance
  from the Korean EAS donor is 0.00039665; from the Vietnamese SEA donor it is
  0.00181922. Thus export rounding erased real posterior-summary differences
  in these cases, but their remaining separation is small.
- Passing precise features through the **unchanged, rounded-feature-trained**
  decoder did not solve the Cambodian error: it remained 100% EAS. For the
  nine selected validation controls, mean TV worsened from 0.26020 to 0.29971.
  Across all 11 controls, it changed from 0.30380 to 0.32900. These selected
  pure-control scores must not be compared with the pedigree leaderboard.
- The largest per-sample prediction TV change was 0.30309, despite maximum
  normalized feature TV change of 0.0025393. For HG02087, SEA prediction fell
  from 93.13% to 62.82%. The frozen decoder is therefore sensitive to small
  feature changes; precise output is not a safe drop-in performance upgrade.

The initial summary reader incorrectly expected an uncompressed global file.
It was corrected to read FLARE's `.global.anc.gz` outputs, without rerunning
inference. The exact original inference script was archived and verified
against its recorded hash; the corrected summary implementation has its own
hash. All input hashes were verified on summary recovery. Preparation now
automatically snapshots its implementation for future runs. All **48 tests**
passed, including compressed global parsing, export-only source transformation,
invalid inputs and archived-script provenance.

The precision-aware decoder retraining screen described below follows this
control, using consistently precise training, validation and pedigree features
alongside a matched rounded arm.
This control establishes quantization loss for the two tested pairs, not that
the recovered information can reliably discriminate EAS/SEA populations or
improve admixture prediction. Existing production files, labels and acceptance
gates remain unchanged.

## Matched precision-aware decoder retraining (2026-10-14)

[run_precision_retraining.sh](run_precision_retraining.sh) launches a
six-model fixed screen: Ridge, row-stochastic constrained mapping and
ExtraTrees, each trained and evaluated on rounded versus precise features.
It refuses existing output directories. To launch a new isolated run:

```bash
PRECISION_RETRAINING_OUTDIR=results_v3_precision_retraining_02 \
  bash run_precision_retraining.sh
```

Screen 01 completed successfully in
[results_v3_precision_retraining_01](results_v3_precision_retraining_01).
The matched comparison is
[matched_precision_screen.tsv](results_v3_precision_retraining_01/matched_precision_screen.tsv);
the summary was produced only after all 17 inference tasks and six model fits
completed.

The screen reuses all existing panel-ready genomes: 108 training families
with 492 descendants and 27 development families with 123 descendants.
It does not rerun Ped-sim or use old mosaic training examples. The audited
307 pure-training and 159 pure-validation donors, plus the known real-user
regression case, are re-inferred with the verified isolated precision jar.
The existing donor-disjointness and label audits must pass before launch.
Pure-validation and development founder donors remain disjoint from the
pure/pedigree decoder-training pool.

There are 17 inference tasks (three pure/real datasets, eleven training
pedigree batches and three development batches). Each task uses a fresh
network-disabled container, six CPUs, 40 GiB memory limit and the frozen
v1 model/settings. At most three inference containers run concurrently.
Existing FLARE models, installed jar and published results are read-only.
For every inferred posterior summary, the two feature arms are:

- **precision**: double-precision global export, normalized;
- **rounded**: the same precise values rounded with `np.rint` to three decimals,
  then normalized, matching the legacy extraction convention.

Both arms therefore share inference, targets, markers and ancestry posteriors;
only feature quantization differs. This is an explicit matched quantization
control, not an assumption that all regenerated features reproduce historical
outputs. Sample IDs/columns are verified exactly, and differences from the
historical normalized features are persisted per task in
`feature_precision_audit.tsv`.

Each of the six fits runs in its own fresh container. Hyperparameters, seed
20261013 and training weights match within every precision pair, reusing the
previous training-screen settings: Ridge alpha 0.001/no intercept, constrained
mapping alpha 1e-5, and ExtraTrees 200 trees/minimum leaf 10/full features.
All fits use 75% source-balanced pure training and 25% family-balanced
pedigree training, without calibration or hyperparameter search. The known
real-user constraints and original unified gates remain unchanged; insufficient
minor-ancestry support is explicitly marked unevaluable.

Every saved model must reproduce its development, pure-validation and
real-user predictions on reload. Input hashes and paired hyperparameters
are verified, and family-level precise-minus-rounded TV differences are
persisted. Models remain development-only and are not automatically promoted.
The cohorts were already inspected, so this is not independent final
validation or microarray/imputation certification. Implementation passed
**50 regression tests**, including quantization, normalization and exact
sample-aligned feature replacement.

### Completed screen 01

All regenerated rounded feature rows reproduced the historical features:
zero differing rows at tolerance 1e-12, maximum absolute difference 1.11e-16.
The three rounded model arms independently reproduced the prior matched
pure-plus-pedigree screen's development and pure-validation predictions to
tolerance 1e-10. Every saved model was independently reloaded in a further
fresh container and reproduced all 123 development, 159 pure-validation and
one real-user prediction to tolerance 1e-12. Coverage and input/output hashes
were verified; see
[persistence_audit.json](results_v3_precision_retraining_01/persistence_audit.json).

| Decoder, pure + pedigree training | Rounded pedigree mean TV | Precise pedigree mean TV | Rounded pure-validation mean TV | Precise pure-validation mean TV |
|---|---:|---:|---:|---:|
| Ridge | 0.163501 | 0.163491 | 0.233636 | 0.233660 |
| constrained linear mapping | 0.182048 | 0.182215 | 0.262431 | 0.262806 |
| ExtraTrees | 0.272250 | 0.270322 | 0.188932 | 0.183048 |

Precision made essentially no difference to the two linear models. ExtraTrees
improved by 0.001927 absolute pedigree TV and 0.005884 pure-validation TV, but
did not approach frozen v1's 0.16318 pedigree / 0.08785 pure-validation TV.
The precise ExtraTrees model's balanced p95 TV is 0.49826, above the unchanged
0.40 limit. All six models failed the known real-user hard gate and at least
one existing error gate. Minor-ancestry support remains 19, below the required
20; no model was promoted.

An independent companion audit persisted
[paired_collision_summary.tsv](results_v3_precision_retraining_01/paired_collision_summary.tsv)
and [pure_source_diagnostics.tsv](results_v3_precision_retraining_01/pure_source_diagnostics.tsv).
Across all 466 pure donors, higher precision eliminated both cross-source
exact-feature collision groups (37 affected rows became zero). Nevertheless,
SEA pure-validation mean TV barely changed for the linear models; ExtraTrees
changed from 0.40182 to 0.38583 across the 22 SEA validation donors.

This establishes that export quantization removes real information, but
recovering it is not a major overall performance improvement for these three
fixed decoder configurations. It does not establish a biological ceiling or
rule out benefits with other representations, learners or training coverage.
Further work should prioritize population-generalizing EAS/SEA discrimination
and representation/mixture consistency rather than another precision-only
hyperparameter sweep. Production files, ancestry labels and gates are unchanged.

## Four-model EAS/SEA representation screen (2026-10-14)

[run_eas_sea_model_screen.sh](run_eas_sea_model_screen.sh) implements the
targeted four-model screen: Ridge and ExtraTrees on each of the frozen
unrefined diagonal K30 and diagonal K36 representations. No new latent
representation is trained, no old mosaic examples are used and there is no
hyperparameter search. Precision remains at the legacy three-decimal level
to reuse existing pure, development and real-user features consistently.

The initial wrapper in
[results_v3_eas_sea_screen_01](results_v3_eas_sea_screen_01) missed the
`fractions` helper import after successful FLARE inference. The scheduler and
remaining containers were stopped before changing code. All eleven completed
K30 inference outputs were checksum-verified and recovered into
[results_v3_eas_sea_screen_02](results_v3_eas_sea_screen_02); incomplete K36
outputs were not reused. Exact old implementation/launcher snapshots and all
old input hashes were verified, and recovery provenance was persisted.
The missing import was fixed and an actual compressed-output export regression
test added. Screen 02 completed the remaining eleven K36 inference tasks,
all four decoder experiments and 124 population-fold fits. Existing results
were not overwritten.

To run a further isolated screen:

```bash
EAS_SEA_SCREEN_OUTDIR=results_v3_eas_sea_screen_03 \
  bash run_eas_sea_model_screen.sh
```

Only the 108-family training cohort requires new inference: 11 existing
panel-ready batches for each representation, 22 tasks total. All simulation
genomes, realized truth, founder assignments and the 27-family development
cohort are reused. Existing frozen artifacts and representation-audit hashes
must verify before inference. Each task uses the stock jar and frozen model,
seed 20261002, six CPUs and a 40 GiB container limit, with at most three
inference tasks in parallel. Preparation, each inference task, each main
model experiment and summary run in fresh network-disabled containers.

Both learners use the previously fixed training-screen settings:
Ridge alpha 0.001/no intercept, or ExtraTrees 200 trees/minimum leaf 10/full
features, seed 20261013. Main fits use all 307 pure-training donors and 492
training descendants, weighted 75% source-balanced pure / 25% family-balanced
pedigree, without calibration. Evaluation uses the same 123 development
descendants, 159 pure-validation donors and one real-user regression case.
Saved models must reproduce every prediction on reload. Original acceptance
gates and the minor-support threshold are unchanged; no model is automatically
promoted.

### Population-purged decoder controls

Each of the four experiments additionally refits its fixed learner once per
EAS/SEA validation fine population. It excludes all pure-training donors
from that population **and every training family containing any such founder**.
The remaining pure/pedigree weights are recomputed with the same 75/25 split.
Held-out validation donors never participate in decoder training. Purged donor
IDs, complete excluded families, remaining populations and model/prediction
files are persisted for every fold. Coverage, source completeness and
prediction reproduction are enforced.

These are **decoder-level population holdouts**: the frozen representation
and FLARE reference are not rebuilt without that population. They therefore
do not constitute fully unseen-population representation validation.
The checks include separate EAS and SEA TV, source-balanced dominant accuracy,
and the smaller pedigree-eligible EAS/SEA subset, which has only four SEA
validation donors.

The completed comparison is
[eas_sea_model_screen.tsv](results_v3_eas_sea_screen_02/eas_sea_model_screen.tsv),
alongside [frozen_baselines.tsv](results_v3_eas_sea_screen_02/frozen_baselines.tsv),
per-source metrics, main and population-held-out predictions, fold provenance
and development models. Summary was produced only after all inference and
model experiments finished. Preparation passed
**54 regression tests**, including compressed-output export, whole-family population purging,
source-balanced diagnostics and unchanged tail/real/minor-support gates.
The existing development cohort has already been inspected; any promising
model still needs independent, broader validation before deployment.

### Completed screen 02

| Model | Pedigree mean TV | Pure-validation mean TV | EAS pure TV | SEA pure TV |
|---|---:|---:|---:|---:|
| frozen v1 | 0.16318 | 0.08785 | 0.01537 | 0.33708 |
| K30 + Ridge, pure/pedigree trained | 0.15635 | 0.19222 | 0.24413 | 0.25963 |
| K30 + ExtraTrees, pure/pedigree trained | 0.30078 | 0.18875 | 0.24607 | 0.25554 |
| K36 + Ridge, pure/pedigree trained | 0.17638 | 0.20977 | 0.31451 | 0.31765 |
| K36 + ExtraTrees, pure/pedigree trained | 0.28920 | 0.18109 | 0.28516 | 0.26769 |

K30 + Ridge improves mean pedigree TV by 0.00682 absolute versus frozen v1,
and by 0.05564 versus K30's historical frozen decoder (0.21199). Its ordinary
pedigree p95 TV is 0.45654. This is a useful model-development lead, not an
overall replacement: pure-validation TV worsens markedly, and lower SEA TV
comes with much worse EAS proportions. Dominant-source accuracy alone does
not capture that loss of proportion quality.

Population-purged EAS/SEA balanced dominant accuracy is 0.81661 for K30/Ridge
and 0.76646 for K36/Ridge. Among the smaller pedigree-eligible subset it is
0.72917 and 0.43750 respectively. K30/Ridge's held-population pure EAS and SEA
mean TV are 0.30409 and 0.28859. These remain decoder-only holdouts, not a fully
unseen-population benchmark; only four eligible SEA validation donors are
available.

All four models fail the existing real-user hard gate and at least one error
gate. K30/Ridge's balanced p95 TV is 0.44467 (limit 0.40), and balanced
off-target mass is 0.14648 (limit 0.14). Minor support is still 19 and remains
explicitly unevaluable. No model was promoted and no threshold was relaxed.

An additional fresh-container persistence audit independently reloaded all
four main models and all 124 population-fold models. It reproduced all
development/reference/real-user predictions, all 51 EAS/SEA population-held-out
predictions per experiment, verified exact excluded donor/family sets, and
checked every recovered source-output hash:
[persistence_audit.json](results_v3_eas_sea_screen_02/persistence_audit.json).
The next development question is whether improved mixture behavior can be
retained without degrading pure-source proportions and other ancestry groups,
rather than treating a lower pedigree mean alone as success.

## Confirmed reference-overlap failure and fixes

Genome-wide controls kept the eight target genomes, frozen FLARE model,
decoder, calibrator and seed fixed. Targets were renamed to distinct study IDs,
as Ped-sim offspring are. The first inference used the original reference;
the second excluded only the five overlapping donors from its active
`ref-panel`. No model was retrained.

| Founder | Predicted OCE with original reference | Predicted OCE with donor-disjoint reference |
| --- | ---: | ---: |
| HGDP00490 (Melanesian) | 0.0% | 95.4% |
| HGDP00491 (Melanesian) | 0.0% | 100.0% |
| HGDP00540 (Papuan) | 0.0% | 100.0% |
| HGDP00545 (held-out Papuan control) | 100.0% | 100.0% |
| HGDP00662 (held-out Melanesian control) | 97.4% | 100.0% |

The three overlapping OCE founders originally had 100% `anc_14` features
(decoded mostly as SSA), changing to the OCE-associated `anc_26` after
reference removal. HGDP00832's predicted NAM changed from 2.3% to 90.2%.
HG00759 changed from mostly SSA to 58.6% EAS / 41.4% SEA; the KHV control
HG02142 remained approximately 88% SEA. All eight target genomes had the same
235,503 input marker records, zero missing genotypes and zero unphased
genotypes. FLARE retained 235,453 versus 235,450 markers.

This is strong experimental evidence that self/reference copying can break
latent-state inference for these founders under the frozen model. It does not
establish a universal FLARE defect or provide a corrected offspring leaderboard.
The 403 descendants from 114 families with none of the five overlapping donors
have mean TV 0.19113 for the historical winner and 0.16211 for frozen v1.
These are diagnostic subsets with a changed scenario/source distribution.

Implemented safeguards:

- Reference-disjoint selection is mandatory. Overlapping candidates are
  explicitly excluded and recorded, not silently relabeled as held-out.
- Train, validation, supplemental and reference-panel fine labels are checked
  against `ancestry_role_audit.tsv`. HG00759 is corrected to `1KG_CDX / EAS`,
  consistently with the training ontology; KHV remains SEA despite the broader
  1000 Genomes EAS superpopulation category. Supplemental OCE fine labels now
  use their actual Papuan/Melanesian labels, not the placeholder `OCE`.
- Actual training-mosaic donor IDs are also excluded, not just pure-training IDs.
- Before any expensive training, founder preflight attempts the complete
  requested plan and fails on insufficient distinct within-family donors.
- As requested, the **reduced** design rejects infeasible tree proposals,
  preserving all nine classes and within-family donor uniqueness. Its 252-family
  preflight has 28 families per class, using 120 eligible reference-disjoint
  donors: EUR 26, EAS 24, SAS 29, SEA 4, GME 8, NAM 4, OCE 2, SSA 23.
  This changes source combinations/dominance coverage; it is not equivalent to
  the full design. The manifest records rejected proposals and terminal-source
  combination counts. `PEDIGREE_DESIGN=full` keeps the original design and fails
  until enough independent donors are supplied.
- Prediction files, metrics and leaderboards identify the design and validated
  founder provenance. Unified selection rejects historical unvalidated
  predictions and mixed full/reduced designs. Launchers and founder preparation
  refuse to overwrite unvalidated historical pedigree results.

Validation passed: **23 unit tests**, the real 252-family preflight, offline
Nextflow previews for both normal and separate-output recovery paths, and a
real chromosome-1 regional Ped-sim smoke test covering all nine classes and
91 nodes. Founder assignments and resulting ancestry truth were verified.
The regional smoke test is an integration check, not an accuracy benchmark. Corrected
genome-wide offspring scores still require a new run. Do not relax the existing
selection gates or tune the decoder against the confounded historical scores.

## Historical pedigree performance audit (superseded interpretation)

The lowest pedigree mean TV is **0.26173354** for
`gmm_diag_unrefined_K30_ef013d39e0__extratrees__ridge__cdeb886004`.
TV is `0.5 * sum(abs(prediction - realized_truth))`, a proportion error, not
classification accuracy. Its pedigree p95 TV is 0.54522493; mean off-target
mass is 0.13744969. Its pure-source validation mean TV is 0.05611477.

The following checks were executed against the saved run:

- Independently reconstructed physical-span ancestry from all **572**
  chromosome/batch breakpoint files for all 2,548 planned individuals, checked
  equal haplotype spans and all **50,512** offspring parental-haplotype
  transmissions. Saved chromosome and genome-wide ancestry matched to floating
  point precision, with complete 22-autosome coverage.
- Compared **1,737,963 phased offspring genotypes** with the assigned original
  donor haplotypes at the recorded breakpoints: **zero mismatches**. This sampled
  batch 1 chromosomes 1, 10 and 22 and batch 26 chromosome 1; it is not an
  exhaustive genotype comparison of every chromosome/batch.
- Reloaded the frozen best decoder and calibrator and reproduced all 1,148
  saved predictions to floating point precision; independently recomputed TV.
- Verified no pedigree founder overlaps the donors named in the actual
  `training_mosaics8_metadata.tsv` donor-weight records. The five previously
  approved FLARE-reference overlaps were not excluded in the historical run.
  The controls above show that this was a major confounder, not merely a
  theoretical benchmark limitation.
- Checked alternative ancestry weights: prediction mean TV is **0.26204242**
  against genetic-map-weighted truth and **0.26148135** against panel-marker-
  weighted truth. Training mosaics use genetic-map-weighted targets whereas
  pedigree evaluation uses physical base pairs. These targets should be aligned
  in future training, but this difference does not explain the observed error.

A decoder control separates simulation from regression behavior. For 412 direct
founder crosses whose two donors have existing validation features, the
offspring FLARE feature vector differs from the ideal average of its founders'
feature vectors by mean TV **0.03639**. Passing those ideal mixed features
through the frozen decoder still gives mean ancestry TV **0.22115**, versus
**0.22026** for the actual offspring. Simply averaging the two already decoded
parent predictions gives **0.07189**. The corresponding 21 available `f1`
families have ideal-feature decoder TV **0.25512** versus actual TV **0.25922**.
This shows residual decoder generalization problems are possible, but the
control excluded reference-overlapping donors and does not explain their
catastrophic failures. It is not a production method or an independent final test.

The 1,008 training mosaics blend a median of **104 distinct donors** per genome
(minimum 15) and contain only two-, three- and four-source mixtures. Pedigree
F1s have two donors; deeper trees retain a small, family-specific founder pool,
and some scenarios contain eight sources. This distribution mismatch is a
strong candidate explanation for poor generalization despite good pure-source
validation. It has not been established as the only cause by retraining.

The aggregate includes intermediate non-founder ancestors, not just the final
descendant of each family. Restricting to the 252 terminal descendants gives
mean TV **0.24053**; the 28 `f1` families give **0.30675**. Thus poor performance
is not solely an artifact of counting intermediate individuals.

Both existing unit-test suites passed (14 model/selection tests and 4 pedigree
tests). No models, predictions, gates or simulation code were changed by this
initial audit. The original conclusion that these aggregate scores established
genuine model performance was too strong: simulation/truth integrity does not
exclude reference-copying failure. These scores are historical diagnostics,
not corrected performance, consumer-array accuracy or an unbiased final test.

## Changes and why

1. **100 trainable FLARE representations, not thousands of expensive Bayesian fits.** Configuration: `config/search_space_100.json`. One exploratory Bayesian DP GMM (K18, `max_iter=1200`, `tol=0.003`, `reg_covar=1e-5`, one initialization). Its `converged_` requirement remains strict; failing to converge is a real strategy failure. The 99 non-Bayesian candidates span diagonal/tied/spherical GMM, PCA/GMM, K-means, PCA/K-means and unrefined control. `gmm_tied_K24_293b7defcb`, the representation behind all five older finalists, is **included unchanged**. The frozen v1 control is also imported separately, so **101 including that fixed control**.
2. **Ped-sim participates in finalist selection.** A broad Stage 3 shortlist is scored with the same TV distance on reference individuals and realized pedigree-ancestry fractions, with macro-averaging over scenario, ancestry combination and family. Additional tail-TV, off-target mass, minor-ancestry recall and real-world regression safeguards remain. Stage 1 still uses reference-based screening; it now retains 24 representations with up to two from each family, and Stage 2 retains up to 80 candidates including at least one per surviving representation.
3. **REFINE_MODEL staging.** `path(config_json, stageAs:'refinement_input.config.json')` is used and the output config is copied from `refinement_input.config.json`, avoiding self-copy failures.
4. **Manifest-based Stage 2 selection.** Successful `FIT_STAGE2_CANDIDATE` output directory paths are concatenated by `collectFile()` into `stage2_result_dirs.txt`. `SELECT_STAGE2_CANDIDATES` stages the **one file**; neither selector stages thousands of candidate directories. Both selectors support `--manifest`, verify each line's `metrics.json` and `candidate.json`, and report a precise error for missing/incomplete entries. This implements the observed workaround; it does not establish the cause of exit 139.
5. **Recovery mode** (`--resume_stage2_manifest`): bypasses model builds, refinement, Stage 1, candidate generation and Stage 2 fits. Uses the existing published `selection/selected_representations.tsv`, `models/*.final.model`, `models/*.config.json`, `features/dev/*.features.tsv` and the supplied exact stage-2 manifest. Resumes selection, Stage 3, Ped-sim and unified finalist selection. It does not mark prior failed Bayesian fits successful.
6. **Preserves earlier outputs.** The corrected sweep defaults to `/workspace/results_v3_sweep_100_disjoint_reduced`. Recovery can read models/features from `RECOVERY_SOURCE_OUTDIR` while writing new pedigree results to `RUN_OUTDIR`. The work/cache directory is part of this project deployment; avoid modifying files in any currently active workspace.

## Composition: 100 trainable representations

| Family | Count |
|:--|--:|
| Diagonal GMM | 32 |
| Tied-covariance GMM | 14 |
| Spherical GMM | 14 |
| PCA + GMM | 12 |
| K-means | 16 |
| PCA + K-means | 8 |
| Unrefined diagonal GMM | 3 |
| Bayesian DP-GMM (exploratory) | **1** |
| **Total** | **100** |

With the frozen v1 baseline, there are 101 total representations **attempted/imported** if every trainable strategy succeeds. Unconverged candidates do not proceed.

## Deploy

Extract this package into a **new directory** on your SNPster host, e.g. `.../flare/ancestry_benchmarking_100/`. It includes a complete `main.nf`, configuration, scripts and tests.

Build the pedigree image from the existing benchmark image:

```bash
docker build -t snpster-ancestry-benchmark:pedigree-100 .
```

The Dockerfile adds a pinned Ped-sim build and its Boost/zlib dependencies while
preserving the existing Nextflow, Java, Python, `bcftools` and `tabix` tools.
The launcher defaults to this image and its `/opt/ped-sim/ped-sim` binary. Set
`PEDSIM_HOST_DIR` only to override it with a host installation.
The pedigree map is the existing four-column
`/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map`; the imputer's
`plink.GRCh38.map` path on this server is a directory, not a map file.
Run the launchers with `bash` if extracted files do not have executable permissions.

### Server founder supplementation

The historical five supplementary donors are retained in the metadata for
traceability, with corrected authoritative fine/source labels. All five are
now excluded because they overlap FLARE's reference. The resulting pool has
only 2 OCE, 4 NAM and 4 SEA donors; use the explicitly reduced design or supply
more independent genomes. Additional metadata is rejected for training overlap,
missing genotypes, unsupported sources or label disagreement. Increasing batch
size or reusing a donor within a family is not a solution to insufficient
within-family founders.

### Static checks

```bash
python3 -m unittest discover -s tests -v
python3 -m unittest discover -s pedigree_benchmark/tests -v
bash -n run_benchmark_v3.sh run_in_container_v3.sh
python3 scripts/generate_model_configs.py --search-space config/search_space_100.json --outdir /tmp/ancestry_100_sanity --max-models 0
```

For a dry run on your actual server after the configs have been generated (run inside Docker with the standard mounts):

```bash
nextflow run main.nf -preview -params-file config/params.json \
  --search_space /workspace/config/search_space_100.json \
  --generated_model_configs '/workspace/config/generated/model_configs/*.json' \
  --outdir /workspace/results_v3_sweep_100_disjoint_reduced \
  --pedigree_design reduced \
  --pedigree_enabled true
```

The live Nextflow preview passed on the SNPster server on 2026-10-09.

### Historical server run (finished; no passing finalists)

Container: `snpster-ancestry-pedigree-100`.
Image: `snpster-ancestry-benchmark:pedigree-100` (Ped-sim 1.4.7).
Settings: 100 trainable representations plus the fixed control, 6,000 decoder
candidates, 252 families, batch size 10, all 22 autosomes, 20 CPUs and 180 GiB RAM.
Output: `results_v3_sweep_100/` in this project. The older benchmark output is
untouched. Python output is unbuffered in the running container.

```bash
docker logs --follow --tail 60 snpster-ancestry-pedigree-100
docker stats --no-stream snpster-ancestry-pedigree-100
```

The regional Ped-sim check validates software integration, not genome-wide accuracy.

### Small end-to-end smoke test

```bash
MAX_MODELS=8 MAX_STAGE2_CANDIDATES=80 \
PEDIGREE_ENABLED=true PEDIGREE_FAMILIES=9 PEDIGREE_BATCH_SIZE=3 \
PEDIGREE_DESIGN=reduced RUN_OUTDIR=/workspace/results_v3_smoke_disjoint_reduced \
PEDIGREE_CHROMOSOMES=1,2 \
./run_in_container_v3.sh
```

Smoke-test only: two chromosomes do not provide full-genome accuracy evidence. This produces tasks under the new output directory by default. Use another output directory or clear the generated configuration safely before running the full sweep: do not mix 8-model and 100-model outputs with different scientific parameters.

### Full new search

```bash
MAX_MODELS=0 MAX_STAGE2_CANDIDATES=6000 \
PEDIGREE_ENABLED=true PEDIGREE_FAMILIES=252 PEDIGREE_BATCH_SIZE=10 \
PEDIGREE_DESIGN=reduced \
PEDIGREE_CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
RUN_OUTDIR=/workspace/results_v3_sweep_100_disjoint_reduced \
./run_in_container_v3.sh
```

The original full design cycles through 28 pairwise ancestry source combinations
and nine classes. The reduced design keeps equal class counts but changes
combination/dominance coverage when proposals cannot be built from the available
founders. Founder reuse **across different families** is on by default; within
each pedigree founders remain distinct. Use a separate output root for smoke
versus full for unambiguous provenance. Ped-sim uses a converted sex-averaged
PLINK map with the **Poisson crossover model** by default; do not claim
sex-specific interference unless compatible maps and parameters are supplied.

## Resuming after a panel merge failure

`PDB_MERGE_PANEL` concatenates chromosome VCFs in numeric chromosome order,
without `--allow-overlaps`. These non-overlapping inputs do not need indexes;
the merged output is still indexed with `tabix`. The former overlap mode
required unstaged input indexes and failed with `could not load index`.

The fix only changes the merge task, leaving completed simulations and model
fits eligible for the existing Nextflow cache. The project is bind-mounted into
the named container, so no image rebuild is needed. With the failed container
stopped, resume its original command and environment:

```bash
docker start snpster-ancestry-pedigree-100
docker logs --follow --tail 60 snpster-ancestry-pedigree-100
```

Do not delete `work/` or `.nextflow/`, or change the scientific settings.
The launcher already uses `-resume`.

The launcher now generates configs in a temporary directory and replaces only
files whose contents changed. It preserves unchanged input timestamps, which
Nextflow includes in its default cache keys. Previously, deleting and recreating
identical configs at startup caused model builds to restart despite `-resume`.
If a prior restart already recreated those files, the launcher fix alone cannot
restore the original cache keys; their original metadata must also be recovered
or the documented stage-2 recovery mode used.

Targeted regression test (requires Nextflow, bcftools and tabix):

```bash
bash tests/test_merge_panel.sh
bash tests/test_cached_configs.sh
```

This tests unindexed inputs supplied in scrambled chromosome order, sample
order, phased genotypes, record counts, and indexed access to the merged VCF.

## If stage 2 fails again

First preserve the exact input set used for the failed selection. A new full run can generate a manifest automatically; otherwise reconstruct it from the failed task's input paths as you did in the verified recovery, **not** from a mixed glob of multiple experiments.

```bash
# INSIDE the correct container and old or new benchmark directory:
# Read existing model/feature artifacts without overwriting old pedigree scores.
RUN_OUTDIR=/workspace/results_v3_sweep_100_disjoint_reduced \
RECOVERY_SOURCE_OUTDIR=/workspace/results_v3_sweep_100 \
PEDIGREE_DESIGN=reduced \
SEARCH_SPACE=config/search_space_100.json \
RECOVERY_STAGE2_MANIFEST=/workspace/results_v3_sweep_100/selection/recovery_stage2_result_dirs.txt \
PEDIGREE_ENABLED=true PEDIGREE_FAMILIES=252 PEDIGREE_BATCH_SIZE=10 \
./run_benchmark_v3.sh
```

The manifest path above is the existing 100-model run's durable manifest on
this server, repointed to verified published outputs before work cleanup.
Preserve it and its referenced completed candidate directories;
on another deployment, supply that run's exact manifest rather than a broad
glob. This command is for a new corrected run; it was previewed, not launched.

For the **October 9 old 27,121-fit recovery**, use
`SEARCH_SPACE=config/search_space.json`, `RECOVERY_SOURCE_OUTDIR=/workspace/results_v3`,
its **exact** stage-2 manifest, a fresh `RUN_OUTDIR`, and the same dataset/decoder
configuration. Do **not** combine old fits with the 100-model search space or a
changed training objective. This mode uses no generated model configs.

## Interpretation / limitations

- The old top five all share the same tied-GMM K24 representation. Despite passing Lonni, none passed the old absolute metric gates; the top had mean TV about `0.1248` and p95 TV about `0.5856`. They are controls, not established winners.
- Full unified model selection continues to use **clean high-density Ped-sim descendants projected to FLARE sites**, NOT the consumer-array → harmonization → Beagle → QC production data path. The Lonni production-domain problem is not yet proven solved.
- `prepare_founders.py` enforces training-mosaic and FLARE-reference disjointness and validates sample/source/fine labels. A separate, independently sourced final test remains necessary to eliminate selection bias and dependence between families.
- The unified scoring code requires every configured scenario and can refuse to choose any model if no candidate meets the declared gates. This is intentional. Tune numerical thresholds only as a documented, pre-specified decision, not to force a passing finalist.
- All chromosomes should be present for the final test, and the VCF/GRCh38 map sample/contig IDs must be verified on the actual server.
- Using many pedigree descendants from a limited donor pool does not create thousands of independent founder ancestries. Balance families and evaluate a separately held-out pool.
