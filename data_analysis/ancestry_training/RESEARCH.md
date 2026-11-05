# Evidence behind the standardized comparison

## Repository lessons carried forward

The completed experiments and rejection criteria are documented in
[the historical sweep record](../FLARE/SNPster_v3_100model_pedigree/README_NEW_SWEEP.md)
and [the v2 benchmark contract](../FLARE/SNPster_v3_100model_pedigree/config/ancestry_benchmark_contract_v2.json).
The new workflow does not rewrite those contracts or their running producers.

- The real chr22 imputation concordance failure was about **0.7711**, not an
  identity/key bug. Do not weaken the production 0.90 gate, substitute reference
  calls for truth, or confuse successful model execution with validity.
- More pure-source accuracy did not guarantee mixed/pedigree accuracy. Rich
  PCA, fine-profile and kernel features did not reliably solve domain/phase
  transfer. Include pure/mixed validation breakdowns and a simple frequency
  control; do not assume a large tree ensemble will beat native FLARE.
- Reference/donor overlap invalidated earlier favorable results. Disjoint sample
  IDs alone are insufficient: families and founder donors must also be purged,
  including against the external imputation reference.
- Calibrating on previously inspected evaluation labels is not independent
  evidence. Reserve reference, base-training, calibration, validation and test
  roles before fitting; keep native versus learned results distinguishable.
- Realized physical ancestry is the estimand. Use native local posteriors and
  chromosome spans, not marker-count averages, nominal pedigree fractions or
  rounded global-export features.

The historical contract specified an imputation-only benchmark and a small
candidate budget. The new user-requested workflow explicitly supports **two named
input-stage tracks**, but does not silently reclassify old evidence or launch a
new large grid. Its full-mode completeness/metric discipline follows the same
physical and equal-person principles.

## Primary-source research

### FLARE first

[FLARE upstream](https://github.com/browning-lab/flare) and its
[Apache-2.0 license](https://github.com/browning-lab/flare/blob/master/LICENSE-2.0.txt)
make it a sensible primary software component when future commercial use
matters. This is **not** legal clearance of reference data, maps, dependencies,
patents or a packaged service.

The [original FLARE paper](https://doi.org/10.1016/j.ajhg.2022.12.010) describes
scalable local ancestry inference. Upstream
[parameter estimation](https://github.com/browning-lab/flare/blob/master/src/admix/ParamEstimator.java)
and [options](https://github.com/browning-lab/flare/blob/master/src/admix/AdmixPar.java)
show why a supplied model does not automatically freeze target fitting:
`em`/generation/copying options still matter. Keep target inference singleton
when comparing co-batch-independent configurations, save emitted native models,
and tune declared FLARE hyperparameters rather than pretending `k` is a
supervised ethnicity-count option.

Upstream
[global ancestry export](https://github.com/browning-lab/flare/blob/master/src/admix/EstimatedGlobalAncProportions.java)
uses marker/haplotype averages. That is not this project's full-span physical
estimand, so aggregate ANP1/ANP2 local posteriors instead.

### Clusters are representations, not automatic ancestries

FLARE2's
[create_model_file.py](https://github.com/browning-lab/flare/blob/master/create_model_file.py)
uses a Gaussian mixture of reference copying profiles. Its `anc_0`, `anc_1`, ...
axes are cluster labels, not a guaranteed biological source order.
See also the [FLARE2 preprint](https://doi.org/10.1101/2025.10.13.681993).
A different `k` changes that representation; it requires an explicit stable
source map or a learned map fitted on donor-purged non-evaluation people.
The implemented PCA/GMM control exercises cluster interpretation, but **is not
an implementation or validation of FLARE2's copying-profile clustering**.
Native FLARE2 fitting is a separate future adapter, not an assumed capability of
the pinned local FLARE JAR.

### Gnomix as research comparator

[Gnomix upstream](https://github.com/AI-sandbox/gnomix) and its
[license](https://github.com/AI-sandbox/gnomix/blob/main/LICENSE.md) restrict the
software to academic research use without a separate commercial arrangement.
Keep it opt-in and in a distinct research-only selection track.

The workflow uses the genuine pinned upstream fast LR+CRF architecture, not a
tree model renamed Gnomix. A shared reference-only LR base is fitted once per
window/context/reference policy. CRF L2 is applied to
`smooth.model.CRF.set_params(c2=...)` before fitting on separate labelled
training haplotypes. This avoids both duplicate bases and a hyperparameter
record that never changes the actual estimator.
The implementation deliberately does not invoke upstream admixture simulation
on an evaluation cohort.

### FLARE plus a learned global decoder

A small nonnegative row-stochastic map is the most transparent first learned
interpretation of latent component proportions. It maps a simplex to the named
source simplex and exposes the fitted mapping. Ridge is a low-capacity control;
RandomForest/XGBoost can learn nonlinear interactions but can also overfit small
calibration cohorts and mask-specific errors. Fit on truly held-out calibration
features, select on validation, and test only the locked finalists. Preserve
native results even when a learned decoder wins.

The implementation also includes PCA/Ridge, PCA/RF and PCA/XGBoost global
controls; they are genuine comparisons, not a claim that PCA or trees should
replace FLARE. A promising future extension is a donor-purged classifier/frequency
feature stack like the strongest historical phase-free development model.
It is **not** silently included here: it needs its own frozen feature contract,
out-of-fold donor-purged training budget and transfer validation.

## Compatibility and limits of this evidence

Research above describes upstream behavior, not an assurance that every new
upstream option is supported by the locally pinned JAR. The implementation
allowlists only supported surfaces, records exact JAR/source hashes, and tests
real task containers. `impute=false` belongs to Beagle phase-only preprocessing;
it must not be passed to the pinned FLARE JAR.

No new large biological experiment is part of this implementation. Synthetic
smoke scores cannot establish the best hyperparameters, the best interpretation,
clinical validity or a commercially deployable model. Use the standardized
validation and frozen independent test to answer those questions later.
