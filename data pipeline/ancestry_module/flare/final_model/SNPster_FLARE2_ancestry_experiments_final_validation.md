# SNPster FLARE2 ancestry-model experiments

**Date:** 2026-10-05  
**Project:** SNPster ancestry inference  
**Purpose:** Develop a robust local/global ancestry model that avoids the pathological Central Asian (CAS) sink seen with the original broad-panel FLARE setup.

> **Frozen production-model status (2026-10-05):** The model-development phase described below is complete. Sections 1-26 are retained as the historical experimental record. The final production decision, metric definitions, validation results, and integration guidance are documented in Sections 27 onward. No further hyperparameter tuning should be performed against the development or final synthetic benchmark sets unless a new reproducible failure mode is identified.

## 1. Problem being addressed

The original SNPster ancestry workflow used FLARE with broad reference panels:

- EUR
- EAS
- SAS
- SEA
- OCE
- SSA
- MID
- CAS
- NAM
- AMR

For some admixed users, `CAS` behaved as a sink because Central Asian references such as Kazakh, Uygur, Tajik, and Turkmen are genetically intermediate/admixture-rich.

### Lonni: known approximately half European / half Filipino

Original broad FLARE chromosome 1 result:

| Broad ancestry | Proportion |
|---|---:|
| CAS | 58.3% |
| SEA | 21.0% |
| EUR | 16.4% |
| OCE | 2.3% |
| EAS | 1.6% |
| SAS | 0.2% |
| MID | 0.1% |

This was biologically implausible.

Removing CAS gave a much more plausible chromosome 1 result:

| Broad ancestry | Proportion |
|---|---:|
| EUR | 45.5% |
| SEA | 31.2% |
| EAS | 16.9% |
| OCE | 2.9% |
| SAS | 2.6% |
| MID | 0.8% |
| NAM | 0.1% |
| AMR | 0.1% |

However, removing CAS is not acceptable for genuine Central Asian users.

## 2. Reference data

The merged ancestry reference contains **4,886 samples**.

### HGDP + 1000 Genomes

- 4,091 samples
- Metadata from `gnomad_meta_updated.tsv`
- Fine populations include CEU, GBR, FIN, IBS, TSI, CHB, CHS, JPT, KHV, PJL, GIH, YRI, ESN, Yoruba, French, Russian, Uygur, Hazara, Papuan, Melanesian, and many others.

### GSE77508

196 Southeast Asian samples including:

- Bajo
- Burmese
- Dusun
- Kankanaey / Igorot
- Lebbo
- Luzon Filipino
- Malay
- Murut
- Vietnamese
- Vizayan Filipino

### GSE73996

381 selected non-Siberian samples including:

- Avar
- British
- Dinka
- Egyptian
- German
- Han South
- Indian
- Indonesian Java
- Iranian
- Italian
- Kazakh
- Korean
- Manchu
- Mari
- Mongol
- Mordva
- Papua New Guinea
- Sri Lankan
- Tajik
- Tibetan
- Turk
- Turkmen
- Uygur
- Vietnamese

### Kazakh PRJEB89820

- 218 samples
- Fine population label: `Kazakh2025_Kazakh`

## 3. Fine reference panel

A fine-grained FLARE panel was created with **113 fine populations**.

Relevant files:

```text
/ref/flare_reference.chr1.vcf.gz
...
/ref/flare_reference.chr22.vcf.gz

/ref/flare_reference.fine.ref.panel
/ref/fine_sample_metadata.tsv
```

The fine panel contains all 4,886 reference samples.

The goal is to keep populations such as Kazakh, Uygur, Filipino, German, CEU, Yoruba, etc. separate at the donor-panel level instead of immediately forcing them into 10 broad ancestry states.

## 4. Direct 113-population FLARE experiment

FLARE was first run directly with the 113 fine populations.

For Lonni, chromosome 1 had large contributions from Luzon Filipino, Uygur, Kazakh, GBR, Malay, German, Mordva, British, French, Russian, and others.

After manually collapsing those fine populations back to broad groups:

| Broad group | Proportion |
|---|---:|
| CAS | 38.3% |
| SEA | 34.0% |
| EUR | 26.7% |
| SAS | 0.3% |
| OCE | 0.6% |
| MID | 0.1% |

**Conclusion:** simply increasing the number of direct ancestry labels did not solve the CAS problem.

## 5. FLARE2 approach

FLARE 0.6.0 supports the FLARE2-style workflow:

1. Run FLARE with `panel-probs=true`
2. Produce a `.panels` file containing reference-panel copying profiles
3. Cluster those copying profiles using `create_model_file.py`
4. Use the learned model with FLARE
5. Refine copying probabilities with `update-p=true` and `em=true`
6. Freeze the learned model for new users with `update-p=false`

The conceptual change is:

> Fine populations remain donor panels, while latent ancestry states are learned from patterns of copying across those panels.

A latent ancestry state therefore does not have to correspond one-to-one with a population such as Kazakh, Filipino, CEU, or Uygur.

## 6. One-person FLARE2 experiment

Lonni alone was initially used to generate FLARE2 panel-copying profiles.

K values of 10, 15, 20, 25, and 30 were explored.

This reduced the CAS problem somewhat, but the inferred states were diffuse.

**Conclusion:** one target individual does not provide enough haplotype diversity for robust FLARE2 clustering.

## 7. Seventeen-person calibration experiment

A set of 17 SNPster users was then used as calibration targets.

Broad primary ancestries included:

- AMR
- CAS
- EAS
- EUR
- MID
- SAS
- SEA
- SSA

Examples:

```text
173      AMR   0.472662
murtaza  CAS   0.508188
278      EAS   0.844461
145      EUR   0.934343
285      MID   0.486351
162      SAS   0.907363
17       SEA   0.689598
57       SSA   0.679964
```

The latent-state structure improved markedly.

However, these labels came from a previous ancestry method and some users had only ~50-70% primary ancestry.

**Conclusion:** more target haplotypes helped, but the calibration labels were partly circular.

## 8. Held-out-reference calibration

The next experiment used known reference individuals themselves as calibration targets.

From the original 4,886 reference samples:

- **552 samples were held out as FLARE targets**
- **4,334 samples remained as the donor reference**

The held-out samples were stratified by fine population approximately as:

```text
>=100 samples -> hold out 10
50-99         -> hold out 8
20-49         -> hold out 5
10-19         -> hold out 3
6-9           -> hold out 1
<6            -> hold out 0
```

Held-out counts by broad group:

| Broad group | Held-out samples |
|---|---:|
| AMR | 39 |
| CAS | 30 |
| EAS | 92 |
| EUR | 88 |
| MID | 36 |
| NAM | 12 |
| OCE | 9 |
| SAS | 95 |
| SEA | 61 |
| SSA | 90 |
| **Total** | **552** |

Relevant files:

```text
/output/ref_holdout/holdout_samples.txt
/output/ref_holdout/reference_samples.txt
/output/ref_holdout/holdout_metadata.tsv
/output/ref_holdout/reference_reduced.ref.panel

/output/ref_holdout/holdout_targets.chr1.vcf.gz
/output/ref_holdout/reference_reduced.chr1.vcf.gz
```

The chromosome 1 FLARE panel-probability run used:

```text
reference samples : 4334
target samples    : 552
markers           : 19344
```

Panel-copying output:

```text
/output/ref_holdout/calibration.panel_probs.chr1.panels
```

## 9. Fast K10 proof-of-concept

The official clustering step became computationally expensive with 552 targets, so an exploratory fast version of `create_model_file.py` was made.

The fast Gaussian-mixture configuration used:

```python
GaussianMixture(
    n_components=nanc,
    covariance_type="diag",
    init_params="k-means++",
    n_init=2,
    max_iter=100,
    random_state=20261001
)
```

A K10 model completed successfully:

```text
/output/ref_holdout/fast.K10.model
```

It was refined using the 552 held-out targets:

```text
/output/ref_holdout/fast.K10.final.chr1.model
```

with:

```text
update-p=true
em=true
```

For independent users the refined model was frozen with:

```text
update-p=false
```

## 10. K10 latent-state interpretation from known holdouts

Mean latent-state use among the known held-out broad groups:

### anc_0

```text
SEA  0.1704
OCE  0.0063
```

Strongly SEA-associated.

### anc_1

```text
EAS  0.9090
SEA  0.7300
CAS  0.3620
```

East / East-Southeast Eurasian axis.

### anc_2

```text
NAM  0.8979
AMR  0.1875
```

Strongly Native-American-associated.

### anc_3

```text
OCE  0.9840
SEA  0.0334
```

Strong Oceanian state.

### anc_4

```text
SAS  0.0457
```

Rare SAS-associated component.

### anc_5

```text
SSA  0.0339
```

Rare SSA-associated component.

### anc_6

```text
SSA  0.6762
MID  0.0441
AMR  0.0420
```

Strong Sub-Saharan African state.

### anc_7

```text
AMR  0.5689
SSA  0.2725
EUR  0.1480
SAS  0.1035
```

Broad admixed-American / mixed state.

### anc_8

```text
EUR  0.8415
MID  0.7882
CAS  0.3411
```

Strong western-Eurasian state. K10 does not cleanly separate Europe from the Middle East.

### anc_9

```text
SAS  0.8235
CAS  0.2869
MID  0.1044
```

Strong South-Asian-oriented state.

## 11. Held-out broad-group hard classification

Assigning each latent state to the broad group with the highest mean produced:

```text
Accuracy: 360 / 552 = 65.2%
```

This is not the main success criterion.

The confusion matrix showed systematic overlaps:

- MID largely used the same western-Eurasian state as EUR
- CAS decomposed among eastern-Eurasian, western-Eurasian, and SAS-oriented states
- SEA shared much of the main eastern-Eurasian state with EAS
- AMR used mixed ancestry states

Most importantly:

> No K10 component behaves as a dedicated CAS sink.

Central Asian held-outs instead use combinations of `anc_1`, `anc_8`, and `anc_9`.

## 12. Independent validation

### Frederik — known European

Chromosome 1:

```text
anc_8 = 1.000
all other states = 0
```

Genome-wide chromosomes 1-22:

```text
anc_8 = 100.00%
all other states = 0%
```

Every chromosome was individually 100% `anc_8`.

This strongly supports `anc_8` as the European/western-Eurasian state for a European individual.

### Lonni — approximately half European / half Filipino

Chromosome 1:

```text
anc_0   0.012
anc_1   0.444
anc_3   0.032
anc_8   0.466
anc_9   0.046
```

The former CAS sink disappeared.

Genome-wide, genetic-length weighted:

| State | Proportion |
|---|---:|
| anc_8 | 44.01% |
| anc_1 | 41.95% |
| anc_9 | 12.00% |
| anc_3 | 1.11% |
| anc_0 | 0.54% |
| anc_6 | 0.29% |
| anc_7 | 0.07% |
| anc_4 | 0.02% |
| anc_2 | 0.01% |
| anc_5 | 0.00% |

This is broadly compatible with an approximately half-European / half-Filipino genome.

The 12% `anc_9` should **not** currently be interpreted literally as 12% South Asian ancestry. It is a latent copying state and may partly absorb Asian haplotypes that K10 cannot resolve more finely.

Examples of chromosome variation:

```text
chr1:
anc_1 = 44.4%
anc_8 = 46.6%
anc_9 = 4.6%

chr8:
anc_1 = 57.2%
anc_8 = 35.1%

chr21:
anc_9 = 79.5%
anc_8 = 16.2%
anc_1 = 4.3%
```

### Murtaza — Central Asian / Afghan-related

Chromosome 1:

```text
anc_1 = 40.4%
anc_8 = 11.0%
anc_9 = 48.4%
anc_7 = 0.3%
```

Instead of collapsing into a ~55% CAS state, the ancestry is decomposed across eastern-Eurasian, western-Eurasian, and SAS-oriented states.

## 13. Main conclusions

### Original broad CAS was the major failure mode

Central Asian references are too genetically intermediate to function cleanly as one direct ancestry state.

### Direct use of 113 fine populations was not enough

More labels alone did not solve the problem.

### FLARE2 latent-state learning is substantially more promising

The model can learn states from copying patterns across many fine donor populations.

### More calibration haplotypes improved the model

Observed progression:

```text
1 target
    -> diffuse / unstable clustering

17 SNPster targets
    -> more structured, but labels partly circular

552 known held-out reference targets
    -> clean and interpretable ancestry axes
```

### Fast K10 is a successful proof of concept

It produced:

- Frederik: 100% western-Eurasian state genome-wide
- Lonni: approximately balanced western vs eastern/Asian ancestry genome-wide
- Murtaza: decomposed Central Asian ancestry rather than a CAS sink

### K10 is likely too coarse for final use

It does not cleanly separate:

- EUR vs MID
- EAS vs SEA
- some Central Asian substructure

Higher K may give better resolution.

### Latent states are not yet consumer-facing ancestry labels

States such as `anc_1`, `anc_8`, and `anc_9` represent recurring haplotype-copying patterns, not necessarily discrete ethnicity percentages.

## 14. Current important paths

Full reference:

```text
/ref/flare_reference.chr1.vcf.gz
...
/ref/flare_reference.chr22.vcf.gz

/ref/flare_reference.fine.ref.panel
/ref/fine_sample_metadata.tsv
```

Genetic maps:

```text
/maps/plink.chr1.GRCh38.map
...
/maps/plink.chr22.GRCh38.map
```

Held-out calibration:

```text
/output/ref_holdout/holdout_samples.txt
/output/ref_holdout/reference_samples.txt
/output/ref_holdout/holdout_metadata.tsv

/output/ref_holdout/holdout_targets.chr1.vcf.gz
/output/ref_holdout/reference_reduced.chr1.vcf.gz
/output/ref_holdout/reference_reduced.ref.panel
```

Panel-copying data:

```text
/output/ref_holdout/calibration.panel_probs.chr1.panels
```

Successful K10 model:

```text
/output/ref_holdout/fast.K10.model
/output/ref_holdout/fast.K10.final.chr1.model
```

Reduced genome-wide references:

```text
/output/ref_holdout/reduced_refs/reference_reduced.chr1.vcf.gz
...
/output/ref_holdout/reduced_refs/reference_reduced.chr22.vcf.gz
```

## 15. Next experiment: heavy K15 / K20 / K25

The next step is to use the original/heavier `create_model_file.py` rather than the exploratory diagonal-covariance version.

Models:

```text
K = 15
K = 20
K = 25
```

Input:

```text
/output/ref_holdout/calibration.panel_probs.chr1.panels
```

Suggested outputs:

```text
/output/ref_holdout/heavy.K15.model
/output/ref_holdout/heavy.K20.model
/output/ref_holdout/heavy.K25.model
```

Each completed model should then be refined using the same 552 held-out targets and 4,334-reference split.

The comparison should focus on:

1. whether EUR and MID separate better
2. whether EAS and SEA separate better
3. whether Central Asian haplotypes remain decomposed rather than forming a pathological sink
4. whether Frederik remains cleanly European
5. whether Lonni remains approximately half western-Eurasian / half eastern-Asian
6. whether Murtaza remains a sensible multi-axis Eurasian mixture
7. genome-wide consistency

The desired outcome is not necessarily a one-to-one match between K states and broad ancestry labels.

The desired outcome is:

> Stable latent ancestry states that reflect real haplotype structure, avoid pathological sink populations, and generalize sensibly to independent users.

## 16. Important caveats

- Current model construction has used chromosome 1 copying profiles.
- K10 used an exploratory diagonal-covariance GMM and is not yet the intended final production model.
- Some reference populations are themselves admixed.
- Broad geographic ancestry groups are human-defined categories and need not correspond to discrete genetic states.
- Local ancestry should remain interpreted at the latent-state level until state mappings are validated more broadly.
- Final production validation should use chromosomes 1-22 and substantially more independent people.

## Current status

The held-out-reference FLARE2 strategy has solved the most obvious failure observed in the original ancestry workflow:

> A half-European / half-Filipino individual is no longer classified primarily as Central Asian, while a European individual is cleanly assigned to the western-Eurasian state and a Central Asian individual is represented as a combination of Eurasian ancestry axes.

The next step is to determine whether full/heavy K15, K20, or K25 clustering gives better ancestry resolution than the successful K10 proof of concept.


---

## Appendix A — exact heavy-model workflow

Run the original `create_model_file.py`, not `create_model_file_fast.py`.

### Build heavy K15, K20, and K25 models

```bash
for K in 15 20 25; do
    echo
    echo "=================================================="
    echo "Building HEAVY K=${K}"
    echo "=================================================="

    rm -f \
      /output/ref_holdout/heavy.K${K}.model \
      /output/ref_holdout/heavy.K${K}.log

    /usr/bin/time -v \
      python3 /output/create_model_file.py \
        "$K" \
        /output/ref_holdout/calibration.panel_probs.chr1.panels \
        /output/ref_holdout/heavy.K${K}

    STATUS=$?

    if [ "$STATUS" -ne 0 ]; then
        echo "K=${K} failed with exit code ${STATUS}"
        break
    fi

    if [ ! -s /output/ref_holdout/heavy.K${K}.model ]; then
        echo "K=${K} produced an empty model file"
        break
    fi

    echo "K=${K} completed successfully"
    ls -lh \
      /output/ref_holdout/heavy.K${K}.model \
      /output/ref_holdout/heavy.K${K}.log
done
```

### Refine every completed heavy model

```bash
for K in 15 20 25; do
    MODEL=/output/ref_holdout/heavy.K${K}.model

    if [ ! -s "$MODEL" ]; then
        echo "Skipping K=${K}: model not present"
        continue
    fi

    echo "===== refining heavy K=${K} ====="

    flare \
      gt=/output/ref_holdout/holdout_targets.chr1.vcf.gz \
      ref=/output/ref_holdout/reference_reduced.chr1.vcf.gz \
      ref-panel=/output/ref_holdout/reference_reduced.ref.panel \
      map=/maps/plink.chr1.GRCh38.map \
      model="$MODEL" \
      update-p=true \
      em=true \
      array=false \
      probs=true \
      nthreads=20 \
      seed=20261001 \
      out=/output/ref_holdout/heavy.K${K}.final.chr1
done
```

### Check completed outputs

```bash
ls -lh /output/ref_holdout/heavy.K*
```

A valid model must be non-empty. A zero-byte `.model` file means model construction did not finish successfully.

### Independent validation

For a new user, use the refined model with:

```text
update-p=false
```

The held-out calibration targets must not be used to update the model during independent validation.

## 19. All-autosome model-selection experiment

Model selection was repeated using copying-profile information from **all 22 autosomes**.

The 552 reference individuals previously held out from the 4,334-person FLARE donor reference were split by fine population into a model-training subset and a completely untouched validation subset. The validation subset did not contribute to `panel-probs`, GMM clustering, or `update-p=true` refinement.

### Validation performance

| Model | Accuracy | Balanced accuracy |
|---|---:|---:|
| diag K10 | 67.38% | 68.70% |
| diag K12 | 63.10% | 72.36% |
| diag K15 | 65.78% | 72.25% |
| diag K15, reg=1e-5 | 71.66% | **78.55%** |
| **diag K18** | **75.94%** | **75.98%** |
| diag K20 | 70.05% | 77.61% |
| diag K20, reg=1e-5 | 66.84% | 74.16% |
| tied K15 | 62.57% | 69.88% |
| tied K20 | 63.64% | 70.85% |

Diagonal covariance generalized better than tied covariance and the earlier full-covariance models. `diag_K18` became the leading model because it had the best ordinary out-of-sample accuracy while maintaining nearly identical balanced accuracy.

## 20. Independent known-person validation

### diag K18

```text
Frederik:
anc_16 100.0%

Lonni:
anc_16  42.6%
anc_2   31.0%
anc_6   14.7%
anc_3    5.7%
anc_8    4.8%
anc_9    1.0%

Murtaza:
anc_6   34.6%
anc_8   28.5%
anc_3   21.9%
anc_16  11.7%
anc_2    3.3%
```

This was the most coherent independent behavior among the shortlisted all-autosome models.

## 21. Fine-population interpretation of diag K18

### `anc_2` — Southeast Asian / Filipino-associated

Examples:

```text
Kankanaey / Igorot   0.980
Indonesia-Java       0.959
Luzon Filipino       0.944
Dusun                 0.869
Malay                 0.860
Cambodian             0.832
Bajo                  0.719
Murut                 0.702
```

### `anc_3` — East Asian-associated

```text
Japanese / JPT       ~0.998
Korean                0.994
Manchu                0.972
Han                   0.964
CHB                   0.944
Tibetan               0.929
```

### `anc_6` — Northern / Central Eurasian-associated

```text
Yakut                 0.903
Kazakh2025            0.710
GSE Kazakh            0.619
Mongol                0.603
Hezhen                0.577
Oroqen                0.538
HGDP Uygur            0.448
Hazara                0.293
GSE Uygur             0.286
```

### `anc_8` — South Asian-associated

```text
ITU                   0.994
Indian                0.989
Sri Lankan            0.989
Brahui                0.955
Sindhi                0.944
Burusho               0.937
Pathan                0.920
GIH                   0.911
```

### `anc_16` — Western Eurasian / European-associated

Most European reference groups are overwhelmingly `anc_16`; Middle-Eastern populations also use it substantially.

```text
French                0.997
Basque                0.997
Italian               0.996
German                0.996
Sardinian             0.992
FIN                   0.992
TSI                   0.992
GBR                   0.971

Turk                  0.821
Bedouin               0.812
Palestinian           0.773
Egyptian              0.768
Iranian               0.679
```

Therefore `anc_16` is best treated as a western-Eurasian state rather than as a universally literal European ethnicity percentage.

## 22. Remaining question from Lonni

The leading K18 model gives Lonni:

```text
anc_16  42.6%  western Eurasian
anc_2   31.0%  Southeast Asian / Filipino
anc_6   14.7%  northern / central Eurasian
anc_3    5.7%  East Asian
anc_8    4.8%  South Asian
anc_9    1.0%  Oceanian
```

The overall western-vs-Asian split is plausible for an approximately half-European / half-Filipino individual.

However, the `anc_6` contribution needs controlled validation because the strongest Filipino reference populations are predominantly `anc_2` and generally have little `anc_6`.

## 23. Synthetic 50/50 admixture benchmark

### Objective

Create artificial diploid individuals whose ancestry is known exactly and test whether K18 reconstructs the expected ancestry mixture without producing spurious sink components.

### Design

Synthetic source individuals are selected **only from the untouched validation subset**. They were not used in GMM clustering, FLARE refinement, or the 4,334-person donor reference.

Each synthetic individual is created by taking one phased haplotype from a validation individual from source population A and the other phased haplotype from a validation individual from source population B across chromosomes 1-22.

This produces an exact first-generation 50/50 admixed diploid individual.

### Priority contrasts

```text
EUR x Filipino/SEA
EUR x EAS
EUR x SAS
EUR x CAS
EUR x MID
EAS x SEA
SAS x CAS
```

The most important first benchmark is `EUR x Filipino`, because it directly tests whether K18 systematically creates a northern/central-Eurasian (`anc_6`) contribution in an individual known to be exactly 50% European and 50% Filipino.

### Metrics

For each synthetic individual:

1. expected broad ancestry
2. raw K18 latent-state vector
3. projected broad ancestry using mappings learned only from the training subset
4. off-target leakage outside the two known source groups
5. 50/50 balance error
6. total-variation distance from the known ancestry mixture

If synthetic EUR x Filipino individuals repeatedly show substantial `anc_6`, that demonstrates a systematic K18 limitation.

If synthetic EUR x Filipino individuals reconstruct approximately 50/50 EUR/SEA while Lonni alone carries substantial `anc_6`, Lonni's signal is more likely to reflect individual-specific haplotype structure rather than an intrinsic model artifact.

## 24. Final K18 multi-seed stability test

The leading configuration is fixed at:

```text
K = 18
covariance_type = diag
n_init = 10
max_iter = 300
reg_covar = 1e-6
```

The existing all-autosome model uses GMM seed `20261002`. Four additional GMM seeds are tested:

```text
42
1234
9999
314159
```

The FLARE seed is held fixed while the GMM seed changes. This isolates the stability of the ancestry-state discovery step rather than combining GMM and FLARE randomness.

Because GMM component labels are arbitrary, ancestry states from each new seed are aligned to the baseline K18 states using the ancestry proportions of the same untouched validation individuals. State matching is performed with a one-to-one assignment maximizing similarity across validation samples.

Primary stability metrics:

- mean matched-state similarity
- median matched-state similarity
- mean per-individual total-variation distance after state alignment
- 95th-percentile total-variation distance
- validation accuracy and balanced accuracy for each seed

The biological structure is considered stable if equivalent states recur across seeds even when their numeric `anc_N` labels change.

## 25. Latent-profile synthetic benchmark

The first synthetic benchmark projected every latent state to a single hard broad-group label. This is useful for simple contrasts such as EUR × Filipino, but it is too crude for populations such as Middle Eastern and Central Asian groups because these populations naturally occupy combinations of latent states.

The final synthetic benchmark therefore compares the **raw latent ancestry vectors directly**.

For each synthetic individual, two expected latent vectors are calculated.

### Donor-derived expectation

If the synthetic individual was generated from held-out validation donors A and B:

```text
expected = 0.5 × latent_profile(A) + 0.5 × latent_profile(B)
```

The synthetic genome uses one parental haplotype from each donor, so this diploid donor average is an expectation rather than an exact locus-by-locus target, but it provides a direct model-specific additivity benchmark.

### Population-derived expectation

A second expectation uses the mean latent profile of each parental fine population in the untouched validation set:

```text
expected =
    0.5 × mean_profile(fine_population_A)
  + 0.5 × mean_profile(fine_population_B)
```

This reduces noise caused by choosing one particular donor pair.

### Metrics

For every synthetic individual and contrast:

- total-variation distance to donor-derived expectation
- total-variation distance to population-derived expectation
- mean observed latent vector
- mean expected latent vector

This benchmark is permutation-invariant within each model because the observed and expected vectors use the same model's latent states.

The goal is no longer to force every latent state into exactly one broad geographic label. Instead, the model is tested on whether ancestry composition behaves approximately linearly when two known parental genomes are combined.

## 26. Production-model decision rule

K18 will be considered ready to move from model research to SNPster integration if it satisfies all of the following:

1. strong untouched-reference validation performance;
2. coherent behavior on independent known individuals;
3. low error on controlled synthetic admixture benchmarks;
4. no systematic EUR × Filipino to Central-Asian sink artifact;
5. reproducible latent population structure across GMM seeds.

After these checks, further hyperparameter searching should stop unless a specific reproducible failure mode remains.

---

## 27. Final ancestry ontology

The original broad labels were revised before the final model was trained. The final consumer-facing ancestry vector contains **eight additive source components**:

```text
EUR  European
EAS  East Asian
SAS  South Asian
SEA  South East Asian
GME  Greater Middle Eastern
NAM  Native American
OCE  Oceanian
SSA  Sub-Saharan African
```

These eight components are treated as additive source proportions and are constrained to sum to 1.0 (100%).

Several labels that are useful as population or cohort descriptors are **not** treated as additive ancestral sources:

```text
HLA  Hispanic / Latin American affinity
AFA  African American / Afro-Caribbean affinity
CAS  Central Asian population affinity
```

The reason is conceptual: HLA, AFA, and CAS describe admixed or clinal population groupings rather than independent ancestral source populations. In particular, Central Asian populations such as Kazakh, Uygur, Tajik, and Turkmen lie on combinations of western-Eurasian, South-Asian, East-Asian, and related axes. Treating CAS as one additive source recreated the original sink problem.

The production ancestry output therefore separates:

```text
additive source ancestry  -> EUR/EAS/SAS/SEA/GME/NAM/OCE/SSA
population affinity       -> CAS/HLA/AFA and future affinity labels
```

The additive ancestry vector is intended for PGS applicability weighting. Affinity labels are separate contextual evidence and do not have to sum to 100%.

## 28. Source-only training and validation split

The 552 held-out reference samples were split by fine population into a model-training subset and an untouched validation subset using seed `20261002`.

Overall split:

```text
training     365
validation   187
```

After separating additive-source samples from affinity-only samples:

```text
source training       307
affinity training      58
source validation     159
affinity validation    28
```

Source-training counts:

| Source | N |
|---|---:|
| EAS | 63 |
| SAS | 62 |
| EUR | 58 |
| SSA | 48 |
| SEA | 39 |
| GME | 22 |
| NAM | 9 |
| OCE | 6 |

The production source model uses only the eight additive-source groups. Affinity-only samples remain available for separate diagnostics.

Important files:

```text
/output/source_model/source_reference.autosomes.vcf.gz
/output/source_model/source_train.autosomes.vcf.gz
/output/source_model/source_validation.autosomes.vcf.gz
/output/source_model/source_reference.ref.panel
/output/source_model/source_train_metadata.tsv
/output/source_model/source_validation_metadata.tsv
```

## 29. Final latent model selection

Source-only diagonal Gaussian-mixture models were compared at K18, K24, and K30. Each model was trained on the 307 source-training individuals and refined with FLARE using:

```text
update-p=true
em=true
array=false
```

Validation was always frozen:

```text
update-p=false
em=false
probs=false
array=false
seed=20261002
```

### Pure-source validation

On the 159 untouched source-validation individuals:

| Model | 3NN accuracy | 3NN balanced accuracy |
|---|---:|---:|
| K18 | 93.08% | 94.23% |
| K24 | 96.23% | 97.34% |
| K30 | 97.48% | 98.19% |

K30 was the best hard pure-source classifier, but hard classification was not the production objective. SNPster requires quantitatively meaningful ancestry mixtures.

### Synthetic latent reconstruction

A controlled benchmark of 336 FLARE-processed mosaics was generated from untouched validation individuals:

```text
2-way mixtures   84
3-way mixtures  112
4-way mixtures  140
total           336
```

The test compared the observed FLARE latent profile with the donor-weighted expected latent profile using total-variation distance.

| Model | Mean latent TV |
|---|---:|
| K18 | **4.34%** |
| K24 | 4.96% |
| K30 | 4.58% |

K18 showed the best latent additivity overall and was therefore retained for ancestry-percentage decoding.

### Final latent-model decision

The production latent model is:

```text
K = 18
covariance_type = diag
n_init = 10
max_iter = 300
reg_covar = 1e-6
GMM seed = 20261002
```

Model path:

```text
/output/source_model/models/diag_K18/train.final.model
```

The Bayesian branch was discontinued. The diagonal GMM is easier to interpret and already showed strong, stable behavior on the relevant source-decomposition task.

## 30. Latent-to-source decoder development

FLARE's K18 output is an 18-dimensional latent ancestry profile. A second-stage decoder maps that latent profile to the eight source proportions.

Several decoder families were evaluated:

- distance-weighted / nearest-neighbor decoders;
- multi-output Ridge regression;
- ExtraTrees regression;
- a small multi-layer perceptron.

The initial Ridge decoder reconstructed mixtures reasonably well but performed poorly on pure-source endpoints. Pure-sample weighting did not solve this, indicating that a single global linear mapping was too restrictive.

### Nonlinear decoder benchmark

ExtraTrees provided the best balance between pure-source and mixture reconstruction. A useful early configuration was:

```text
min_samples_leaf = 5
pure_weight = 50
threshold = 5%
```

with approximately:

```text
pure TV       10.26%
mixture TV    11.28%
off-target     6.27%
```

However, the largest improvement came from calibrating the decoder with **actual FLARE-processed mosaics**, rather than only idealized linear combinations of latent profiles.

## 31. Actual FLARE-mosaic calibration

Three independent training-mosaic sets were generated from the 307 training donors:

```text
336 mosaics per run
3 runs
1008 total actual FLARE training mosaics
```

These mosaics were processed through the frozen K18 FLARE model. They were then added to the ExtraTrees training set with explicit sample weighting.

The final base decoder configuration is:

```text
ExtraTreesRegressor
n_estimators = 200
min_samples_leaf = 10
max_features = 1.0
pure_weight = 50
actual_mosaic_weight = 100
random_state = 20261005
```

The actual-mosaic calibration strongly improved ancestry reconstruction. Across the decoder sweep, increasing actual-mosaic weight reduced average synthetic TV from approximately 11.3% without actual mosaics to approximately 7.1% at high actual-mosaic weight, while pure-source TV remained approximately stable.

Frozen base-decoder path:

```text
/output/source_model/crossfit_calibration/k18_extratrees_base.joblib
```

## 32. Global linear calibration

The raw ExtraTrees output still contained systematic cross-loading between ancestry sources, particularly EAS and SEA. A final linear calibration layer was therefore learned using out-of-fold predictions from the three broad training-mosaic runs.

The calibrator is a multi-output Ridge regression followed by projection back onto the probability simplex.

Final configuration:

```text
Ridge
alpha = 0.0001
fit_intercept = true
```

Frozen calibrator path:

```text
/output/source_model/crossfit_calibration/k18_linear_calibrator.joblib
```

On the development validation set the calibration changed performance approximately as follows:

| Metric | Raw ExtraTrees | Globally calibrated |
|---|---:|---:|
| Pure-source TV | 11.69% | **8.75%** |
| Broad synthetic TV | 10.22% | **6.58%** |
| Broad synthetic p95 TV | 19.16% | **15.06%** |
| Off-target mass | 8.49% | **3.15%** |
| Dominant-source accuracy | 93.08% | **93.08%** |

A separate rare-aware calibrator was investigated but rejected. Giving 1-10% matched ancestry examples more weight improved low-level incremental recovery only modestly while consistently worsening the main global ancestry metrics. The globally calibrated model is therefore the production output.

## 33. Low-level ancestry and negative-control experiments

Low-level ancestry was studied explicitly using matched backbones. For each family, an intact validation individual was used as the 0% baseline and 1%, 2%, 5%, or 10% of the genome was replaced by tracts from a donor belonging to another source group.

The raw decoder showed strong within-family monotonicity:

```text
exactly monotonic families       91.7%
monotonic within +/-0.5%         98.2%
mean within-family correlation    0.981
median within-family correlation  0.992
```

The raw decoder systematically under-recovered small ancestry components on average:

| True inserted ancestry | Mean incremental raw response | Mean recovery |
|---:|---:|---:|
| 1% | 0.42% | 40.7% |
| 2% | 0.93% | 46.2% |
| 5% | 3.03% | 60.8% |
| 10% | 7.15% | 71.7% |

The global Ridge calibration intentionally shrinks small components further, which improves total global ancestry accuracy but reduces sensitivity to trace ancestry.

### Same-source negative control

A second benchmark replaced 1%, 2%, 5%, or 10% of a person's genome with DNA from a *different individual belonging to the same source group*. The true ancestry therefore did not change. This measured how much decoder components can move merely because of within-source donor heterogeneity.

For the raw decoder, source-specific discrimination of true cross-source ancestry from same-source perturbation was:

| True inserted ancestry | Component AUC | Sensitivity above 95th-percentile component null |
|---:|---:|---:|
| 1% | 0.929 | 78.0% |
| 2% | 0.955 | 86.3% |
| 5% | 0.973 | 94.1% |
| 10% | 0.981 | 97.0% |

Pair-specific 95th-percentile sensitivity was approximately:

```text
1%   81.0%
2%   88.7%
5%   92.9%
10%  96.4%
```

The much stricter family-wise maximum-null statistic was not used as the primary quantitative resolution criterion because it compares the correct minor source against the largest movement among all seven alternative ancestry components and becomes increasingly conservative as the replaced fraction grows.

Practical interpretation:

```text
<2%    trace signal; below reliable quantitative resolution
2-5%   low-confidence minor ancestry
>=5%   generally meaningful source component
>=10%  strong source component
```

No hard cutoff is applied internally. Continuous calibrated proportions are retained for PGS applicability calculations.

## 34. Final untouched synthetic generalization test

After the model architecture was frozen, one final synthetic test was generated with a new seed and ancestry templates that had not been used for development.

Final-test templates:

```text
2-way:
65 / 35
85 / 15
95 / 5

3-way:
50 / 35 / 15
70 / 20 / 10

4-way:
40 / 30 / 20 / 10
55 / 25 / 15 / 5
```

The benchmark contained **672 newly generated genomes**. The frozen K18 model, frozen ExtraTrees decoder, and frozen Ridge calibrator were applied without refitting.

### Pre-declared acceptance criteria

Before examining final results, the acceptance criteria were fixed as:

```text
calibrated mean TV     <= 8%
calibrated p95 TV      <= 18%
off-target ancestry    <= 5%
```

Stronger optional targets were:

```text
mean TV                <= 7%
p95 TV                 <= 16%
```

### Final overall result

| Metric | Raw decoder | Calibrated production output |
|---|---:|---:|
| Mean TV | 12.37% | **7.90%** |
| Median TV | 12.05% | **7.30%** |
| p95 TV | 21.72% | **15.98%** |
| Maximum TV | 33.17% | **29.76%** |
| Dominant-source accuracy | 96.58% | **98.21%** |
| Off-target ancestry mass | 10.06% | **4.21%** |

Calibration improved the individual-level TV error for **90.6%** of the final-test genomes. Mean TV improved by approximately **4.47 percentage points**.

### Performance by number of ancestral sources

| Sources | N | Mean TV | Median TV | p95 TV | Off-target | Dominant accuracy |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 168 | **5.29%** | 4.26% | 13.16% | 2.69% | 99.41% |
| 3 | 224 | **7.93%** | 7.37% | 15.09% | 4.91% | 98.21% |
| 4 | 280 | **9.45%** | 8.61% | 18.62% | 4.56% | 97.50% |

Four-way mixtures are predictably harder than simple two-way mixtures. This should be considered when interpreting uncertainty.

### Performance by mixture template

| Template | Mean TV | p95 TV | Off-target |
|---|---:|---:|---:|
| 65 / 35 | 5.92% | 15.22% | 2.82% |
| 85 / 15 | 5.75% | 13.62% | 2.53% |
| 95 / 5 | **4.21%** | **9.84%** | 2.72% |
| 50 / 35 / 15 | 8.43% | 16.70% | 5.17% |
| 70 / 20 / 10 | 7.42% | 12.92% | 4.66% |
| 40 / 30 / 20 / 10 | 9.62% | 18.11% | 4.55% |
| 55 / 25 / 15 / 5 | 9.28% | 19.03% | 4.58% |

### Source-specific MAE

| Source | Calibrated MAE | Mean bias |
|---|---:|---:|
| EUR | 1.69% | +0.11% |
| EAS | 3.14% | -0.55% |
| SAS | 1.62% | -0.53% |
| SEA | 3.21% | +1.42% |
| GME | 1.82% | +0.42% |
| NAM | 1.24% | +0.60% |
| OCE | 1.43% | -0.49% |
| SSA | 1.65% | -0.98% |

EAS and SEA remain the hardest source pair to partition precisely. This is a documented limitation rather than a reason to continue tuning against the final test set.

### Final pass/fail

```text
Mean TV <= 8%:       PASS   7.902%
P95 TV <= 18%:       PASS  15.975%
Off-target <= 5%:    PASS   4.212%

Strong mean <= 7%:   NO
Strong p95 <= 16%:   YES

FINAL STATUS: PASS
```

No model changes are made after observing these final-test results.

## 35. Metric definitions

This section defines the metrics used throughout the model-development and validation experiments.

### 35.1 Total-variation distance (TV)

For a true ancestry vector `p` and predicted ancestry vector `q`, both summing to 1:

```text
TV(p, q) = 0.5 * sum_i |p_i - q_i|
```

TV lies between 0 and 1.

```text
TV = 0      -> prediction exactly matches truth
TV = 0.05   -> 5% of total ancestry probability mass is assigned incorrectly
TV = 0.10   -> 10% of total ancestry probability mass is assigned incorrectly
```

Example:

```text
truth:
EUR  50%
SAS  30%
GME  20%

prediction:
EUR  55%
SAS  27%
GME  16%
EAS   2%
```

Absolute differences are:

```text
EUR  5%
SAS  3%
GME  4%
EAS  2%
-------------
sum 14%
```

Therefore:

```text
TV = 14% / 2 = 7%
```

The factor of one-half prevents counting the same misplaced ancestry mass twice: mass removed from one source must appear in another source.

TV is the primary ancestry-composition error metric because it measures the whole probability vector at once and is easy to interpret as misplaced ancestry mass.

### 35.2 Mean TV

Mean TV is the arithmetic average of individual TV errors:

```text
mean_TV = mean(TV_1, TV_2, ..., TV_n)
```

It describes typical aggregate ancestry-composition error across the evaluation cohort.

### 35.3 Median TV

Median TV is the 50th percentile of individual TV errors. Half of individuals have lower TV and half have higher TV. It is less sensitive to a small number of very difficult individuals than mean TV.

### 35.4 p95 TV

The 95th-percentile TV is the error level below which 95% of tested individuals fall.

Example:

```text
p95 TV = 0.16
```

means that 95% of test individuals have TV <= 16%, while the worst 5% have larger errors.

This is useful for describing the error tail rather than only average behavior.

### 35.5 Maximum TV

Maximum TV is the largest TV error observed in the benchmark. It identifies the single worst tested individual but is highly sensitive to one outlier, so it is not used alone for model selection.

### 35.6 Source-specific mean absolute error (MAE)

For source `s`:

```text
MAE_s = mean(|predicted_s - true_s|)
```

Example:

```text
true SEA      40%
predicted SEA 34%
absolute error 6 percentage points
```

MAE measures how precisely one particular ancestry source is estimated. Unlike TV, it does not consider the complete ancestry vector simultaneously.

### 35.7 Bias

For source `s`:

```text
bias_s = mean(predicted_s - true_s)
```

Interpretation:

```text
positive bias -> source is systematically overestimated
negative bias -> source is systematically underestimated
bias near 0   -> no strong directional error on average
```

A low bias does not necessarily imply low error because over- and underestimates can cancel. MAE should therefore be considered alongside bias.

### 35.8 Dominant-source accuracy

For each individual:

```text
true dominant source      = argmax(true ancestry vector)
predicted dominant source = argmax(predicted ancestry vector)
```

Dominant-source accuracy is the fraction of individuals for whom these two labels agree.

It measures broad classification correctness but does not measure whether the ancestry percentages themselves are accurate. For that reason TV remains the primary quantitative metric.

### 35.9 Balanced accuracy

Balanced accuracy is the mean recall across ancestry classes:

```text
balanced_accuracy = mean(class-specific recall)
```

It gives each source class equal importance even when the validation sample counts are unequal. This was useful during latent-model selection because rare classes such as OCE and NAM otherwise contribute little to ordinary accuracy.

### 35.10 Off-target ancestry mass

In a synthetic genome, the true contributing ancestry sources are known exactly. Off-target mass is the total predicted probability assigned to sources whose true proportion is zero:

```text
off_target = sum(predicted_i for all i where true_i = 0)
```

Example:

```text
truth:
EUR 60%
SEA 40%
all others 0%

prediction:
EUR 55%
SEA 37%
EAS  5%
SAS  3%
```

Then:

```text
off-target = EAS + SAS = 8%
```

This metric is especially useful for detecting spurious ancestry components.

### 35.11 Recovery fraction

For matched low-level ancestry experiments:

```text
recovery = observed incremental ancestry / true inserted ancestry
```

Example:

```text
true inserted ancestry = 5%
model increase          = 3%
recovery                = 3 / 5 = 60%
```

Recovery measures sensitivity of the decoder to incremental ancestry. It is not used as the primary production composition metric because a real user does not have an observable 0%-ancestry counterfactual baseline.

### 35.12 AUC in the same-source negative-control benchmark

AUC measures how well a true cross-source ancestry signal is ranked above changes caused by replacing DNA with another donor from the same ancestry source.

```text
AUC = 0.5   -> no discrimination
AUC = 1.0   -> perfect discrimination
```

The raw decoder achieved component AUC values of approximately 0.93, 0.96, 0.97, and 0.98 for true 1%, 2%, 5%, and 10% inserted ancestry respectively.

### 35.13 Component-null sensitivity

For each same-source replacement fraction, a null distribution is built from changes in ancestry components when the true source does not change.

`Sensitivity above 95th-percentile component null` is:

```text
fraction of true cross-source test cases
where the correct minor-source increase
exceeds the 95th percentile of same-source component changes
```

This gives an empirical detection measure for small ancestry components without assuming a parametric noise distribution.

### 35.14 Family-wise maximum-null sensitivity

A stricter null statistic compares the correct ancestry signal against the **largest** false ancestry movement among all alternative source components in each same-source genome.

This controls a much broader discovery question and is therefore substantially more conservative than the component-specific null. It was retained as a diagnostic but was not used as the primary quantitative ancestry-resolution criterion.

### 35.15 Latent reconstruction TV

Before source decoding, the K18 model itself was evaluated in latent-state space.

For a synthetic mosaic assembled from known validation donors, the expected latent vector is the genetic-length-weighted average of the donor latent profiles. FLARE's observed latent vector is then compared to this expectation with TV.

K18 achieved approximately:

```text
mean latent reconstruction TV = 4.34%
```

This showed that much of the remaining final ancestry error came from latent-to-source decoding/calibration rather than from FLARE's latent representation itself.

## 36. Frozen production architecture

As of 2026-10-05, the SNPster ancestry model is frozen at:

```text
Genotype / imputed autosomal VCF
        |
        v
FLARE 0.6.0
K18 diagonal GMM
update-p=false
em=false
        |
        v
18-dimensional latent ancestry vector
        |
        v
ExtraTrees source decoder
200 trees
min_samples_leaf=10
pure_weight=50
actual_mosaic_weight=100
        |
        v
8 raw source proportions
        |
        v
Ridge global calibration
alpha=0.0001
        |
        v
simplex projection
        |
        v
EUR / EAS / SAS / SEA / GME / NAM / OCE / SSA
```

Frozen artifacts:

```text
FLARE model:
/output/source_model/models/diag_K18/train.final.model

ExtraTrees decoder:
/output/source_model/crossfit_calibration/k18_extratrees_base.joblib

Ridge calibrator:
/output/source_model/crossfit_calibration/k18_linear_calibrator.joblib

Reference:
/output/source_model/source_reference.autosomes.vcf.gz
/output/source_model/source_reference.ref.panel

Map:
/output/model_sweep_allchr/plink.autosomes.GRCh38.map
```

For new users FLARE must remain frozen:

```text
update-p=false
em=false
probs=false
array=false
```

No development benchmark should be used for additional hyperparameter tuning unless a new reproducible model failure is identified.

## 37. Production interpretation for SNPster

The calibrated vector is the production ancestry vector and should be used for PGS applicability weighting.

Store both calibrated and raw values if convenient:

```text
raw_EUR ... raw_SSA
cal_EUR ... cal_SSA
```

The raw vector is useful for diagnostics and low-level ancestry sensitivity analyses. The calibrated vector has substantially better overall ancestry-composition accuracy and lower off-target leakage and is therefore the consumer-facing / PGS-facing estimate.

Do **not** hard-zero all values below 5%. The calibrated vector should remain continuous and sum to 1.0.

Suggested interpretation language:

```text
<2%    trace / below reliable quantitative resolution
2-5%   low-confidence minor component
>=5%   generally meaningful source component
>=10%  strong source component
```

These are interpretation bands, not hard computational thresholds.

### PGS usage

The ancestry vector should be used to evaluate how relevant ancestry-specific PGS evaluation evidence is for the user. It should not be used to naively blend unrelated published polygenic scores as, for example:

```text
0.7 * PGS_EUR + 0.3 * PGS_EAS
```

unless the underlying scoring method explicitly supports ancestry- or local-ancestry-specific combination. The primary role of the ancestry vector is applicability, confidence, normalization context, and evidence weighting.

## 38. Final production status

The ancestry model-development phase is complete.

The frozen model passed all pre-declared final-test acceptance criteria:

```text
mean TV <= 8%       PASS
p95 TV <= 18%       PASS
off-target <= 5%    PASS
```

The final calibrated test performance was:

```text
mean TV              7.90%
median TV            7.30%
p95 TV              15.98%
dominant accuracy   98.21%
off-target           4.21%
```

Known remaining limitations:

1. EAS and SEA are harder to partition than the other source groups.
2. Four-way mixtures have wider errors than two-way mixtures.
3. Components below approximately 2% should not be interpreted quantitatively.
4. The final synthetic benchmark used untouched synthetic realizations but the donor pool had previously been used during model-development validation; therefore this is a strong synthetic generalization test rather than a fully external donor-independent validation cohort.
5. Population-affinity outputs such as CAS remain separate from the additive source vector and require their own validated reporting layer.

The next work is **integration and real-user sanity checking**, not further hyperparameter search.

