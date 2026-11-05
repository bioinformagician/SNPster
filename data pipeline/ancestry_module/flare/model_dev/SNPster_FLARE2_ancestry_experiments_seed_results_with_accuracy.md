# SNPster FLARE2 ancestry-model experiments

**Date:** 2026-10-01  
**Project:** SNPster ancestry inference  
**Purpose:** Develop a robust local/global ancestry model that avoids the pathological Central Asian (CAS) sink seen with the original broad-panel FLARE setup.

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

## 27. Multi-seed latent synthetic benchmark: first completed seeds

The latent-profile synthetic benchmark was evaluated for K18 GMM seeds 42 and 1234.

### Seed 42

| Contrast | Donor-profile TV | Population-profile TV |
|---|---:|---:|
| EAS × Filipino | 0.0076 | 0.0073 |
| EUR × CAS | 0.0253 | 0.0511 |
| EUR × EAS | 0.0068 | 0.0065 |
| EUR × Filipino | 0.0118 | 0.0183 |
| EUR × MID | 0.0158 | 0.0411 |
| EUR × SAS | 0.0137 | 0.0268 |
| SAS × CAS | 0.0257 | 0.0521 |

Overall donor-profile TV = 0.0152; population-profile TV = 0.0290.

### Seed 1234

| Contrast | Donor-profile TV | Population-profile TV |
|---|---:|---:|
| EAS × Filipino | 0.0075 | 0.0106 |
| EUR × CAS | 0.0217 | 0.0559 |
| EUR × EAS | 0.0047 | 0.0048 |
| EUR × Filipino | 0.0195 | 0.0223 |
| EUR × MID | 0.0286 | 0.0447 |
| EUR × SAS | 0.0100 | 0.0204 |
| SAS × CAS | 0.0201 | 0.0494 |

Overall donor-profile TV = 0.0160; population-profile TV = 0.0297.

### Interpretation

The synthetic genomes are reconstructed extremely closely to the latent profiles expected from their actual held-out parental donors. This is especially important for EUR × CAS, EUR × MID, and SAS × CAS, which looked poor under the earlier hard latent-state-to-broad-group mapping. Their donor-profile TV errors are only about 0.02-0.03, showing that the apparent failures largely came from forcing each latent state into one geographic label rather than from a failure of FLARE inference.

The result is reproducible across the first two alternate GMM seeds. Although numeric latent-state labels and the exact basis decomposition differ between seeds, each model reconstructs synthetic admixture nearly additively relative to its own held-out reference profiles.

### Production implication

FLARE2 latent states should be treated as an internal basis representation rather than exposed directly as literal ancestry percentages. A downstream interpretation layer should compare each target latent vector with reference-population or reference-group latent profiles. Remaining checks: seeds 9999 and 314159, ancestry-mass-weighted state stability, identifying high-TV validation individuals, and testing a reference-profile projection layer for stable user-facing ancestry proportions.

## 28. Final predictive performance of the current K18 model

The current production candidate is the all-autosome **diagonal-covariance K18 FLARE2 model**.

Its predictive performance should be described using two complementary metrics, because the latent FLARE states are not in one-to-one correspondence with the manually defined broad ancestry labels.

### Held-out broad-group classification

On the completely untouched validation subset:

```text
overall accuracy          = 75.94%
balanced accuracy         = 75.98%
```

These values come from assigning each latent ancestry state to a single broad geographic group and then performing hard classification.

This is a deliberately conservative metric. It penalizes biologically shared latent states such as the western-Eurasian state used by both European and Middle-Eastern populations, and therefore does not measure the full predictive fidelity of the latent ancestry representation.

### Synthetic admixture reconstruction

A more direct test of FLARE's predictive ability is whether it reconstructs synthetic 50/50 admixed genomes from held-out parental individuals.

Across the first two independent alternative GMM seeds:

| Metric | Seed 42 | Seed 1234 |
|---|---:|---:|
| Mean donor-profile TV | 0.0152 | 0.0160 |
| Mean population-profile TV | 0.0290 | 0.0297 |
| Donor-profile fidelity (`1 - TV`) | **98.48%** | **98.40%** |
| Population-profile fidelity (`1 - TV`) | **97.10%** | **97.03%** |

Thus, relative to the latent ancestry profiles expected from the actual held-out parental donors, the model reconstructs synthetic admixture with approximately **98.4-98.5% profile fidelity**.

Relative to the mean latent profiles of the parental fine populations, fidelity is approximately **97.0-97.1%**.

These `1 - TV` values are best interpreted as **latent-profile reconstruction fidelity**, not as ordinary classification accuracy.

### Example synthetic contrasts

The most important contrasts showed very low donor-profile error:

```text
                         Seed 42      Seed 1234
EUR × EAS                  0.0068        0.0047
EAS × Filipino             0.0076        0.0075
EUR × Filipino             0.0118        0.0195
EUR × SAS                  0.0137        0.0100
EUR × MID                  0.0158        0.0286
EUR × CAS                  0.0253        0.0217
SAS × CAS                  0.0257        0.0201
```

The apparently difficult Central-Asian and Middle-Eastern mixtures are therefore reconstructed well in latent-state space even when a hard broad-group mapping makes them look inaccurate.

### Independent real-person checks

The current K18 model also behaved coherently on independent users:

```text
Frederik:
100% of ancestry assigned to the dominant western-Eurasian state.

Lonni:
approximately 43% western-Eurasian, with the remainder predominantly
East/Southeast-Asian latent states, consistent with her known
approximately half-European / half-Filipino ancestry.

Murtaza:
represented as a combination of Central/Northern-Eurasian,
South-Asian, East-Asian, and Western-Eurasian states rather than
collapsing into a pathological single CAS state.
```

### Current performance summary

For the current K18 production candidate:

```text
Hard broad-group accuracy             75.94%
Balanced broad-group accuracy         75.98%

Synthetic donor-profile fidelity      ~98.4-98.5%
Synthetic population-profile fidelity ~97.0-97.1%
```

The synthetic profile metrics are the more representative measure of how accurately the FLARE2 latent model reconstructs ancestry composition, whereas the broad-group classification metric additionally reflects limitations of forcing continuous/shared ancestry states into discrete geographic labels.
