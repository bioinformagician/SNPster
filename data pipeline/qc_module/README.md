# Imputation References and Beagle Output QC

## References in the pipeline

The workflow uses two distinct reference sets:

- **Beagle imputation:** `imputation_module` selects a reference panel and a PLINK genetic map for the chromosome in the merged, harmonized VCF. In the imputer container these are mounted under `/data/beagle_references` and `/data/plink.GRCh38.map`. Reference panels may be `.vcf.gz`, `.bref3`, or `.bcf`; map files are `.map`. Both chromosome-specific filenames must contain a matching chromosome token such as `chr1.` for the reference or `chr1...map` for the map. If either chromosome-specific asset is missing, imputation stops with a mapping error. Beagle is run with `gt`, `ref`, `map`, and `gp=true`, producing genotype probabilities (`GP`) in its compressed VCF output.
- **Ancestry (separate from imputation):** the orchestration mounts the HGDP+1kGP reference VCFs and the 10-population panel for ancestry estimation. It uses the configured reference filename pattern and chromosomes 1, 2, 21, and 22. These are not the Beagle imputation references.

## QC flow and requirements

Nextflow splits the Beagle VCF by chromosome, then runs `qc_module` on each split. The QC reader expects one sample column with `GT:DS:GP` fields and three diploid genotype probabilities in `GP` order `0/0,0/1,1/1`. Rows with missing or unparsable required values fail the numeric comparisons.

The default thresholds are `GP_MIN=0.90` and `DS_TOL=0.05`. A variant must satisfy both the maximum of its three GP values being at least `GP_MIN` and a genotype/dosage consistency rule:

- `0/0` or `0|0`: `DS <= DS_TOL`
- heterozygous (`0/1`, `1/0`, or phased equivalents): `GP_01 >= GP_MIN` and `DS` within `1 +/- DS_TOL`
- `1/1` or `1|1`: `DS >= 2 - DS_TOL`

SNP-only and biallelic-only filters are also available. When enabled, they require single-base REF/ALT alleles and no comma-separated ALT alleles, respectively. The QC module’s standalone defaults enable both, but the integrated `imputation_qc` Docker image sets `SNPS_ONLY=False` and `BIALLELIC_ONLY=False`; therefore the current Nextflow pipeline does not apply those two filters unless the container environment is changed.

`DR2` and `AF` are read from Beagle INFO fields and preserved in the QC output when present, but no DR2 or AF cutoff is currently applied. An empty input or a file with no passing variants fails QC. Passing variants are written as a compressed, tabix-indexed, single-sample VCF named with `_QC.vcf.gz`; its sample column is emitted as `GT:DS:GP`.