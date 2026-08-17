# HGDP+1KG Reference Panel Setup

This directory contains tools for downloading and preparing the HGDP+1000 Genomes reference panel for imputation and ancestry inference.

## Overview

The HGDP+1KG reference panel combines:
- **HGDP (Human Genome Diversity Project)**: ~900 samples from diverse global populations
- **1000 Genomes Project**: ~3,200 samples from 26 populations

This combined panel provides better representation across global populations for imputation and ancestry estimation.

## Files

- `prepare_hgdp_1kg_reference_panel.nf` - Nextflow workflow to download and prepare the reference panel
- `convert_hgdp_1kg_metadata_to_panel.py` - Python script to convert metadata to panel format

## Quick Start

### 1. Download and Prepare Reference Panel (Nextflow)

```bash
cd "/home/frederik/github_projects/SNPster/data pipeline/util_dev"

# Download raw BCF files and metadata, then prepare VCF files
nextflow run prepare_hgdp_1kg_reference_panel.nf \
  --download_raw true \
  --input_dir /srv/dependencies/imputation_runner/imputer/beagle_references_hgdp \
  --output_dir /srv/dependencies/imputation_runner/imputer/beagle_references_hgdp/prepared_vcf \
  --panel_output_dir /srv/dependencies/imputation_runner

# Or if files are already downloaded, just prepare them
nextflow run prepare_hgdp_1kg_reference_panel.nf \
  --download_raw false \
  --input_dir /srv/dependencies/imputation_runner/imputer/beagle_references_hgdp
```

### 2. Manual Panel File Creation (Python)

If you need to regenerate just the panel file:

```bash
cd "/home/frederik/github_projects/SNPster/data pipeline/util_dev"
source /home/frederik/github_projects/SNPster/.venv-1/bin/activate

# Convert metadata to panel format
python convert_hgdp_1kg_metadata_to_panel.py \
  /srv/dependencies/imputation_runner/imputer/beagle_references_hgdp/gnomad_meta_updated.tsv \
  /tmp/hgdp_1kg_panel.txt

# Copy to the correct location (requires sudo)
sudo cp /tmp/hgdp_1kg_panel.txt /srv/dependencies/imputation_runner/
```

## Panel File Format

The generated panel file follows the standard format expected by ancestry inference tools:

```
sample      pop         super_pop   gender
HGDP00001   Brahui      SAS         male
HGDP00003   Brahui      SAS         male
HG00096     GBR         EUR         male
```

### Columns:
- **sample**: Sample ID (e.g., HGDP00001, HG00096)
- **pop**: Detailed population code (e.g., Brahui, GBR, YRI)
- **super_pop**: Superpopulation (EUR, AFR, EAS, SAS, AMR)
- **gender**: male/female/unknown

## Superpopulation Mapping

The script uses the **standard 5-population model** compatible with **PGS-Calc** and most polygenic scoring resources:

| Code | Superpopulation | Description |
|------|----------------|-------------|
| EUR  | European       | Non-Finnish European, Finnish, **Middle Eastern** |
| AFR  | African        | African populations |
| EAS  | East Asian     | East Asian, Oceania |
| SAS  | South Asian    | South/Central Asian populations |
| AMR  | American       | American/Latino populations |

### Why Middle East → EUR?

Middle Eastern populations map to **EUR** (not SAS) because:
1. **Genetic clustering**: Middle Eastern populations cluster with Europeans in PCA space (West Eurasia)
2. **PGS-Calc compatibility**: This matches the 1000 Genomes 5-population model used by PGS-Calc
3. **PGS performance**: Polygenic scores trained on European cohorts generally perform better in Middle Eastern populations than South Asian ones
4. **Standard convention**: This is the widely accepted mapping in population genetics research

**Note**: If you need a separate Middle Eastern category (6 or 7-population model), you can uncomment the alternative mappings in the Python script.

### Detailed Mappings:

**From gnomAD population inference:**
- `nfe`, `fin`, `mid` → EUR (European)
- `afr` → AFR (African)
- `eas` → EAS (East Asian)
- `sas` → SAS (South Asian)
- `amr` → AMR (American)
- `oth` → Excluded from panel (ambiguous ancestry)

**From HGDP/TGP genetic regions:**
- `EUR`, `MID` → EUR (European/Middle Eastern)
- `AFR` → AFR (African)
- `EAS`, `OCE` → EAS (East Asian/Oceania)
- `CSA` → SAS (Central/South Asian)
- `AMR` → AMR (American)

## Sample Counts (Expected)

After filtering:
- **Total**: ~4,040 samples
- **AFR**: ~996 samples
- **EUR**: ~901 samples (including Middle Eastern)
- **EAS**: ~808 samples
- **SAS**: ~788 samples
- **AMR**: ~547 samples

## Data Sources

- **Reference VCF files**: 
  - URL: `https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/`
  - Format: BCF (Binary VCF)
  - Chromosomes: 1-22 (autosomes)

- **Metadata file**:
  - URL: `https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv`
  - Contains: Sample IDs, populations, sex, quality metrics

## VCF Preparation

The workflow processes BCF files to:
1. Convert chromosome names from `chr1` to `1` (numeric)
2. Keep only **biallelic SNPs** (removes indels and multiallelic sites)
3. Create indexed VCF.gz files suitable for imputation

## Usage in Ancestry Module

To use this panel with the ancestry module, update the configuration:

```python
# data pipeline/ancestry_module/config.py

REFERENCE_VCF_DIR = "/srv/dependencies/imputation_runner/imputer/beagle_references_hgdp/prepared_vcf"
REFERENCE_VCF_PATTERN = "hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"
POPULATION_PANEL_FILE = "/srv/dependencies/imputation_runner/hgdp_1kg_panel.txt"
REFERENCE_PANEL = "HGDP+1KG"
K_POPULATIONS = 5  # EUR, AFR, EAS, SAS, AMR
```

## Compatibility with PGS-Calc

This panel format and superpopulation mapping is **fully compatible with PGS-Calc** (pgscatalog/pgsc_calc):

- Uses the same HGDP+1KG reference data (gnomAD v3.1)
- Maps to the standard 5 superpopulations: EUR, AFR, EAS, SAS, AMR
- Compatible with PGS-Calc's `--run_ancestry` feature for ancestry adjustment
- Can be used as a custom reference panel for PGS-Calc

See [PGS-Calc documentation on genetic ancestry](https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html) for more details.

## Troubleshooting

### Panel file has missing gender values

If you see many "unknown" gender values, check which sex column is being used:
- Preferred: `project_meta.sex` (most complete)
- Alternative: `sex_imputation.sex_karyotype` (inferred from genotypes)
- Fallback: `sex` or `bergstrom.sex` (may have missing values)

### Metadata conversion fails

Ensure pandas is installed in your Python environment:
```bash
source /home/frederik/github_projects/SNPster/.venv-1/bin/activate
pip install pandas
```

### Permission denied when copying files

The `/srv/dependencies/` directory requires sudo access. Copy to `/tmp/` first:
```bash
python convert_hgdp_1kg_metadata_to_panel.py metadata.tsv /tmp/panel.txt
sudo cp /tmp/panel.txt /srv/dependencies/imputation_runner/
```

## References

- **gnomAD HGDP+1KG dataset**: https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
- **HGDP**: Bergström et al. (2020) "Insights into human genetic variation and population history from 929 diverse genomes" Science
- **1000 Genomes**: The 1000 Genomes Project Consortium (2015) Nature

## Notes

- The panel excludes samples with ambiguous ancestry (`OTH` superpopulation)
- Middle Eastern populations are grouped with European for the 5-population model
- Oceania populations are grouped with East Asian for proximity
- Sex-chromosome variants are excluded (autosomes only, chr 1-22)
