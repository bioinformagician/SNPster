# Population Mapping Rationale for HGDP+1KG Panel

## Summary

The HGDP+1KG panel uses the **5-population model** (EUR, AFR, EAS, SAS, AMR) that is standard for:
- PGS-Calc (pgscatalog/pgsc_calc)
- 1000 Genomes Project
- Most polygenic scoring resources

## Mapping Details

### Middle Eastern → EUR (European)

**Rationale:**
1. **Genetic clustering**: Middle Eastern populations cluster with Europeans in principal component analysis (PCA), forming "West Eurasia"
2. **Shared ancestry**: Europeans and Middle Easterners share more recent common ancestry than Middle Easterners and South Asians
3. **Migration history**: Ancient Near Eastern farmers migrated to Europe, contributing significantly to European ancestry
4. **LD patterns**: Linkage disequilibrium patterns are more similar between Middle Eastern and European populations
5. **PGS performance**: Polygenic scores trained on European datasets show better transferability to Middle Eastern populations than to South Asian populations

**Scientific basis:**
- Lazaridis et al. (2014) Nature: "Ancient human genomes suggest three ancestral populations for present-day Europeans"
- Hofmanová et al. (2016) PNAS: "Early farmers from across Europe directly descended from Neolithic Aegeans"
- Population structure studies consistently show Middle East clustering with Europe in West Eurasia

### Oceania → EAS (East Asian)

**Rationale:**
1. **Geographic proximity**: Oceania is closest to East Asia
2. **Migration routes**: Austronesian expansion from East Asia to Pacific islands
3. **Limited PGS data**: Few Oceanian populations in PGS training datasets
4. **5-pop model constraint**: Must assign to one of five groups

**Alternative**: Could be mapped to `OTH` or kept separate in expanded models.

### Finnish → EUR (European)

**Rationale:**
1. **Geographic location**: Finland is in Europe
2. **Founder effects**: While Finns have unique population structure due to bottlenecks, they're still European
3. **PGS applications**: Finnish cohorts (e.g., FinnGen) are typically grouped with European studies

## Alternative Models

### 6-Population Model (with separate Middle East)
EUR, AFR, EAS, SAS, AMR, **MID**

**Pros:**
- Recognizes Middle Eastern populations as distinct
- More granular ancestry resolution
- Better for regional studies

**Cons:**
- Not standard in most PGS resources
- Requires PGS-Calc custom reference panel setup
- Limited PGS training data for Middle Eastern populations

### 7-Population Model (with Middle East + Oceania)
EUR, AFR, EAS, SAS, AMR, **MID**, **OCE**

**Pros:**
- Maximum granularity
- Recognizes all major geographic regions

**Cons:**
- Reduced sample sizes per group
- Compatibility issues with existing tools
- Very limited PGS data for some groups

## Recommendations

### For General Use (PGS-Calc, Standard Pipelines)
✅ **Use 5-population model**: EUR, AFR, EAS, SAS, AMR
- Best compatibility
- Largest training datasets
- Standard in literature

### For Research with Middle Eastern Focus
Consider 6-population model (EUR, AFR, EAS, SAS, AMR, MID) but:
- Document clearly in methods
- May need custom PGS normalization
- Check PGS-Calc compatibility

### For Your SNPster Pipeline
✅ **Current 5-population model is optimal** because:
1. Compatible with PGS-Calc ancestry adjustment
2. Matches your PGS panel libraries (blood, cancer, cardio, etc.)
3. Standard in clinical genomics
4. Maximizes sample sizes for reference populations

## Implementation in Code

The Python script includes both mappings - uncomment lines to switch:

```python
# Default: 5-population model (recommended)
POPULATION_INFERENCE_MAPPING = {
    'mid': 'EUR',  # Middle East → Europe
    # ... other mappings
}

# Alternative: 6-population model
# POPULATION_INFERENCE_MAPPING = {
#     'mid': 'MID',  # Keep Middle East separate
#     # ... other mappings
# }
```

## References

1. **PGS-Calc documentation**: https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html
2. **1000 Genomes populations**: https://www.internationalgenome.org/category/population/
3. **Population descriptors framework**: National Academies (2023) "Using Population Descriptors in Genetics and Genomics Research"
4. **Genetic structure**: Rosenberg et al. (2002) Science "Genetic Structure of Human Populations"
5. **PGS transferability**: Martin et al. (2019) Nature Genetics "Clinical use of current polygenic risk scores may exacerbate health disparities"

## Conclusion

**The 5-population model with Middle East → EUR is scientifically sound and optimal for your polygenic scoring pipeline.**

If you need to justify this choice:
- It's the standard model used by PGS-Calc and the broader field
- It's based on genetic clustering, not geography or labels
- It maximizes compatibility with existing PGS resources
- Middle Eastern populations are genetically part of West Eurasia
