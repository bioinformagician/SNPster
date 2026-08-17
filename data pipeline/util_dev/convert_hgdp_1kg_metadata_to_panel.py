#!/usr/bin/env python3
"""
Convert HGDP+1KG gnomAD metadata to standard panel format.

Converts the gnomAD HGDP+1KG metadata file to the standard panel format expected
by ancestry inference tools (sample, pop, super_pop, gender).

Maps genetic regions to standard 5 superpopulations:
- AFR (African), AMR (American), EAS (East Asian), EUR (European), SAS (South Asian)
"""

import sys
import pandas as pd


# Mapping from gnomAD population codes to standard superpopulations
# Note: This uses the standard 5-population model (EUR, AFR, EAS, SAS, AMR)
# as used by PGS-Calc and most polygenic scoring resources.
POPULATION_INFERENCE_MAPPING = {
    'afr': 'AFR',  # African
    'amr': 'AMR',  # American/Latino
    'eas': 'EAS',  # East Asian
    'sas': 'SAS',  # South Asian
    'nfe': 'EUR',  # Non-Finnish European
    'fin': 'EUR',  # Finnish (European)
    'mid': 'EUR',  # Middle Eastern (mapping to EUR for 5-pop model - genetically part of West Eurasia)
    'oth': 'OTH',  # Other (will be filtered or kept as-is)
}

# Alternative 6-population model (uncomment to use separate Middle Eastern category)
# POPULATION_INFERENCE_MAPPING = {
#     'afr': 'AFR', 'amr': 'AMR', 'eas': 'EAS', 'sas': 'SAS',
#     'nfe': 'EUR', 'fin': 'EUR', 'mid': 'MID', 'oth': 'OTH',
# }

# Mapping from HGDP/TGP genetic regions to standard superpopulations
GENETIC_REGION_MAPPING = {
    'AFR': 'AFR',  # Africa
    'EUR': 'EUR',  # Europe
    'EAS': 'EAS',  # East Asia
    'AMR': 'AMR',  # Americas
    'CSA': 'SAS',  # Central/South Asia
    'MID': 'EUR',  # Middle East (mapping to EUR for 5-pop model - genetically closer to Europeans)
    'OCE': 'EAS',  # Oceania (mapping to EAS for proximity, or could use OTH)
}

# Alternative 7-population model (uncomment to use separate Middle Eastern category)
# GENETIC_REGION_MAPPING = {
#     'AFR': 'AFR', 'EUR': 'EUR', 'EAS': 'EAS', 'AMR': 'AMR',
#     'CSA': 'SAS', 'MID': 'MID', 'OCE': 'OCE',
# }


def convert_metadata(input_file, output_file, use_inference=True, keep_other=False):
    """
    Convert gnomAD metadata to panel format.
    
    Parameters
    ----------
    input_file : str
        Path to gnomad_meta_updated.tsv
    output_file : str
        Path to output panel file
    use_inference : bool
        If True, use population_inference.pop (preferred); 
        if False, use hgdp_tgp_meta.Genetic.region
    keep_other : bool
        If True, keep samples with 'OTH' superpopulation;
        if False, exclude them
    """
    print(f"Reading metadata from {input_file}...")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    print(f"Total samples in metadata: {len(df)}")
    
    # Extract sample ID
    sample_col = 's'
    
    # Extract population (detailed population)
    pop_col = 'population' if 'population' in df.columns else 'hgdp_tgp_meta.Population'
    
    # Extract superpopulation based on preference
    if use_inference and 'population_inference.pop' in df.columns:
        print("Using population_inference.pop for superpopulation assignment...")
        super_pop_col = 'population_inference.pop'
        mapping = POPULATION_INFERENCE_MAPPING
    else:
        print("Using hgdp_tgp_meta.Genetic.region for superpopulation assignment...")
        super_pop_col = 'hgdp_tgp_meta.Genetic.region'
        mapping = GENETIC_REGION_MAPPING
    
    # Extract sex - try columns in order of preference
    sex_col = None
    for col in ['project_meta.sex', 'sex_imputation.sex_karyotype', 'sex', 'bergstrom.sex', 'project_meta.v2_sex']:
        if col in df.columns and df[col].notna().any():
            sex_col = col
            print(f"Using {col} for sex information")
            break
    
    if sex_col is None:
        print("Warning: No valid sex column found, setting all to 'unknown'")
        df['sex_temp'] = 'unknown'
        sex_col = 'sex_temp'
    
    # Create panel dataframe
    panel_df = pd.DataFrame({
        'sample': df[sample_col],
        'pop': df[pop_col],
        'super_pop_raw': df[super_pop_col],
        'gender': df[sex_col]
    })
    
    # Map to standard superpopulations
    panel_df['super_pop'] = panel_df['super_pop_raw'].map(mapping)
    
    # Handle unmapped values
    unmapped = panel_df['super_pop'].isna()
    if unmapped.any():
        print(f"\nWarning: {unmapped.sum()} samples with unmapped superpopulation:")
        print(panel_df[unmapped]['super_pop_raw'].value_counts())
        # Set unmapped to 'OTH'
        panel_df.loc[unmapped, 'super_pop'] = 'OTH'
    
    # Filter out 'OTH' if requested
    if not keep_other:
        before = len(panel_df)
        panel_df = panel_df[panel_df['super_pop'] != 'OTH']
        print(f"Filtered out {before - len(panel_df)} samples with OTH superpopulation")
    
    # Standardize gender to lowercase (convert to string first to handle NaN)
    panel_df['gender'] = panel_df['gender'].astype(str).str.lower()
    # Replace 'nan' string with actual missing indicator
    panel_df.loc[panel_df['gender'] == 'nan', 'gender'] = 'unknown'
    panel_df.loc[panel_df['gender'] == 'na', 'gender'] = 'unknown'
    # Standardize M/F to male/female
    panel_df['gender'] = panel_df['gender'].replace({'m': 'male', 'f': 'female'})
    # Handle karyotype values (XX = female, XY = male)
    panel_df['gender'] = panel_df['gender'].replace({
        'xx': 'female',
        'xy': 'male',
        'x': 'female',  # Turner syndrome - typically classified as female
        'xxy': 'male',  # Klinefelter syndrome - typically classified as male
        'xyy': 'male',  # XYY syndrome - typically classified as male
        'ambiguous': 'unknown'
    })
    
    # Select final columns
    panel_df = panel_df[['sample', 'pop', 'super_pop', 'gender']]
    
    # Summary statistics
    print(f"\nPanel summary:")
    print(f"  Total samples: {len(panel_df)}")
    print(f"\n  Superpopulation counts:")
    print(panel_df['super_pop'].value_counts().to_string())
    print(f"\n  Gender distribution:")
    print(panel_df['gender'].value_counts().to_string())
    
    # Check for any missing values
    missing = panel_df.isna().sum()
    if missing.any():
        print(f"\nWarning: Missing values detected:")
        print(missing[missing > 0].to_string())
    
    # Save to file
    print(f"\nWriting panel to {output_file}...")
    panel_df.to_csv(output_file, sep='\t', index=False)
    print("Done!")
    
    return panel_df


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python convert_hgdp_1kg_metadata_to_panel.py <input_tsv> <output_panel>")
        print("\nExample:")
        print("  python convert_hgdp_1kg_metadata_to_panel.py gnomad_meta_updated.tsv hgdp_1kg_panel.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_metadata(input_file, output_file, use_inference=True, keep_other=False)
