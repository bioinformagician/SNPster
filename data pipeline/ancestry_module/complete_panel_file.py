#!/usr/bin/env python3
"""
Create a complete panel file with all samples from the reference VCF.
Fills in missing samples using gnomAD metadata.
"""

import pandas as pd
import subprocess
import sys

# Paths
CHR22_VCF = "/srv/dependencies/ancestry/hgdp_1kg_references/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"
EXISTING_PANEL = "/srv/dependencies/imputation_runner/hgdp_1kg_panel.txt"
GNOMAD_META = "/srv/dependencies/ancestry/hgdp_1kg_references/gnomad_meta_updated.tsv"
OUTPUT_PANEL = "/tmp/hgdp_1kg_panel_complete_10pop.txt"

# Population mapping from gnomAD to 10-pop system
POPULATION_MAP = {
    # 1000 Genomes populations
    'ACB': 'SSA',  # African Caribbean in Barbados
    'ASW': 'SSA',  # African Ancestry in Southwest US
    'BEB': 'SAS',  # Bengali in Bangladesh
    'CDX': 'SEA',  # Chinese Dai in Xishuangbanna
    'CEU': 'EUR',  # Utah residents with Northern and Western European ancestry
    'CHB': 'EAS',  # Han Chinese in Beijing
    'CHS': 'EAS',  # Han Chinese South
    'CLM': 'AMR',  # Colombian in Medellin
    'ESN': 'SSA',  # Esan in Nigeria
    'FIN': 'EUR',  # Finnish in Finland
    'GBR': 'EUR',  # British in England and Scotland
    'GIH': 'SAS',  # Gujarati Indian in Houston
    'GWD': 'SSA',  # Gambian in Western Division
    'IBS': 'EUR',  # Iberian populations in Spain
    'ITU': 'SAS',  # Indian Telugu in the UK
    'JPT': 'EAS',  # Japanese in Tokyo
    'KHV': 'SEA',  # Kinh in Ho Chi Minh City
    'LWK': 'SSA',  # Luhya in Webuye, Kenya
    'MSL': 'SSA',  # Mende in Sierra Leone
    'MXL': 'AMR',  # Mexican Ancestry in Los Angeles
    'PEL': 'AMR',  # Peruvian in Lima
    'PJL': 'SAS',  # Punjabi in Lahore
    'PUR': 'AMR',  # Puerto Rican in Puerto Rico
    'STU': 'SAS',  # Sri Lankan Tamil in the UK
    'TSI': 'EUR',  # Toscani in Italia
    'YRI': 'SSA',  # Yoruba in Ibadan, Nigeria
    
    # HGDP populations
    'Adygei': 'EUR',
    'Balochi': 'SAS',
    'BantuKenya': 'SSA',
    'BantuSouthAfrica': 'SSA',
    'Basque': 'EUR',
    'Bedouin': 'MID',
    'BergamoItalian': 'EUR',
    'Biaka': 'SSA',
    'Brahui': 'SAS',
    'Burusho': 'SAS',
    'Cambodian': 'SEA',
    'Colombian': 'AMR',
    'Dai': 'SEA',
    'Daur': 'EAS',
    'Druze': 'MID',
    'French': 'EUR',
    'Han': 'EAS',
    'Hazara': 'SAS',
    'Hezhen': 'EAS',
    'Japanese': 'EAS',
    'Kalash': 'SAS',
    'Karitiana': 'NAM',
    'Lahu': 'EAS',
    'Makrani': 'SAS',
    'Mandenka': 'SSA',
    'Maya': 'NAM',
    'Mbuti': 'SSA',
    'Melanesian': 'OCE',  # Oceanian!
    'Miao': 'EAS',
    'Mongolian': 'EAS',
    'Mozabite': 'MID',
    'Naxi': 'EAS',
    'NorthernHan': 'EAS',
    'Orcadian': 'EUR',
    'Oroqen': 'EAS',
    'Palestinian': 'MID',
    'Papuan': 'OCE',  # Oceanian!
    'Pathan': 'SAS',
    'Pima': 'NAM',
    'Russian': 'EUR',
    'San': 'SSA',
    'Sardinian': 'EUR',
    'She': 'EAS',
    'Sindhi': 'SAS',
    'Surui': 'NAM',
    'Tu': 'EAS',
    'Tujia': 'EAS',
    'Tuscan': 'EUR',
    'Uygur': 'CAS',
    'Xibo': 'EAS',
    'Yakut': 'EAS',
    'Yi': 'EAS',
    'Yoruba': 'SSA',
}

print("Reading existing panel file...")
panel_df = pd.read_csv(EXISTING_PANEL, sep='\t')
print(f"  Existing panel has {len(panel_df)} samples")

print("\nGetting sample list from chr22 VCF...")
result = subprocess.run(
    ['bcftools', 'query', '-l', CHR22_VCF],
    capture_output=True, text=True, check=True
)
vcf_samples = set(result.stdout.strip().split('\n'))
print(f"  chr22 VCF has {len(vcf_samples)} samples")

existing_samples = set(panel_df['sample'])
missing_samples = vcf_samples - existing_samples
print(f"  {len(missing_samples)} samples missing from panel")

if missing_samples:
    print("\nReading gnomAD metadata for missing samples...")
    # Read only necessary columns to save memory
    meta_df = pd.read_csv(
        GNOMAD_META,
        sep='\t',
        usecols=['s', 'hgdp_tgp_meta.Population', 'sex'],
        dtype=str
    )
    meta_df.columns = ['sample', 'pop', 'sex']
    
    # Filter to missing samples only
    missing_df = meta_df[meta_df['sample'].isin(missing_samples)].copy()
    
    if len(missing_df) == 0:
        print(f"  ERROR: None of the missing samples found in gnomAD metadata!")
        sys.exit(1)
    
    print(f"  Found {len(missing_df)} missing samples in gnomAD metadata")
    
    # Map population to super_pop
    missing_df['super_pop'] = missing_df['pop'].map(POPULATION_MAP)
    
    # Check for unmapped populations
    unmapped = missing_df[missing_df['super_pop'].isna()]
    if len(unmapped) > 0:
        print(f"\n  WARNING: {len(unmapped)} samples have unmapped populations:")
        print(unmapped[['sample', 'pop']].drop_duplicates('pop'))
        # Default to EUR for unmapped (conservative choice)
        missing_df['super_pop'] = missing_df['super_pop'].fillna('EUR')
    
    # Map sex (M/F to male/female)
    missing_df['gender'] = missing_df['sex'].map({'M': 'male', 'F': 'female', 'NA': 'unknown'})
    missing_df['gender'] = missing_df['gender'].fillna('unknown')
    
    # Select columns to match panel format
    missing_df = missing_df[['sample', 'pop', 'super_pop', 'gender']]
    
    # Show population distribution of missing samples
    print("\n  Population distribution of missing samples:")
    for pop, count in missing_df['super_pop'].value_counts().items():
        print(f"    {pop}: {count}")
    
    # Combine with existing panel
    complete_panel = pd.concat([panel_df, missing_df], ignore_index=True)
else:
    complete_panel = panel_df

print(f"\nComplete panel has {len(complete_panel)} samples")
print(f"  Writing to {OUTPUT_PANEL}")
complete_panel.to_csv(OUTPUT_PANEL, sep='\t', index=False)

# Show final population distribution
print("\nFinal population distribution:")
for pop, count in complete_panel['super_pop'].value_counts().sort_index().items():
    print(f"  {pop}: {count}")

print("\n✓ Complete panel file created!")
print(f"\nTo activate:")
print(f"  sudo cp {OUTPUT_PANEL} /srv/dependencies/ancestry/hgdp_1kg_references/hgdp_1kg_panel_10pop.txt")
print(f"  sudo cp {OUTPUT_PANEL} /srv/dependencies/imputation_runner/hgdp_1kg_panel.txt")
