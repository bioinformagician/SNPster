import os

# Path to 1000 Genomes reference VCFs (your existing Beagle reference files)
REFERENCE_VCF_DIR = os.getenv(
    "REFERENCE_VCF_DIR",
    "/srv/dependencies/imputation_runner/imputer/beagle_references"
)

REFERENCE_VCF_PATTERN = os.getenv(
    "REFERENCE_VCF_PATTERN",
    "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.numericCHR.SNV_INDEL_biallelic.vcf.gz"
)


POPULATION_PANEL_FILE = os.getenv(
    "POPULATION_PANEL_FILE",
    "/srv/dependencies/imputation_runner/hgdp_1kg_panel.txt"
)

# Chromosome to use for ancestry inference (chr22 is fastest, smallest autosome)
DEFAULT_CHROMOSOME = int(os.getenv("ANCESTRY_CHROMOSOME", "22"))

# Number of ancestral populations - flexible based on reference panel
# 5 = classic super populations (EUR, AFR, EAS, SAS, AMR)
# 10 = enhanced stratification (adds MID, SSA, SEA, CAS, NAM, OCE)
K_POPULATIONS = int(os.getenv("K_POPULATIONS", "5"))

ANCESTRY_METHOD = os.getenv("ANCESTRY_METHOD", "ADMIXTURE")  # Options: ADMIXTURE, PCA, or other methods
REFERENCE_PANEL = os.getenv("REFERENCE_PANEL", "1000G 30x HGDP")  # Options: 1000G, HGDP, or other panels

# Full population labels (10 categories)
ANCESTRY_LABELS = {
    'EUR': 'European',
    'EAS': 'East Asian',
    'SAS': 'South Asian',
    'SSA': 'Sub-Saharan African',
    'MID': 'Greater Middle Eastern',
    'AMR': 'Hispanic or Latin American',
    'SEA': 'South East Asian',
    'CAS': 'Central Asian',
    'NAM': 'Native American',
    'OCE': 'Oceanian',
    # Legacy support
    'AFR': 'African (legacy)'
}

# Backward compatibility - map old codes to new codes
LEGACY_POPULATION_MAP = {
    'AFR': 'SSA'  # Old generic African → Sub-Saharan African
}

# Expected database columns (10 populations)
DB_ANCESTRY_COLUMNS = ['eur', 'eas', 'sas', 'ssa', 'mid', 'amr', 'sea', 'cas', 'nam', 'oce']

# Population code validation
VALID_POPULATION_CODES = set(ANCESTRY_LABELS.keys())

