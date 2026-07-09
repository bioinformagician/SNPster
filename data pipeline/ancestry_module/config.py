import os

# Path to 1000 Genomes reference VCFs (your existing Beagle reference files)
REFERENCE_VCF_DIR = os.getenv(
    "REFERENCE_VCF_DIR",
    "/srv/dependencies/imputation_runner/imputer/beagle_references"
)

# Path to 1000 Genomes population panel file
# Download from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
POPULATION_PANEL_FILE = os.getenv(
    "POPULATION_PANEL_FILE",
    "/srv/dependencies/imputation_runner/1000G_panel.txt"
)

# Chromosome to use for ancestry inference (chr22 is fastest, smallest autosome)
DEFAULT_CHROMOSOME = int(os.getenv("ANCESTRY_CHROMOSOME", "22"))

# Number of ancestral populations (5 = EUR, AFR, EAS, SAS, AMR)
K_POPULATIONS = int(os.getenv("K_POPULATIONS", "5"))

ANCESTRY_METHOD = os.getenv("ANCESTRY_METHOD", "ADMIXTURE")  # Options: ADMIXTURE, PCA, or other methods
REFERENCE_PANEL = os.getenv("REFERENCE_PANEL", "1000G 30x")  # Options: 1000G, HGDP, or other panels
