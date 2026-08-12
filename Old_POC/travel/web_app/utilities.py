import pandas as pd
import numpy as np
from scipy.stats import norm
from pathlib import Path
import os

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent
# Define paths relative to the script directory  
PROJECT_ROOT = SCRIPT_DIR.parent

# Allow environment variables to override default paths (optional)
if 'GENOMICS_PROJECT_ROOT' in os.environ:
    PROJECT_ROOT = Path(os.environ['GENOMICS_PROJECT_ROOT'])

# Define commonly used directories
GENOMES_DIR = Path(os.environ.get('GENOMICS_GENOMES_DIR', PROJECT_ROOT / "genomes"))
BACKGROUND_DATA_DIR = PROJECT_ROOT / "background_data"
WWW_DIR = SCRIPT_DIR / "www"

custom_css = """
.main-title {
    color: #2c3e50;
    margin-bottom: 2rem;
    text-align: center;
    padding: 1rem;
    background: linear-gradient(to right, #f8f9fa, #e9ecef);
    border-radius: 8px;
}
body {
    margin: 0;
    padding: 0;
}

/* Make carousel full width and adjust spacing */
.carousel {
    width: 100vw; /* Full viewport width */
    margin-left: calc(-50vw + 50%); /* Negative margin trick for full width */
    margin-right: calc(-50vw + 50%);
    margin-top: -1rem; /* Remove top spacing */
    position: relative;
    left: 50%;
    right: 50%;
    transform: translateX(-50%);
}

.carousel-item {
    width: 100%;
    height: 600px; /* You can adjust this value */
    overflow: hidden;
}

.carousel-item img {
    width: 100%;
    height: 600px; /* Match the container height */
    object-fit: cover;
    object-position: center;
}

.carousel-caption {
    background: rgba(0, 0, 0, 0.6);
    padding: 20px;
    border-radius: 5px;
}

/* Container modifications */
.container-fluid {
    padding-left: 0 !important;
    padding-right: 0 !important;
}

.analysis-card {
    background-color: transparent;
    border: none; /* Explicitly remove borders */
    border-radius: 0px; /* Keeps the rounded corners, if needed */
    padding: 1.5rem;
    margin: 1rem 0;
    box-shadow: 0 0px 0px rgba(0, 0, 0, 0.1); /* Optional: Retain or remove for a clean look */
}

.consent-text {
    background-color: #f8f9fa;
    padding: 1rem;
    border-radius: 5px;
    margin: 1rem 0;
    font-size: 0.9rem;
    color: #495057;
}

.overview-card {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 2rem;
    margin: 1rem 0;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.card-title {
    color: #2c3e50;
    border-bottom: 2px solid #e9ecef;
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
}

.info-card {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 2rem;
    margin: 1rem 0;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.info-text {
    line-height: 1.6;
    color: #495057;
}

.feature-list {
    list-style-type: none;
    padding-left: 0;
    margin-top: 1rem;
}

.feature-list li {
    margin-bottom: 0.5rem;
    padding-left: 1.5rem;
    position: relative;
}

.feature-list li:before {
    content: "•";
    color: #2c3e50;
    position: absolute;
    left: 0;
}

.required-field {
    color: #dc3545;
    font-size: 0.8rem;
    margin-top: 0.25rem;
}

.centered-container {
    color: #2c3e50;
    border-bottom: 2px solid #e9ecef;
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
    max-width: 66%;
    margin: 0 auto;
}

.upload-instruction {
    margin-bottom: 0;
    line-height: 1.2;
}

.required-field {
    color: #dc3545;
    margin-top: 0.2rem;
}

:root {
    --primary-color: #2C3E50;
    --secondary-color: #3498DB;
    --accent-color: #E74C3C;
    --background-color: #ECF0F1;  /* Change this line */
    --text-color: #2C3E50;
    --card-background: #FFFFFF;
}

.scrollable-text {
    height: 300px;  /* Adjust this value to control the height of the scrollable area */
    overflow-y: auto;
    padding: 15px;
    border: 1px solid #e0e0e0;
    border-radius: 5px;
    background-color: rgba(255, 255, 255, 0.8);
}

.transparent-card {
        background-color: rgba(255, 255, 255, 0.2);  /* White with 20% opacity */
        backdrop-filter: blur(10px);  /* Optional: adds a blur effect */
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin-bottom: 2rem;
        border: none;
        transition: transform 0.2s;
    }

    /* Optional hover effect */
    .transparent-card:hover {
        background-color: rgba(255, 255, 255, 0.3);  /* Slightly more opaque on hover */
    }
.navbar {
        position: sticky;
        top: 0;
        z-index: 1000;  /* Ensures navbar stays on top of other content */
        background-color: var(--surface-color);  /* Use your navbar color */
    }
"""

# Use relative path for GWAS data - with environment variable override
default_gwas_path = BACKGROUND_DATA_DIR/ "gwas_data_filtered_phased_low_filt.csv"
gwas_data_path = Path(os.environ.get('GENOMICS_GWAS_DATA', default_gwas_path))


def prepare_data(gwas_data_path, personal_data):
    ###############################################################################
    # 1. Read in the GWAS data
    ###############################################################################
    gwas_data = pd.read_csv(gwas_data_path)


    ###############################################################################
    # 4. Merge on 'rsid'
    ###############################################################################
    merged_data = pd.merge(
        gwas_data,
        personal_data,
        on="rsid",
        how="inner"
    )
    
    #print nrow 
    # Print the number of rows after dropping duplicates
    print(f"Number of rows after dropping duplicates: {merged_data.shape[0]}")
    merged_data = merged_data.drop_duplicates() #the primary key for gwas_data is rsid + trait, but as personal_data doesnt have trait, there will be dublicates 
    # Drop duplicate rows where 'rsid' and 'trait' are the same
    merged_data = merged_data.drop_duplicates(subset=['rsid', 'trait'], keep='first')

    
    print(f"Number of rows after dropping duplicates: {merged_data.shape[0]}")

    merged_data.rename(
        columns={
            "genotype": "personal_genotype",
        },
        inplace=True,
    )




    ###############################################################################
    # 5. Group by trait and count
    ###############################################################################
    counts_for_each_trait = (
        merged_data
        .groupby("trait")
        .size()
        .reset_index(name="count")
    )

    ###############################################################################
    # 6. Row-wise risk_status determination
    ###############################################################################
    def compute_risk_status(row):
        """
        Returns:
          2 if the genotype is homozygous for the risk allele,
          1 if the genotype is heterozygous for the risk allele,
          0 otherwise
        """
        genotype = row["personal_genotype"]
        risk_allele = row["risk_allele"]
        # Count how many times the risk_allele appears in the genotype
        if genotype.count(risk_allele) == 2:
            return 2
        elif risk_allele in genotype:
            return 1
        else:
            return 0

    merged_data["risk_status"] = merged_data.apply(compute_risk_status, axis=1)



    ###############################################################################
    # 5. Get rare snp df
    ###############################################################################

    #pick 5 rarest snps by looking for min in the risk.allele.frequency
    merged_data["other_allele_frequency"] = 1 - merged_data["risk_allele_frequency"]

    # Filter rows where "-" is NOT in personal_genotype
    merged_data_no_minus = merged_data[~merged_data['personal_genotype'].str.contains('-', na=False)]

    #if risk_stats = 1, multiply risk_allele_frequency by other_allele_frequency, if = 0 other_allele_frequency^2, if 2 risk_allele_frequency^2

    merged_data_no_minus['genotype_probability'] = np.where(merged_data_no_minus["risk_status"] == 1, merged_data_no_minus["risk_allele_frequency"] * merged_data_no_minus["other_allele_frequency"], np.where(merged_data_no_minus["risk_status"] == 0, merged_data_no_minus["other_allele_frequency"]**2, merged_data_no_minus["risk_allele_frequency"]**2))
    
    #if '-' in the personal_genotype, set the genotype_probability to NA

    rarest_genotypes = merged_data_no_minus.nsmallest(5, 'genotype_probability')

    rare_snps_no_risk = merged_data_no_minus[merged_data_no_minus['risk_status']==0].nsmallest(5, 'other_allele_frequency')


    rare_snps_risk = merged_data_no_minus[merged_data_no_minus['risk_status']>0].nsmallest(5, 'risk_allele_frequency')

    
    rare_snps = pd.concat([rare_snps_risk, rare_snps_no_risk, rarest_genotypes])

    #remove first column


    #REMOVE confidence_interval, or, trait_direction, population_score, chromosome, position, risk_status
    rare_snps = rare_snps.drop(columns=["Unnamed: 0","confidence_interval", "or", "trait_direction", "population_score", "chromosome", "position", "risk_status", "pmid","odds_ratio", "trait"])


    ###############################################################################
    # 7. add status category
    ###############################################################################


    def status_category(x):
        # If risk_status <= 1 -> "good", else -> "bad"
        return "good" if x <= 1 else "bad"

    merged_data["risk_status_category"] = merged_data["risk_status"].apply(status_category)


    ###############################################################################
    # 8. Calculate final_risk and remove NaNs
    ###############################################################################

    #For the traits (where or = False, i.e something that has a unit increase or decrease), adjust unit direction
    merged_data.loc[(merged_data["or"] == False) & (merged_data["trait_direction"] == "decrease"), "odds_ratio"] *= -1

    merged_data["final_risk"] = (
        np.where(
            merged_data["or"],  # Apply log transformation only where `or == True`
            np.log(merged_data["odds_ratio"]),
            merged_data["odds_ratio"]  # Use raw odds_ratio where `or == False`
        ) * merged_data["risk_status"]
        - merged_data["population_score"]
    )

  

    merged_data_pruned = merged_data.dropna(subset=["final_risk"])

    

    ###############################################################################
    # 9. Aggregate by 'trait'
    ###############################################################################
    #make pmid string
    merged_data_pruned["pmid"] = merged_data_pruned["pmid"].astype(str)
    
    aggregated_true = (
        merged_data_pruned
        .groupby("trait", as_index=False)
        .agg({                       # Count unique rsid
            "final_risk":        "sum",                            # Sum final_risk                       # Count unique pmid
            "chromosome":        "nunique",                        # Count unique chromosome
            "or":                "first",                          # Keep the first or
            "population_score":  "sum",                            # Sum population_score
            "pmid":       lambda x: ",".join(map(str, set(x))),  # Concatenate pmid values
            "rsid":       lambda x: ",".join(map(str, set(x))),
            "Category":         "first"   # Concatenate rsid values
        })
    )

    #for the rows where or == False, set the final_risk to be ((population_risk+final_risk)/population_risk)*100

    #aggregated_true["final_risk"] = np.where(aggregated_true["or"] == False, np.log(((aggregated_true["population_score"] + aggregated_true["final_risk"]) / aggregated_true["population_score"])), aggregated_true["final_risk"])

    # Calculate percentile ~ pnorm(final_risk)
    aggregated_true["percentile"] = norm.cdf(aggregated_true["final_risk"]) * 100

    # Filter to only include rows where 'or' == True
    

    # Sort final_data by percentile
    final_data = aggregated_true
    final_data_or = final_data.sort_values("percentile")

    # Return any objects of interest
    return {"polygenomic_risk_score_data": final_data_or, "rare_snps": rare_snps, "full_data": merged_data_pruned}

def get_genome_file_path(filename):
    """
    Get the full path to a genome file in the genomes directory.
    
    Args:
        filename (str): The name of the genome file
        
    Returns:
        Path: Full path to the genome file
    """
    return GENOMES_DIR / filename

def list_available_genome_files():
    """
    List all available genome files in the genomes directory.
    
    Returns:
        dict: Dictionary mapping file names to string paths
    """
    if GENOMES_DIR.exists():
        paths = [f.name for f in GENOMES_DIR.glob("*.txt")]
        paths.append("")
        paths.sort()
        paths = {f: str(get_genome_file_path(f)) for f in paths if f}
        
        paths[''] = ''
        return paths
    
    return {}

def get_background_data_path(filename):
    """
    Get the full path to a background data file.
    
    Args:
        filename (str): The name of the background data file
        
    Returns:
        Path: Full path to the background data file
    """
    return BACKGROUND_DATA_DIR / filename

sports_df = pd.DataFrame({
    'sport': ['Running', 'Swimming', 'Power Lifting', 'Sprinting', 'Endurance'],
    'genetic_score': [85, 92, 78, 88, 95],
    'description': [
        'You have genetic variants associated with endurance running performance.',
        'Your genetics suggest good potential for swimming activities.',
        'Genetic markers indicate strength in power-based activities.',
        'You have favorable genetic variants for sprint performance.',
        'High genetic predisposition for endurance-based activities.'
    ]
})