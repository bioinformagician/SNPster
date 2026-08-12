from shiny import App, render, ui, reactive
from shiny.ui import navbar_options
import pandas as pd
import plotly.express as px
import numpy as np
from htmltools import css
from scipy.stats import norm
from shinywidgets import output_widget, render_widget  
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.stats import norm
from pathlib import Path
import os
#from utilities import *

"""Quick fix so posit can use the functions from utilities.py"""

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
    merged_data_no_minus = merged_data[~merged_data['personal_genotype'].str.contains('-', na=False)].copy()

    #if risk_stats = 1, multiply risk_allele_frequency by other_allele_frequency, if = 0 other_allele_frequency^2, if 2 risk_allele_frequency^2

    merged_data_no_minus.loc[:, 'genotype_probability'] = np.where(
        merged_data_no_minus["risk_status"] == 1, 
        merged_data_no_minus["risk_allele_frequency"] * merged_data_no_minus["other_allele_frequency"], 
        np.where(
            merged_data_no_minus["risk_status"] == 0, 
            merged_data_no_minus["other_allele_frequency"]**2, 
            merged_data_no_minus["risk_allele_frequency"]**2
        )
    )
    
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

    # Handle potential log of invalid values (0 or negative numbers)
    with np.errstate(invalid='ignore', divide='ignore'):
        merged_data["final_risk"] = (
            np.where(
                merged_data["or"],  # Apply log transformation only where `or == True`
                np.where(
                    merged_data["odds_ratio"] > 0,  # Only log positive values
                    np.log(merged_data["odds_ratio"]),
                    0  # Use 0 for invalid values
                ),
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

""" Doc string finished"""

# from IPython.display import HTML  # Not needed - using ui.HTML instead
#from PIL import Image
import base64
import io
import os
from pathlib import Path


def load_image_as_base64(image_path):
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode("utf-8")

image4_base64 = load_image_as_base64(WWW_DIR / "blood5.jpg")
image3_base64 = load_image_as_base64(WWW_DIR / "sports2.jpg")
image1_base64 = load_image_as_base64(WWW_DIR / "gene2.webp")
image2_base64 = load_image_as_base64(WWW_DIR / "pexels-pixabay-531880.jpg")

def get_image_with_validation():
    """Return a validated background image, with fallback options"""
    # Try to return one of the loaded images, with fallbacks
    if image1_base64:
        return image1_base64
    elif image2_base64:
        return image2_base64
    elif image3_base64:
        return image3_base64
    elif image4_base64:
        return image4_base64
    else:
        return None


GENOME_PATHS = list_available_genome_files()  # This will be a list of paths to genome files
#imported objects: gwas_data_path, prepare_data function and custom_css

app_ui = ui.page_navbar(
    
    ui.nav_panel(
        "Genomic Analysis",
        # Method 1: Using a web URL directly
        ui.tags.style("""
    /* root variables */
    :root {
        --primary-color: #2C3E50;
        --secondary-color: #3498DB;
        --accent-color: #E74C3C;
        --background-color: whitesmoke;
        --text-color: #2C3E50;
        --card-background: #FFFFFF;

        /* If you want a surface-color variable, define it here: */
        --surface-color: #b1b1b3;
    }

    body {
        background-color: #F4FCFF; #change background here
        color: var(--text-color);
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }

    /* Navbar styling */
    .navbar {
        position: sticky;
        top: 0;
        z-index: 1000;
        background-color: black !important;  /* Overwrite default with !important */
        padding: 1rem 1.5rem !important;
        min-height: 70px;
    }

    /* Adjust other navbar elements if needed */
    .navbar-brand {
        font-size: 1.5rem !important;
        padding: 0.5rem 1rem !important;
    }

    .nav-link {
        font-size: 1.2rem !important;
        padding: 0.7rem 1rem !important;
    }

    .nav-link.active {
        font-weight: bold;
        border-bottom: 3px solid var(--primary-color);
    }

    .navbar-default {
        background-color: black !important;
    }
"""),
        ui.page_fluid(
            ui.tags.style(custom_css),
            ui.tags.head(
                                ui.tags.script(src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"),
                                ui.tags.link(
        rel="stylesheet",
        href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"
    ),
                                
                            ),
            ui.tags.script("""
    document.addEventListener('click', function(e) {
        if (e.target.classList.contains('card-header') || e.target.parentElement.classList.contains('card-header')) {
            var content = e.target.closest('.card').querySelector('.collapse');
            bootstrap.Collapse.getOrCreateInstance(content).toggle();
        }
    });"""),
    ui.tags.script(
    """
    document.addEventListener('DOMContentLoaded', function () {
    // Introduction section
    const collapseIntro = document.getElementById('collapseIntro');
    const toggleIconIntro = document.getElementById('toggle-icon-1'); // matches above

    if (collapseIntro) {
        collapseIntro.addEventListener('show.bs.collapse', function () {
            toggleIconIntro.classList.remove('fa-plus');
            toggleIconIntro.classList.add('fa-minus');
        });
        collapseIntro.addEventListener('hide.bs.collapse', function () {
            toggleIconIntro.classList.remove('fa-minus');
            toggleIconIntro.classList.add('fa-plus');
        });
    }

    // Consent section
    const collapseConsent = document.getElementById('collapseConsent');
    const toggleIconConsent = document.getElementById('toggle-icon-2'); // matches above

    if (collapseConsent) {
        collapseConsent.addEventListener('show.bs.collapse', function () {
            toggleIconConsent.classList.remove('fa-plus');
            toggleIconConsent.classList.add('fa-minus');
        });
        collapseConsent.addEventListener('hide.bs.collapse', function () {
            toggleIconConsent.classList.remove('fa-minus');
            toggleIconConsent.classList.add('fa-plus');
        });
    }

    // etc...
});
    """
),
            #ui.h1("SNPSTER - Genomic Analysis Platform", class_="main-title"),
ui.tags.div(
    {"id": "imageCarousel", "class": "carousel slide", "data-bs-ride": "carousel"},
    # Indicators
    ui.tags.div(
        {"class": "carousel-indicators"},
        ui.tags.button(
            {"type": "button", "data-bs-target": "#imageCarousel", 
             "data-bs-slide-to": "0", "class": "active", "aria-current": "true"},
        ),
        ui.tags.button(
            {"type": "button", "data-bs-target": "#imageCarousel", 
             "data-bs-slide-to": "1"},
        ),
        ui.tags.button(
            {"type": "button", "data-bs-target": "#imageCarousel", 
             "data-bs-slide-to": "2"},
        ),
        ui.tags.button(
            {"type": "button", "data-bs-target": "#imageCarousel", 
             "data-bs-slide-to": "3"},
        ),
    ),
    # Slides
    ui.tags.div(
        {"class": "carousel-inner"},
        # First slide
        ui.tags.div(
            {"class": "carousel-item active"},
            ui.tags.img(
                {
                    "src": f"data:image/jpeg;base64,{image1_base64}",
                    "class": "d-block w-100",
                    "alt": "First slide"
                }
            ),
            ui.tags.div(
                {"class": "carousel-caption"},
                ui.h5("Analyze Your Data"),
                ui.p("Identify your most unique genotypes and biggest predispositions")
            )
        ),
        # Second slide
        ui.tags.div(
            {"class": "carousel-item"},
            ui.tags.img(
                {
                    "src": f"data:image/jpeg;base64,{image3_base64}",
                    "class": "d-block w-100",
                    "alt": "Second slide"
                }
            ),
            ui.tags.div(
                {"class": "carousel-caption"},
                ui.h5("View your popular traits"),
                ui.p("See if you are at risk for red haired children or have the power gene")
            )
        ),
        # Third slide
        ui.tags.div(
            {"class": "carousel-item"},
            ui.tags.img(
                {
                    "src": f"data:image/jpeg;base64,{image2_base64}",
                    "class": "d-block w-100",
                    "alt": "Third slide"
                }
            ),
            ui.tags.div(
                {"class": "carousel-caption"},
                ui.h5("Get Insights"),
                ui.p("Understand your genetic profile")
            )
        ),
        # Fourth slide (Blood Panel)
        ui.tags.div(
            {"class": "carousel-item"},
            ui.tags.img(
                {
                    "src": f"data:image/jpeg;base64,{image4_base64}",
                    "class": "d-block w-100",
                    "alt": "Fourth slide"
                }
            ),
            ui.tags.div(
                {"class": "carousel-caption"},
                ui.h5("Blood Panel"),
                
                ui.p("Explore detailed insights into your blood biomarkers."),
                # Navigation Button
                
            )
        )
    ),
    # Controls
    ui.tags.button(
        {"class": "carousel-control-prev", "type": "button", 
         "data-bs-target": "#imageCarousel", "data-bs-slide": "prev"},
        ui.tags.span({"class": "carousel-control-prev-icon", "aria-hidden": "true"}),
        ui.tags.span({"class": "visually-hidden"}, "Previous")
    ),
    ui.tags.button(
        {"class": "carousel-control-next", "type": "button", 
         "data-bs-target": "#imageCarousel", "data-bs-slide": "next"},
        ui.tags.span({"class": "carousel-control-next-icon", "aria-hidden": "true"}),
        ui.tags.span({"class": "visually-hidden"}, "Next")
    ),
),
       
    ui.hr(),
    ui.div(
    {"style": "margin-top: 6rem;"},  # Adds space above the next div if needed
    # Additional content for the next div
),

    ui.div(
    {"style": "text-align: center; margin: 2rem 0 4rem 0;"},  # Adjusted margin-bottom to 4rem
    ui.tags.div(
        {"style": "display: inline-block; position: relative;"},
        # Line above the text
        ui.tags.div(
            {"style": """
                position: absolute; 
                top: -12px; 
                left: 50%; 
                transform: translateX(-50%); 
                width: 200px; 
                height: 2px; 
                background-color: #6abf91;  /* Match the green color */
            """}
        ),
        # Title text
        ui.tags.h3(
            "SNPster Genetic Analysis Journey",
            style="""
                display: inline-block; 
                font-family: 'Arial', sans-serif; 
                font-weight: normal; 
                font-size: 50px; 
                color: #333;
            """
        ),
        # Line below the text
        ui.tags.div(
            {"style": """
                position: absolute; 
                bottom: -12px; 
                left: 50%; 
                transform: translateX(-50%); 
                width: 200px; 
                height: 2px; 
                background-color: #6abf91;  /* Match the green color */
            """}
        ),
    )
),
ui.div(
    {"style": "margin-top: 6rem;"},  # Adds space above the next div if needed
    # Additional content for the next div
),

    ui.div(
    {
        "style": "display: flex; justify-content: center; gap: 4rem; margin: 2rem 0;"
    },  # Flexbox container for horizontal alignment and spacing
    ui.div(
            {"style": "text-align: center; margin-bottom: 3rem;"},  # Added margin-bottom for spacing
            ui.tags.i({"class": "fas fa-dna", "style": "font-size: 90px; color: #3498DB;"}),  # Larger DNA icon
            ui.tags.h4("Advance Science", style="margin-top: 1rem; color: #2C3E50;"),  # Title beneath the icon
            ui.tags.div(
                {"style": "font-size: 16px; color: #555; line-height: 1.5; text-align: center; margin: 0;"},
                ui.tags.ul(
                    {"style": "list-style: none; padding: 0; margin: 0;"},
                    ui.tags.li("✔ Support this websites development and accuracy."),
                    ui.tags.li("✔ Help advance research on population genetics"),
                    ui.tags.li("✔ Maintain complete privacy and anonymity"),
                )
            )
        ),
        ui.div(
            {"style": "text-align: center; margin-bottom: 3rem;"},  # Added margin-bottom for spacing
            ui.tags.i({"class": "fas fa-heartbeat", "style": "font-size: 90px; color: #E74C3C;"}),  # Larger Heartbeat icon
            ui.tags.h4("Health Risks", style="margin-top: 1rem; color: #2C3E50;"),  # Title beneath the icon
            ui.tags.ul(
                {"style": "list-style: none; padding: 0; margin: 0;"},
                ui.tags.li("✔ Learn about your health predisposition"),
                ui.tags.li("✔ Find litterature on your genetic risks"),
                ui.tags.li("? Compare yourself to friends and family"),
            ),
        ),
        ui.div(
            {"style": "text-align: center; margin-bottom: 3rem;"},  # Added margin-bottom for spacing
            ui.tags.i({"class": "fas fa-running", "style": "font-size: 90px; color: #2ECC71;"}),  # Larger Running icon
            ui.tags.h4("Traits", style="margin-top: 1rem; color: #2C3E50;"),  # Title beneath the icon
            ui.tags.ul(
                {"style": "list-style: none; padding: 0; margin: 0;"},
                ui.tags.li("✔ Get an overview of your athletic SNPs"),
                ui.tags.li("✔ Check your genome for popular traits"),
                ui.tags.li("✔ See your rarest SNPs"),
            ),
        ),
    ),
    # Inserted collapsible card
ui.div(
    {"class": "container", "style": "max-width: 1600px; margin: 2rem auto;"},  # Center the card and limit width
    ui.div(
        {"class": "card mb-4", "style": "border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"},  # Card styling
        # Header Section
        ui.div(
            {
                "class": "card-header d-flex justify-content-between align-items-center",
                "style": "cursor: pointer; background-color: #fff; border-bottom: 1px solid #ddd;",  # White header
                "data-bs-toggle": "collapse",
                "data-bs-target": "#collapseIntro",
                "aria-expanded": "false",
                "aria-controls": "collapseIntro"
            },
            ui.div(
                ui.h4("1. Introduction", class_="mb-0", style="color: #333;"),
                ui.p("Click to expand for detailed information", class_="text-muted mb-0"),
            ),
            ui.tags.i({"class": "fas fa-plus toggle-icon",
           "id": "toggle-icon-1",
           "style": "color: #3498DB;"})
        ),
        # Collapsible Content
        ui.div(
            {
                "id": "collapseIntro",
                "class": "collapse",
                "aria-labelledby": "headerTraits"
            },
            ui.div(
                {"class": "card-body", "style": "background-color: #fff; padding: 1.5rem; border-top: 1px solid #ddd;"},
                ui.p(
                                """Our platform provides genetic analysis using polygenic risk score(PGS) to help you understand your health predispositions 
                                and genetic traits. Using advanced genomic analysis techniques, we analyze your unique genetic data to provide 
                                insights about your:"""
                            ),
                            ui.tags.ul(
                                {"class": "feature-list"},
                                ui.tags.li("Disease risk assessments and health predispositions"),
                                ui.tags.li("Inherited traits and characteristics"),
                                ui.tags.li("Athletic performance potential"),
                                ui.tags.li("Uniqueness compared to other users"),
                            ),
                            ui.p(
                                """Simply upload your genetic data file, and our platform will analyze various genetic markers to 
                                provide you with detailed insights about your genetic profile. All analysis is performed with strict 
                                privacy measures in place, meaning that your data cannot be traced back to you. Your genetic data
                                holds a lot of information about you, but in reality the amount of your genome contained in your 
                                file from 23andme, ancestry, etc. Is actually just around 0.0002% of your genome - meaning information about 
                                600.000 base pairs out of the 3.100.000.000 it takes to create exactly you."""
                            ),
                            ui.p(
                                """SNPSTER's PGS approach is to an extend inspired by the previous open source service impute.me which
                                you can learn more about in detail following this PMID 32714365 on pubmed. The formula used for the PGS
                                calculation is given as such:"""
                            ),
                            ui.p([
                                "Population SNP Score: ",
                                ui.tags.span("$$\\text{Population score}_{snp} = \\text{frequency}_{snp} \\times 2 \\times \\text{beta}_{snp}$$"),
                                "Zero-centered Score: ",
                                ui.tags.span("$$\\text{Zero-centered score} = \\sum{\\text{Beta}_{snp} \\times \\text{Effect-allele-count}_{snp} - \\text{Population score}_{snp}}$$"),
                                "Z-score: ",
                                ui.tags.span("$$\\text{Z-score} = \\frac{\\text{Zero-centered score}}{\\text{Standard deviation}_{population}}$$")
                            ]),
                            ui.p(
                                """The data used for the PGS calculation is mostly pulled from the gwas catalog at https://www.ebi.ac.uk/gwas/ 
                                the 20th of november 2024 with some minor contribution from data on SNPEDIA. The 1000 genomes project phase 3 
                                30x dataset (GRCh38) is also being used to impute missing genotypes and phasing for correcting reported alleles"""
                            )
            )
        )
    ),
),

ui.div(
    {"class": "container", "style": "max-width: 1600px; margin: 2rem auto;"},  # Center the card and limit width
    ui.div(
        {"class": "card mb-4", "style": "border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"},  # Card styling
        # Header Section
        ui.div(
            {
                "class": "card-header d-flex justify-content-between align-items-center",
                "style": "cursor: pointer; background-color: #fff; border-bottom: 1px solid #ddd;",  # White header
                "data-bs-toggle": "collapse",
                "data-bs-target": "#collapseConsent",
                "aria-expanded": "true",  # Start open
                "aria-controls": "collapseConsent"
            },
            ui.div(
                ui.h4("2. Consent", class_="mb-0", style="color: #333;"),
                ui.p("Click to expand for detailed information", class_="text-muted mb-0"),
            ),
            ui.tags.i({"class": "fas fa-minus toggle-icon",
           "id": "toggle-icon-2",
           "style": "color: #3498DB;"})  # Start with minus icon
        ),
        # Collapsible Content
        ui.div(
            {
                "id": "collapseConsent",
                "class": "collapse show",  # Add "show" to make it open by default
                "aria-labelledby": "headerTraits"
            },
            ui.div(
                {"class": "card-body", "style": "background-color: #fff; padding: 1.5rem; border-top: 1px solid #ddd;"},
                ui.div(
                    {"class": "consent-text"},
                    ui.h4("Privacy and Consent"),
                    ui.p(
                        "By uploading your data, you consent to our analysis of your genetic information. ",
                        "Your data will be processed securely and confidentially with no personally identifiable information attached to your information.",
                        " For the purpose of further developing this platform's accuracy and features, your raw file input will be saved in a database with",
                        " no attachment to your identity. This also means that it will be impossible to remove your data on request, as it cannot be traced",
                        " back to you. For any questions, feedback, or requests, please contact snpster@gmail.com."
                    ),
                    ui.div({"class": "required-field"}, "* Required for analysis"),
                    ui.input_checkbox("consent", "I agree to the privacy, use, and consent terms", width='100%'),
                ),
            )
        )
    ),
),


                ui.card(
                    {
        "class": "centered-container",
        "style": "max-width: 1000px; margin: 2rem auto; padding: 1.5rem; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); border-radius: 8px; background-color: #fff;"
    },
                    
                    
                    ui.h3("3. Upload Your Genetic Data", class_="card-title"),
                    
                    ui.input_file(
                        "user_file",
                        "",
                        button_label="Browse Files",
                        placeholder="Upload Your File Here",
                        multiple=False,
                        accept=".txt,.csv,.tsv",
                        width="100%"
                    ),
                    
                    ui.div(
                    ui.p("Or choose a template genome", class_="upload-instruction"),
                    ),
                    ui.input_select(
                        "sample_genome",
                        label="Select From Existing Sample Files",
                        choices=list(GENOME_PATHS.keys()),
                        selected=""),
                    
                    ui.div(
                            ui.div({"class": "required-field"}, "* Required for analysis")
                        ),
                    
                        ui.h4("Select your superpopulation", class_="card-title"),
                        ui.layout_column_wrap(
                        
                        ui.input_radio_buttons(
                            "user_superpopulation",
                            "",  
                            choices=["European", "East Asian", "African", "South Asian", "Admixed American"], 
                            selected="European",
                            width="100%",
                            inline=False
                        ),
                        output_widget("world_map")

                ),

                    ui.input_task_button("run_analysis", "Analyze My Genome!", class_="btn btn-primary", width="100%",style="margin-top: 0.5rem;"),
                ),

                ui.panel_conditional( "input.consent",

                            ui.hr(),
                            ui.div(
                                    {"style": "background-color: #F4FCFF; padding: 1.5rem;"},
                                    #center title
                                    ui.div(
                        {"style": "text-align: center; margin: 2rem 0 4rem 0;"},  # Adjusted margin-bottom to 4rem
                        ui.tags.div(
                            {"style": "display: inline-block; position: relative;"},
                            # Line above the text
                            ui.tags.div(
                                {"style": """
                                    position: absolute; 
                                    top: -12px; 
                                    left: 50%; 
                                    transform: translateX(-50%); 
                                    width: 200px; 
                                    height: 2px; 
                                    background-color: #6abf91;  /* Match the green color */
                                """}
                            ),
                            # Title text
                            ui.tags.h3(
                                "Your genetic analysis results",
                                style="""
                                    display: inline-block; 
                                    font-family: 'Arial', sans-serif; 
                                    font-weight: normal; 
                                    font-size: 30px; 
                                    color: #333;
                                """
                            ),
                            # Line below the text
                            ui.tags.div(
                                {"style": """
                                    position: absolute; 
                                    bottom: -12px; 
                                    left: 50%; 
                                    transform: translateX(-50%); 
                                    width: 200px; 
                                    height: 2px; 
                                    background-color: #6abf91;  /* Match the green color */
                                """}
                            ),
                        )
                    ),

                        ui.layout_column_wrap(
                            ui.card(
                                {"class": "analysis-card"},
                                ui.h3("Overall disease category predisposition", class_="card-title",style="text-align: center;"),
                                output_widget("disease_category_radar_plot")
                            ),
                            ui.card(
                                {"class": "analysis-card"},
                                ui.h3("Your SNP Sunburst Overview", class_="card-title", style="text-align: center;"),
                                output_widget("pie_plot")
                            )
                        ),
                        ui.card(
                {"style": "height: 1000px;", "class": "analysis-card"},
                ui.div(
                    ui.h3("Disease Risk Overview", class_="card-title",style="text-align: center;"),
                    ui.div(
                        {"style": "display: flex; align-items: center; justify-content: center; gap: 40px;"},
                        ui.div(
                            {"style": "display: flex; align-items: center; color: #28a745;"},
                            ui.tags.i({"class": "fas fa-arrow-left", "style": "font-size: 24px; margin-right: 10px;"}),
                            ui.tags.span("Decreased risk", style="font-weight: bold;")
                        ),
                        ui.div(
                            {"style": "display: flex; align-items: center; color: #dc3545;"},
                            ui.tags.span("Increased risk", style="font-weight: bold; margin-right: 10px;"),
                            ui.tags.i({"class": "fas fa-arrow-right", "style": "font-size: 24px;"})
                        )
                    )
                ),
                output_widget('full_disease_overview')
            ),
                                        
                        ui.card(
                            {"class": "analysis-card"},
                            ui.h3("Detailed Disease Risk Analysis & litterature", class_="card-title",style="text-align: center;"),
                            ui.input_selectize(
                                "traits", 
                                label="Select condition for detailed analysis", 
                                multiple=False, 
                                choices=["Melanoma", "Type 2 Diabetes", "Heart Disease", "Alzheimer's"], 
                                selected="Melanoma"
                            ),
                            output_widget("disease_risk_plot"),
                            ui.hr(),  # Add a horizontal line to separate plot and description
                            ui.card(
                                {"class": "analysis-card", "style": "margin-top: 1rem;"},
                                ui.h4("Condition Description & Litterature", class_="card-title"),
                                ui.output_text("trait_description"),  # This will be populated from the server
                                
                                ui.output_ui("top_results_plot_2")),
                            )
                        )
                    )
    )),
    ui.nav_panel(
    "Sports Performance",
    ui.page_fluid(
        # TOP IMAGE SECTION
         # TOP IMAGE SECTION
        ui.div(
            {
                "style": """
                    text-align: center; 
                    margin: 0;  /* Remove default margins if desired */
                    padding: 0; /* Remove default padding if desired */
                """
            },
            ui.tags.img(
    {
        "src": f"data:image/jpeg;base64,{image3_base64}",  # Base64-encoded image
        "class": "d-block w-100",  # Bootstrap classes for full width and block display
        "alt": "Cropped Image",
        "style": """
            width: 100%;
            height: 70vh;        /* Full browser viewport height */
            object-fit: cover;    /* Ensures aspect ratio and crops excess */
            object-position: center; /* Centers the cropped portion */
            display: block;
            overflow: hidden;
        """
    }
),

        ),


        ui.div(
     {
        "style": """
            display: flex;
            justify-content: space-between; 
            align-items: center; 
            background-color: #E8F5E9; /* Light green tint */
            padding: 0px 0px; 
            position: relative;
            padding-top: 40px; /* Padding between image and this div */
        """
    },
    # Left Icon (Dumbbell with "Power")
    ui.div(
        {
            "style": "text-align: center; margin-left: 200px;"  # Adjust margin to move icon further in
        },
        ui.tags.i(
            {
                "class": "fas fa-dumbbell",  # Font Awesome Dumbbell icon
                "style": """
                    font-size: 80px; /* Adjust size of the icon */
                    color: #6abf91; 
                """
            }
        ),
        ui.tags.p(
            "Power",
            {
                "style": """
                    margin-top: 10px; /* Space between icon and text */
                    font-size: 18px; /* Text size */
                    color: #333; /* Text color */
                    font-weight: 500; /* Slightly bold */
                """
            }
        )
    ),
    # Title in the Center
    ui.div(
        
        {
            "style": "text-align: center;"
        },
        # Line above the text
        ui.tags.div(
            {
                "style": """
                    margin: 0 auto;
                    width: 220px; 
                    height: 2px; 
                    background-color: #6abf91;  
                """
            }
        ),
        # Title text
        ui.tags.h2(
            "Sports Traits",
            style="""
                font-family: 'Arial', sans-serif; 
                font-weight: 600; 
                font-size: 36px; 
                color: #333;
                margin: 10px 0; /* Reduced margin to bring title closer to lines */
            """
        ),
        # Line below the text
        ui.tags.div(
            {
                "style": """
                    margin: 0 auto;
                    width: 220px; 
                    height: 2px; 
                    background-color: #6abf91;  
                """
            }
        )
    ),
    # Right Icon (Running Person with "Endurance")
    ui.div(
        {
            "style": "text-align: center; margin-right: 200px;"  # Adjust margin to move icon further in
        },
        ui.tags.i(
            {
                "class": "fas fa-running",  # Font Awesome Running Person icon
                "style": """
                    font-size: 80px; /* Adjust size of the icon */
                    color: #6abf91; 
                """
            }
        ),
        ui.tags.p(
            "Endurance",
            {
                "style": """
                    margin-top: 10px; /* Space between icon and text */
                    font-size: 18px; /* Text size */
                    color: #333; /* Text color */
                    font-weight: 500; /* Slightly bold */
                """
            }
        )
    )
),

        ui.div(
            {
                "style": """
                    background-color: #E8F5E9; /* Light green tint */
                    text-align: center;
                    padding: 10px 10px; 
                    position: relative;
                """
            },
        

            #add space
            ui.hr(),
        
        ui.tags.style(custom_css),
        ui.output_ui("sports_cards"),

        ui.hr(),

        # You can keep the rest of your UI elements below as desired
        
    ),
    ),
),
    ui.nav_panel(
        "Blood Panel",
        {"id": "BloodPanel"},  # Add an explicit ID
        ui.page_fluid(

            # TOP IMAGE SECTION (same as sports)
            ui.div(
            
                {
                    "style": """
                        text-align: center; 
                        margin: 0;  
                        padding: 0; 
                    """
                },
                ui.tags.img(
                    {
                        "src": f"data:image/jpeg;base64,{image4_base64}",
                        "class": "d-block w-100",
                        "alt": "Blood Panel Image",
                        "style": """
                            width: 100%;
                            height: 70vh;
                            object-fit: cover;
                            object-position: center;
                            display: block;
                            overflow: hidden;
                        """
                    }
                )
            ),

            # Mild red background section
            ui.div(
                
                {
                    "style": """
                        background-color: #FFF5F5; /* Light red tint */
                        text-align: center;
                        padding: 40px 20px;
                        position: relative;
                        margin-bottom: 40px;
                    """
                },
                ui.div(
    {
        "style": """
            position: relative; /* Ensures child elements are positioned relative to this container */
            text-align: center; /* Center-aligns the title */
            padding: 20px 0; /* Space around the title and lines */
        """
    },
    # Decorative line above
    ui.tags.div(
        {
            "style": """
                position: absolute; 
                top: 0;  /* Position at the top of the container */
                left: 50%; 
                transform: translateX(-50%); 
                width: 220px; 
                height: 2px; 
                background-color: #6abf91; /* Green color for the line */
            """
        }
    ),
    # Title text
    ui.tags.h2(
        "Blood Panel",
        style="""
            display: inline-block; 
            font-family: 'Arial', sans-serif; 
            font-weight: 600; 
            font-size: 36px; 
            color: #333;
            margin: 20px 0; /* Space between the title and the lines */
        """
    ),
    # Decorative line below
    ui.tags.div(
        {
            "style": """
                position: absolute; 
                bottom: 0; /* Position at the bottom of the container */
                left: 50%; 
                transform: translateX(-50%); 
                width: 220px; 
                height: 2px; 
                background-color: #6abf91; /* Green color for the line */
            """
        }
    )
),

            ui.hr(),

            # Add the first section above the top HR



ui.div(
    {
        "style": """
            display: flex;
            justify-content: space-between; /* Distribute items to opposite sides */
            align-items: center; /* Vertically center items */
            padding: 10px 20px;
            gap: 20px; /* Space between elements */
        """
    },
    # Left Section: Icon and Text
    ui.div(
        {
            "style": """
                display: flex;
                align-items: flex-start; /* Align items at the start */
                gap: 20px; /* Space between icon and text */
            """
        },
        # Icon on the left
        ui.tags.i(
            {
                "class": "fas fa-user-md",  # Font Awesome Doctor icon
                "style": """
                    font-size: 200px;
                    color: #E74C3C; /* Red for emphasis */
                """
            }
        ),
        # Text Container
        ui.div(
            {
                "style": """
                    display: flex;
                    flex-direction: column; /* Stack title and text vertically */
                    align-items: flex-start; /* Align text to the left */
                    gap: 10px; /* Space between title and text */
                """
            },
            # Title
            ui.tags.p(
                "Cardiovascular disease risks",
                {
                    "style": """
                        font-size: 30px;
                        font-weight: 600;
                        color: #333;
                        margin: 0;
                    """
                }
            ),
            # Subtext below the title
            ui.tags.p(
                """
                Cardiovascular diseases (CVD) are the leading cause of mortality 
                globally. Conditions such as heart attacks, strokes, and 
                peripheral artery disease arise from underlying risk factors 
                including high blood pressure, elevated cholesterol, diabetes, 
                and smoking. In general, hypertension and blood lipids are some 
                of the most important risk factors for CVD and are often
                targeted for prevention and management. This panel provides
                some insights into your potential predisposition to these markers
                amongst others.

                
                """,
                {
                    "style": """
                        font-size: 16px;
                        color: #555;
                        margin: 0;
                        max-width: 400px; /* Limit text width for readability */
                        line-height: 1.5; /* Improve readability with line spacing */
                    """
                }
            )
        )
    ),
    # Right Section: Output Widget
    ui.div(
        {
            "style": """
                display: flex;
                justify-content: flex-end; /* Push widget to the right */
                max-width: 400px; /* Optional: Set a max width for balance */
            """
        },
        output_widget('cvd_pie_plot')
    )
),



# First HR
ui.hr(),

# Add the second section between the two HRs
ui.div(
    {
        "style": """
            display: flex;
            justify-content: space-between; /* Push text/icon and plot to opposite sides */
            align-items: center; /* Vertically center items */
            padding: 10px 20px;
        """
    },
    # Left Section: Plot
    ui.div(
        {
            "style": """
                display: flex;
                justify-content: flex-start; /* Push plot to the left */
            """
        },
        output_widget('blood_hypertension_plot')
    ),
    # Right Section: Text and Icon
    ui.div(
        {
            "style": """
                display: flex;
                justify-content: flex-end; /* Align content to the right */
                align-items: center; /* Vertically center text and icon */
                gap: 20px; /* Space between text and icon */
            """
        },
        # Text on the left
        ui.div(
            {
                "style": """
                    display: flex;
                    flex-direction: column; /* Stack title and text vertically */
                    align-items: flex-start; /* Align text to the left */
                    gap: 10px; /* Space between title and text */
                """
            },
            # Title
            ui.tags.p(
                "Blood Pressure",
                {
                    "style": """
                        font-size: 30px;
                        font-weight: 600;
                        color: #333;
                        margin: 0;
                        text-align: left; /* Align title to the left */
                    """
                }
            ),
            # Subtext below the title
            ui.tags.p(
                """
                Blood pressure is a vital marker of cardiovascular health,
                and is the strongest predictor of cardiovascular events. 
                High blood pressure, or hypertension, increases the risk of 
                heart disease, stroke, and other complications. Luckily,
                lifestyle changes and medications can help manage blood pressure

                """,
                {
                    "style": """
                        font-size: 16px;
                        color: #555;
                        margin: 0;
                        max-width: 400px; /* Limit text width for readability */
                        line-height: 1.5; /* Improve readability with line spacing */
                        text-align: left; /* Align text to the left */
                    """
                }
            )
        ),
        # Icon on the right
        ui.tags.i(
            {
                "class": "fas fa-heartbeat",  # Font Awesome Heartbeat/Blood Pressure icon
                "style": """
                    font-size: 200px;
                    color: #E74C3C; /* Blood red */
                """
            }
        )
    )
),



# Second HR
ui.hr(),


ui.div(
    {
        "style": """
            display: flex;
            justify-content: space-between; /* Push text/icon and plot to opposite sides */
            align-items: center; /* Vertically center items */
            padding: 10px 20px;
        """
    },
    # Left Section: Icon and Text
    ui.div(
        {
            "style": """
                display: flex;
                align-items: flex-start; /* Align content at the start */
                gap: 20px; /* Space between icon and text */
            """
        },
        # Blood Icon
        ui.tags.i(
            {
                "class": "fas fa-tint",  # Font Awesome Drop/Vial icon
                "style": """
                    font-size: 200px;
                    color: #E74C3C; /* Blood red */
                """
            }
        ),
        # Text Container
        ui.div(
            {
                "style": """
                    display: flex;
                    flex-direction: column; /* Stack title and text vertically */
                    align-items: flex-start; /* Align text and title to the left */
                    gap: 10px; /* Space between title and text */
                """
            },
            # Title
            ui.tags.p(
                "Blood Lipids",
                {
                    "style": """
                        font-size: 30px;
                        font-weight: 600;
                        color: #333;
                        margin: 0;
                        text-align: left; /* Align title to the left */
                    """
                }
            ),
            # Subtext below the title
            ui.tags.p(
                """
                Blood lipids, such as cholesterol and triglycerides, 
                play a critical role in cardiovascular health. 
                Imbalances or too large quantities
                in these lipids are major contributors to the 
                development of cardiovascular disease (CVD), which 
                encompasses conditions like heart attacks, strokes, 
                and peripheral artery disease. Too high cholesterol
                can lead to accelerated plaque buildup in arteries,
                which can eventually lead to artherosclerosis and
                heart disease.


                """,
                {
                    "style": """
                        font-size: 16px;
                        color: #555;
                        margin: 0;
                        max-width: 600px; /* Limit the width of the text block */
                        line-height: 1.5; /* Improve readability with line spacing */
                        text-align: left; /* Align text to the left */
                    """
                }
            )
        )
    ),
    # Right Section: Plot
    ui.div(
        {
            "style": """
                display: flex;
                justify-content: flex-end; /* Push plot to the right */
            """
        },
        output_widget('blood_lipid_plot')
    )
),


            ui.hr(),
            ui.h2('Other Blood markers'),
            # Keep just the two tables (Rare Genotypes & Rare Mutations)
            output_widget("other_blood_markers_plot"),
            ui.div(
                {
                    "style": """
                    padding: 40px 20px; 
                    position: relative;
                    """
                },
                ui.card(
                    {"class": "analysis-card"},
                    ui.h3("Rare variants panel", class_="card-title"),
                    # Replace ui.navset_tab(...) with ui.navs_tab(...)
                    ui.navset_tab(
                    # Replace ui.nav("Title", ...) with ui.nav_panel("Title", ...)
                    ui.nav_panel(
                        "Rare Genotypes",
                        ui.div(
                        {"style": "max-height: 400px; overflow-y: auto;"},
                        ui.output_ui("rare_genotypes")
                        )
                    ),
                    ui.nav_panel(
                        "Rare Mutations",
                        ui.div(
                        {"style": "max-height: 400px; overflow-y: auto;"},
                        ui.output_ui("rare_snps")
                        )
                    )
                    )
                )
            )

        )
    ),
    ),
    ui.nav_panel("Commpare with Friends",
                 ui.div(
            {"style": """
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                background-color: rgba(0, 0, 0, 0.7);
                padding: 2rem;
                border-radius: 10px;
                color: white;
                font-size: 2rem;
                font-weight: bold;
                text-shadow: 2px 2px 4px rgba(0,0,0,0.5);
                z-index: 2;
            """},
            ui.h1("Under Construction")
        )),
    title="SNPster",
    navbar_options=navbar_options(bg="#f8f9fa")
)















































def server(input, output, session):
    user_data = reactive.value(None)
    user_prs_data = reactive.value(None)
    full_merged_data = reactive.value(None)
    input_traits = reactive.value(None)
    rare_user_snps = reactive.value(None)

    

    @output
    @render.ui
    def background_image():
        base64_img = get_image_with_validation()
        if base64_img is None:
            return ui.div(
                {"style": """
                    position: fixed;
                    top: 0;
                    left: 0;
                    width: 100vw;
                    height: 100vh;
                    background-color: #f0f0f0;
                    z-index: 1;
                """},
                ui.p("Failed to load background image")
            )
            
        return ui.div(
            {"style": f"""
                position: fixed;
                top: 0;
                left: 0;
                width: 100vw;
                height: 100vh;
                background-image: url(data:image/jpeg;base64,{base64_img});
                background-size: cover;
                background-position: center;
                background-repeat: no-repeat;
                z-index: 1;
            """}
        )
    

    @reactive.effect
    @reactive.event(input.run_analysis)
    def run_analysis():
        
        if input.consent() is False:
                m = ui.modal(
                    "Please read and accept the terms of use to have your genetic data analyzed.",
                    title="Terms of use",
                    easy_close=True,
                    footer=None,
                )
                ui.modal_show(m)
                return
        
        if input.user_file() is not None:
            
            
            file_data = input.user_file()
            
            try:
                first_file = file_data[0]
                user_data_pandas = pd.read_csv(first_file['datapath'], sep="\t", comment="#")
                user_data_pandas.columns = ['rsid', 'chromosome', 'position', 'genotype']
                user_data.set(user_data_pandas)
            except Exception as e:
                return ui.notification_show("An error occured analyzing the data", type = 'error')

            return
            
        
        elif input.sample_genome() != "":
            file_data = input.sample_genome()
            
            
            
            try:
                user_data_pandas = pd.read_csv(
                    GENOME_PATHS[file_data], 
                    sep="\t", 
                    comment="#",
                    low_memory=False,
                    dtype={'rsid': 'string', 'chromosome': 'string', 'position': 'Int64', 'genotype': 'string'}
                )
                user_data_pandas.columns = ['rsid', 'chromosome', 'position', 'genotype']
                user_data.set(user_data_pandas)
            except Exception as e:
                print(e)
                return ui.notification_show (f"An error occured analyzing the data: {e}", type = 'error')

            return

        else:
            return ui.notification_show("No data uploaded yet, either upload data or pick a sample genome", type="warning")

            

        
        #render analysis_content



    @reactive.effect
    def calc_user_data():
        df = user_data.get()
        if df is None:
            return None
        
        data = prepare_data(gwas_data_path, df)

        polygenomic_risk_score_data = data["polygenomic_risk_score_data"]
        rare_snps_data = data["rare_snps"]
        user_prs_data.set(polygenomic_risk_score_data)
        rare_user_snps.set(rare_snps_data)
        full_merged_data.set(data["full_data"])
        new_traits = list(polygenomic_risk_score_data["trait"][polygenomic_risk_score_data["or"] == True].unique())
        traits_to_remove = [
            # These are just too weird or too specific
            'Liver injury in combined anti-retroviral and anti-tuberculosis drug-treated HIV with tuberculosis',
            'Pre-treatment pain in head and neck squamous cell carcinoma',
            'Recurrence of mild malaria attacks',
            'Recurrence of malaria infection (mild or asymptomatic)',
            'Kidney disease (late stage) in type 1 diabetes',
            'Renal function and chronic kidney disease',
            'Chronic kidney disease (severe chronic kidney disease vs normal kidney function) in type 1 diabetes',
            'Acute-on-chronic liver failure in hepatitis B',
            'Decreased low contrast letter acuity in multiple sclerosis',
            
            # Newly added traits
            'Straight vs curly hair',
            'Decreased fine motor function in Charcot-Marie-Tooth disease 1A (eating with utensils)',
            'Impaired insulin sensitivity in response to n-3 PUFA supplementation',
            'Hepatitis B surface antigen seroclearance in chronic hepatitis B infection',
            'Allergic sensitization',
            'Breast cancer in childhood cancer survivors treated with more than 10 gray radiotherapy',
            'Morning vs. evening chronotype',
            'Myocardial infarction in darapladib-treated cardiovascular disease (time to event)',
            'Lung function in never smokers (high FEV1 vs average FEV1)',
            'Lung function in never smokers (low FEV1 vs high FEV1)',
            'Angioedema in response to angiotensin-converting enzyme inhibitor and/or angiotensin receptor blocker',
            'Alzheimer’s disease and/or vascular dementia (clinical subgroup VaD+)',
            'Low HDL-cholesterol levels',
            'High fasting blood glucose',
            'Thick vs thin eyebrows',
            'Decreased sensory function in Charcot-Marie-Tooth disease type 1A',
            'Hearing loss in Charcot-Marie-Tooth disease 1A'
        ]
        
        filtered_traits = [trait for trait in new_traits if trait not in traits_to_remove]

        
        input_traits.set(filtered_traits)
    
    @reactive.effect
    def update_traits():
        traits = input_traits.get()
        if traits is None:
            return None
        
        
        ui.update_selectize("traits", choices=traits, selected='Melanoma')
    
    @output
    @render_widget
    def world_map():

        # Coordinates for superpopulation regions
        SUPERPOPULATION_REGIONS = {
            "European": {
                "lat": [71.0, 71.0, 36.0, 36.0, 71.0],
                "lon": [-10.0, 40.0, 40.0, -10.0, -10.0],
                "color": "rgba(0, 123, 255, 0.5)",  # Blue with transparency
            },
            "East Asian": {
                "lat": [55.0, 55.0, 10.0, 10.0, 55.0],
                "lon": [100.0, 140.0, 140.0, 100.0, 100.0],
                "color": "rgba(255, 0, 0, 0.5)",  # Red with transparency
            },
            "African": {
                "lat": [37.0, 37.0, -35.0, -35.0, 37.0],
                "lon": [-17.0, 51.0, 51.0, -17.0, -17.0],
                "color": "rgba(0, 255, 0, 0.5)",  # Green with transparency
            },
            "South Asian": {
                "lat": [35.0, 35.0, 5.0, 5.0, 35.0],
                "lon": [60.0, 90.0, 90.0, 60.0, 60.0],
                "color": "rgba(255, 165, 0, 0.5)",  # Orange with transparency
            },
            "Admixed American": {
                "lat": [48.0, 48.0, -55.0, -55.0, 48.0],
                "lon": [-125.0, -30.0, -30.0, -125.0, -125.0],
                "color": "rgba(128, 0, 128, 0.5)",  # Purple with transparency
            },
        }
        fig = go.Figure()

        # Get the selected superpopulation
        selected = input.user_superpopulation()
        if selected in SUPERPOPULATION_REGIONS:
            region = SUPERPOPULATION_REGIONS[selected]
            # Add the colored region for the selected superpopulation
            fig.add_trace(
                go.Scattergeo(
                    lat=region["lat"],
                    lon=region["lon"],
                    mode="lines",
                    fill="toself",
                    fillcolor=region["color"],
                    line=dict(color=region["color"], width=1),
                    name=selected,
                )
            )

        # Update map appearance
        fig.update_geos(
            projection_type="natural earth",
            showland=True,
            landcolor="rgb(217, 217, 217)",
            showcountries=True,
            countrycolor="rgb(204, 204, 204)",
        )
        fig.update_layout(
            height=400,
            margin={"r": 0, "t": 0, "l": 0, "b": 0},
        )


        return fig
    
    @output
    @render_widget
    def disease_category_radar_plot():
        prs_data = user_prs_data.get()
        if prs_data is None:
            return None
        
        prs_data = prs_data[prs_data["or"] == True]

        #group by Category and get mean of percentile
        grouped_df = prs_data.groupby("Category", dropna=True)["percentile"].mean().reset_index()

        #make radart plot

        fig = go.Figure(
            data=go.Scatterpolar(
                r=grouped_df['percentile'],
                theta=grouped_df['Category'],
                fill='toself',
                fillcolor='rgba(255,165,0,0.5)',  # Orange fill with 50% opacity
        line=dict(color='orange')         # Orange line color
            )
        )

        fig.update_layout(
            polar=dict(
                bgcolor='rgba(0,0,0,0)',  # Transparent radar background
                radialaxis=dict(
                    visible=True,
                    range=[0, 100],
                    gridcolor="black",  # Make radial grid lines black
                    linecolor="black"  # Make radial axis line black
                ),
                angularaxis=dict(
                    gridcolor="black",  # Make angular grid lines black
                    linecolor="black"  # Make angular axis line black
                )
            ),
            showlegend=False,
            margin=dict(t=40, l=0, r=0, b=0),
            paper_bgcolor="rgba(0,0,0,0)",  # Transparent figure background
            plot_bgcolor="rgba(0,0,0,0)"    # Transparent plot area
        )

        return fig




    @output
    @render_widget
    def pie_plot():

        pie_data = full_merged_data.get()
        if pie_data is None:
            return None
        
        pie_data = pie_data[pie_data['or'] == True].copy()


        # Handle missing categories
        pie_data.loc[:, "Category"] = pie_data["Category"].fillna("Unknown")

        # Aggregate data
        grouped_df = (
            pie_data
            .groupby(["risk_status_category", "Category"], dropna=False)
            .size()
            .reset_index(name="count")
        )

        good_bad_counts =(
            pie_data
            .groupby(["risk_status_category"], dropna=False)
            .size()
            .reset_index(name="count")
        )

        # Define colors
        inner_colors = {"good": "#009E73", "bad": "#D55E00"}  # Inner ring colors
        outer_colors = {
            "Autoimmune and Inflammatory Diseases": "#76c893",
            "Cancer": "#619b8a",
            "CVD and Chronic Conditions": "#4d908e",
            "Infectious Diseases": "#43aa8b",
            "Metabolic Disorders": "#3e3e5b",
            "Neurological Disorders": "#2a9d8f",
            "Other Conditions": "#e76f51",
            "Developmental and Structural Conditions": "#f4a261"
        }

        #intersect between Category and risk_status_category
        grouped_df['concat'] = grouped_df['risk_status_category'] +' - '+ grouped_df['Category']

        #change 'bad - ' to (-) and 'good - ' to (+)

        grouped_df['concat'] = grouped_df['concat'].str.replace('bad - ', '(-)')
        grouped_df['concat'] = grouped_df['concat'].str.replace('good - ', '(+)')
        
        labels = grouped_df['concat'].tolist()

        #add 'bad' and 'good' to the list

        labels.extend(['Predisposition', 'Protective', 'SNPs'])



        parents = []
        for element in labels:
            if element == 'Predisposition':
                parents.append('SNPs')
                continue
            elif element == 'Protective':
                parents.append('SNPs')
                continue
            if element == 'SNPs':
                parents.append('')
                continue
            
            #check if good in element
            if '(+)' in element:
                parents.append('Protective')
            elif '(-)' in element:
                parents.append('Predisposition')

        values = grouped_df['count'].tolist()
        summed_counts = good_bad_counts['count'].tolist()
        sum_of_summed = sum(summed_counts)

        #extend values
        summed_counts.extend([sum_of_summed])

        # Then, extend `values` with the updated `summed_counts`
        values.extend(summed_counts)

        # Define category colors
        category_colors = {
            "Autoimmune and Inflammatory Diseases": "#4CAF50",  # Green
            "Cancer": "#F44336",                                # Red
            "CVD and Chronic Conditions": "#2196F3",           # Blue
            "Infectious Diseases": "#FF9800",                  # Orange
            "Metabolic Disorders": "#9C27B0",                  # Purple
            "Neurological Disorders": "#3F51B5",               # Deep Blue
            "Other Conditions": "#607D8B",                     # Grayish Blue
            "Developmental and Structural Conditions": "#00BCD4" # Cyan
        }

        # Debugging: print grouped_df to ensure replacement worked

        # Assign colors to labels
        outer_color = {}
        for label in labels:
            if "(-)" in label:
                # Remove "(-)" and handle category
                category = label.replace('(-)', '').strip()
            elif "(+)" in label:
                # Remove "(+)" and handle category
                category = label.replace('(+)','').strip()
            else:
                # Use the label directly as the category
                category = label.strip()
            
            # Handle missing categories with a default color
            outer_color[label] = category_colors.get(category, "#000000")  # Default to black for missing categories

            # Handle specific cases for "Predisposition," "Protective," and "SNPs"
            if label == "Predisposition":
                outer_color[label] = "#D55E00"  # Orange-red
            elif label == "Protective":
                outer_color[label] = "#009E73"  # Green
            elif label == "SNPs":
                outer_color[label] = "#EEEEEE"  # Light Gray


        # Assign colors dynamically
        colors = [outer_color[label] for label in labels]



        # Create the sunburst plot
        fig = go.Figure(go.Sunburst(
            labels=labels,
            parents=parents,
            values=values,
            marker=dict(colors=colors),  # Assign the colors
            branchvalues="total"  # Ensure values propagate to parent nodes
        ))

        fig.update_layout(
            title="Protective: (+), Predisposition: (-)",
            margin=dict(t=40, l=0, r=0, b=0),
            paper_bgcolor="rgba(0,0,0,0)",  # Transparent background for the entire figure
            plot_bgcolor="rgba(0,0,0,0)"  # Transparent background for the plot area
        )




        return fig

    
        

    
    @output
    @render_widget
    def disease_risk_plot():
        prs_data = user_prs_data.get()
        selected_traits = input.traits()

        if prs_data is None or selected_traits is None:
            return

        prs_data = prs_data[prs_data["trait"].isin([selected_traits])]
        prs_data = prs_data[prs_data["or"] == True]
        

        percentile_value = prs_data["percentile"].values[0]

        # -------------------------------------------------------------------------
        # 2. Convert that percentile to a z-score
        # -------------------------------------------------------------------------
        z_score = prs_data["final_risk"].values[0]

        # -------------------------------------------------------------------------
        # 3. Create x-values and y-values for the standard normal distribution
        # -------------------------------------------------------------------------
        x = np.linspace(-4, 4, 1000)
        y = norm.pdf(x)
        df_all = pd.DataFrame({"x": x, "y": y})

        # We'll split the data into two parts for coloring:
        #   left side (x <= z_score) and right side (x >= z_score)
        df_left = df_all[df_all["x"] <= z_score]
        df_right = df_all[df_all["x"] >= z_score]

        # -------------------------------------------------------------------------
        # 4. Plot with Plotly Express
        # -------------------------------------------------------------------------
        # First area: left side (red)
        fig = px.area(
            df_left,
            x="x",
            y="y",
            labels={"x": "Z-score", "y": "Density"},
            color_discrete_sequence=["red"],    # force it to be red
            title=f"Your risk compared to the general population — {percentile_value:.1f}th Percentile"
        )

        # Add the second area: right side (green) as another trace
        right_trace = px.area(
            df_right,
            x="x",
            y="y",
            color_discrete_sequence=["green"]  # force it to be green
        ).data[0]

        # Update hoverinfo for both traces to skip hover interaction
        fig.data[0].update(hoverinfo="skip", hovertemplate=None)  # For the first trace (red)
        right_trace.update(hoverinfo="skip", hovertemplate=None)  # For the second trace (green)

        # Add the second trace to the figure
        fig.add_trace(right_trace)

        # Add a vertical dashed line at the z-score
        fig.add_vline(x=z_score, line_width=2, line_dash="dash", line_color="black")

        # Optionally, add an annotation showing the exact percentile
        fig.add_annotation(
            x=z_score,
            y=np.max(y) * 0.9,   # place the annotation near the top of the distribution
            text=f"{percentile_value:.1f}%",
            showarrow=True,
            arrowhead=2
        )

        # Adjust layout as needed
        fig.update_layout(
            xaxis=dict(range=[-4, 4]),
            yaxis=dict(range=[0, max(y) * 1.1]),  # small padding on top
        )

        # Remove y-axis and background
        fig.update_yaxes(visible=False, showticklabels=False)
        fig.update_layout(
            paper_bgcolor="rgba(0, 0, 0, 0)",  
            plot_bgcolor="rgba(0, 0, 0, 0)",
            hovermode=False  # Disable hover globally
        )

        return fig
    
    @output
    @render_widget
    def full_disease_overview():
        prs_data = user_prs_data.get()
        if prs_data is None:
            return None
        

        traits_to_remove = [
            # These are just too weird or too specific
            'Liver injury in combined anti-retroviral and anti-tuberculosis drug-treated HIV with tuberculosis',
            'Pre-treatment pain in head and neck squamous cell carcinoma',
            'Recurrence of mild malaria attacks',
            'Recurrence of malaria infection (mild or asymptomatic)',
            'Kidney disease (late stage) in type 1 diabetes',
            'Renal function and chronic kidney disease',
            'Chronic kidney disease (severe chronic kidney disease vs normal kidney function) in type 1 diabetes',
            'Acute-on-chronic liver failure in hepatitis B',
            'Decreased low contrast letter acuity in multiple sclerosis',
            
            # Newly added traits
            'Straight vs curly hair',
            'Decreased fine motor function in Charcot-Marie-Tooth disease 1A (eating with utensils)',
            'Impaired insulin sensitivity in response to n-3 PUFA supplementation',
            'Hepatitis B surface antigen seroclearance in chronic hepatitis B infection',
            'Allergic sensitization',
            'Breast cancer in childhood cancer survivors treated with more than 10 gray radiotherapy',
            'Morning vs. evening chronotype',
            'Myocardial infarction in darapladib-treated cardiovascular disease (time to event)',
            'Lung function in never smokers (high FEV1 vs average FEV1)',
            'Lung function in never smokers (low FEV1 vs high FEV1)',
            'Angioedema in response to angiotensin-converting enzyme inhibitor and/or angiotensin receptor blocker',
            'Alzheimer’s disease and/or vascular dementia (clinical subgroup VaD+)',
            'Low HDL-cholesterol levels',
            'High fasting blood glucose',
            'Thick vs thin eyebrows',
            'Decreased sensory function in Charcot-Marie-Tooth disease type 1A',
            'Hearing loss in Charcot-Marie-Tooth disease 1A'
        ]

        prs_data = prs_data[~prs_data["trait"].isin(traits_to_remove)]




        # Normalize the percentiles to center at 50
        prs_data = prs_data[prs_data["or"] == True]

        prs_data["percentile_normalized"] = prs_data["percentile"] - 50



        # Create the bar plot
        fig = px.bar(
            prs_data,
            x="percentile_normalized",  # Normalized percentile on the x-axis
            y="trait",                  # Trait on the y-axis
            orientation="h",            # Horizontal bars
            labels={
                "percentile_normalized": "Deviation from Population Average",
                "trait": "Trait",

            },
            title="Disease Overview by Percentile",
            text="percentile",          # Show the original percentile as text on the bars
        )

        # Format the layout for better readability
        fig.update_layout(
            xaxis=dict(
                title="Deviation from Population Average",
                zeroline=True, zerolinewidth=2, zerolinecolor="black",
            ),
            yaxis=dict(title="Trait"),
            title=dict(x=0.5, font=dict(size=18, family="Arial")),  # Center title
            plot_bgcolor="rgba(0, 0, 0, 0)",    # Transparent background
            paper_bgcolor="rgba(0, 0, 0, 0)",   # Transparent paper
        )

        # Add custom styling for bar text
        fig.update_traces(
            texttemplate='%{text:.1f}%',  # Display original percentile with 1 decimal
            textposition="outside",      # Position text outside the bars
            marker_color=prs_data["percentile_normalized"].apply(
                lambda x: "red" if x > 0 else "green"  # Green for above 50, red for below
            ),
        )


        return fig
    
    @output
    @render_widget
    def full_trait_overview():
        prs_data = user_prs_data.get()
        if prs_data is None:
            return None
        
        # Normalize the percentiles to center at 50
        prs_data = prs_data[prs_data["or"] == False]

        # Create the bar plot
        fig = px.bar(
            prs_data,
            x="final_risk",  # Normalized percentile on the x-axis
            y="trait",                  # Trait on the y-axis
            orientation="h",            # Horizontal bars
            labels={
                "final_riskl": "Deviation from 50th Percentile (Right: Increase, Left: Decrease)",
                "trait": "Trait"
            },
            title="Trait Overview by Percentile",
            text="percentile",          # Show the original percentile as text on the bars
        )

        # Format the layout for better readability
        fig.update_layout(
            xaxis=dict(
                title="Deviation from 50th Percentile (Right: Increase, Left: Decrease)",
                zeroline=True, zerolinewidth=2, zerolinecolor="black"
            ),
            yaxis=dict(title="Trait"),
            title=dict(x=0.5, font=dict(size=18, family="Arial")),  # Center title
            plot_bgcolor="rgba(0, 0, 0, 0)",    # Transparent background
            paper_bgcolor="rgba(0, 0, 0, 0)",   # Transparent paper
        )

        # Define colors for increase (blue) and decrease (orange)
        colors = prs_data["final_risk"].apply(
            lambda x: "blue" if x > 0 else "orange"  # Blue for increase, Orange for decrease
        )

        # Add custom styling for bar text and colors
        fig.update_traces(
            texttemplate='%{text:.1f}%',  # Display original percentile with 1 decimal
            textposition="outside",      # Position text outside the bars
            marker_color=colors,         # Set colors based on the direction
        )

        return fig
    
    @output
    @render_widget
    def power_vs_endurance():
        # Sample data
        sprint = 75
        endurance = 50

        fig = go.Figure()

        fig.add_trace(go.Indicator(
            mode="gauge+number",
            value=sprint,
            title={"text": "Sprint Proclivity"},
            gauge={
                "axis": {"range": [0, 100]},
                "bar": {"color": "#FF5733"},
                "steps": [
                    {"range": [0, 50], "color": "#ffe6e6"},
                    {"range": [50, 100], "color": "#ffcccc"}
                ],
            }
        ))

        fig.add_trace(go.Indicator(
            mode="gauge+number",
            value=endurance,
            title={"text": "Endurance Proclivity"},
            gauge={
                "axis": {"range": [0, 100]},
                "bar": {"color": "#33C1FF"},
                "steps": [
                    {"range": [0, 50], "color": "#e6f7ff"},
                    {"range": [50, 100], "color": "#ccefff"}
                ],
            },
            domain={"x": [0.5, 1], "y": [0, 1]}  # Position the second gauge
        ))

        fig.update_layout(
            title="",
            grid={"rows": 1, "columns": 2}
        )

        return fig

        

    @output
    @render.ui
    def sports_cards():
        # Obtain the user data (already merged with sports_variants)
        full_user_upload = user_data.get()
        
        # Load additional sports variant data
        sports_variants = pd.read_csv(
            get_background_data_path("sports_snps.txt"),
            sep="\t"
        )
        

        # Merge your user data with the sports variant info
        full_user_upload = pd.merge(
            full_user_upload, 
            sports_variants, 
            left_on="rsid", 
            right_on="SNP", 
            how="inner"
        )
        
        # Create a list to store the card elements
        cards = []
        for _, row in full_user_upload.iterrows():
            card = ui.card(
                {"class": "sport-card"},
                
                # 1) Title: Gene
                ui.h3(row['gene'], class_="sport-title"),
                
                # 2) Below Title: Personal Genotype
                ui.div(
                    {"class": "genotype"},
                    f"Genotype: {row.get('genotype', 'N/A')}"
                ),
                
                # 3) Same Line: Power Allele and Endurance Allele
                ui.div(
                    {"style": "margin-top: 10px; margin-bottom: 10px;"},
                    ui.span(f"Power Allele: {row.get('power_allele', 'N/A')}"),
                    ui.span(" \u00A0|\u00A0 "),  # Non-breaking space + vertical bar
                    ui.span(f"Endurance Allele: {row.get('endurance_allele', 'N/A')}")
                ),
                
                # 4) Variant Description (Comment)
                ui.p(row.get('Comment', ''))
            )
            
            cards.append(card)
        
        # Return the cards in a responsive column wrap (3 cards per row)
        return ui.layout_column_wrap(
            width=1/3,
            *cards
        )

    @output
    @render.ui

    def top_results_plot():
        if user_data.get() is None:
            return None
        
        data = user_data.get()
        
        """ This function shows these five interesting snps 
        i3003626(rs333) HIV resistance (CCR5) -.
        rs1815739 Sprint performance (ACTN3) C.
        rs12913832 Blue eyes G .
        rs1042522 Increased Lifespan (TP53). C
        rs53576 Empathy and social behavior (OXTR). A
        rs1805007 Red hair T
        rs1805008 Red har T
        rs7221412 nightowl chronotype G"""

        #change data i3003626 to rs333
        data = data.replace("i3003626", "rs333")

        rsid = ["rs333", "rs1815739", "rs12913832", "rs1042522", "rs53576", "rs1805007", "rs1805008", "rs7221412"]
        descriptions = ["HIV resistance (CCR5)", "Sprint performance (ACTN3)", "Blue eyes", "Cancer risk modulation (TP53)", "Empathy and social behavior (OXTR)","Red hair color","Red hair color", "Nightowl chronotype"]
        effect_allele = ["-", "C", "G", "C", "A", "T", "T", "G"]


        
        
        mutation_df = pd.DataFrame(columns=["snps", "descriptions", "effect_allele"])
        mutation_df["rsid"] = rsid
        mutation_df["descriptions"] = descriptions
        mutation_df["effect_allele"] = effect_allele

        #merge data and mutation_df on snps and rsid

        mutation_df_merged = pd.merge(mutation_df, data, on="rsid", how="inner")

        #count number of occurences of the effect_allele in genotype

        mutation_df_merged['occurences'] = mutation_df_merged.apply(
            lambda row: row['genotype'].count(row['effect_allele']), axis=1
        )
        #remove snps column

        mutation_df_merged = mutation_df_merged.drop("snps", axis=1)
        


        def color_rows(row):
            count = row['occurences']  # Use the calculated occurrences
            if count == 2:
                return ["background-color: #D4EDDA; color: #155724;"] * len(row)  # Soft green
            elif count == 1:
                return ["background-color: #FFF3CD; color: #856404;"] * len(row)  # Soft yellow
            elif count == 0:
                return ["background-color: #F8D7DA; color: #721C24;"] * len(row)  # Soft red
            return [""] * len(row)

        # Style the DataFrame
        styled_table = (
            mutation_df_merged.style
            .apply(color_rows, axis=1)  # Apply row-wise coloring based on occurrences
            .set_properties(**{"text-align": "center", "font-size": "14px"})  # Modern font size
            .hide(axis="index")  # Hide the index
        )

        # Return the sty
        return ui.HTML(styled_table.to_html())
    
    @output
    @render.ui
    def top_results_plot_2():
        trait = input.traits()
        user_data = full_merged_data.get()

        if user_data is None or trait is None:
            return
        if trait not in user_data['trait'].values:
            return

        trait_data = user_data[user_data['trait'] == trait]

        columns_to_remove = ['Unnamed: 0', 'or', 'trait_direction', 'population_score', 
                            'Category', 'chromosome', 'position', 'other_allele_frequency', 
                            'risk_status_category', 'final_risk', 'risk_status']
            
        order = ['trait', 'odds_ratio', 'confidence_interval', 'rsid', 'pmid', 
                'gene', 'risk_allele', 'personal_genotype', 'risk_allele_frequency']
            
        trait_data = trait_data.drop(columns=columns_to_remove)
        trait_data = trait_data[order]
        # Order by odds_ratio
        trait_data = trait_data.sort_values(by="odds_ratio", ascending=False)

        # Apply conditional formatting
        def highlight_risk_status(row):
            # Count occurrences of risk_allele in personal_genotype
            risk_allele_count = row['personal_genotype'].count(row['risk_allele'])
            if risk_allele_count == 2:
                return ['background-color: #f8d7da'] * len(row)  # Mild red
            elif risk_allele_count == 1:
                return ['background-color: #fff3cd'] * len(row)  # Mild yellow
            else:
                return [''] * len(row)

        styled_table = (
            trait_data.style
            .apply(highlight_risk_status, axis=1)
            .set_table_styles([
                # Table-wide styles
                {'selector': 'table', 'props': [
                    ('width', '100%'),
                    ('border-collapse', 'collapse'),
                    ('margin', '0 auto'),
                    ('font-family', 'Arial, sans-serif')
                ]},
                # Header styles
                {'selector': 'thead', 'props': [
                    ('background-color', '#4CAF50'),
                    ('color', 'white'),
                    ('font-weight', 'bold'),
                    ('text-transform', 'uppercase')
                ]},
                # Header cell styles
                {'selector': 'th', 'props': [
                    ('padding', '12px 15px'),
                    ('text-align', 'left'),
                    ('border-bottom', '2px solid #ddd')
                ]},
                # Data cell styles
                {'selector': 'td', 'props': [
                    ('padding', '10px 15px'),
                    ('border-bottom', '1px solid #ddd')
                ]},
                # Alternating row colors
                {'selector': 'tbody tr:nth-child(even)', 'props': [
                    ('background-color', '#f8f9fa')
                ]},
                # Hover effect
                {'selector': 'tbody tr:hover', 'props': [
                    ('background-color', '#f2f2f2'),
                    ('transition', 'background-color 0.3s ease')
                ]}
            ])
            .set_properties(**{
                'text-align': 'left',
                'font-size': '14px',
                'white-space': 'nowrap'
            })
        )

        # Convert to HTML and wrap in ui.HTML
        return ui.HTML(styled_table.to_html())



        
    @output
    @render.ui
    def rare_snps():
        if rare_user_snps.get() is None:
            return None
        
        data = rare_user_snps.get()

                # Reset the index if it's not unique
        data = data.reset_index(drop=True)

        # Ensure columns are unique
        data.columns = pd.Index(data.columns).drop_duplicates()

        data = data.reset_index(drop=True)

        data['genotype_probability'] = data['genotype_probability']*100
        #make string with %

        data['genotype_probability'] = data['genotype_probability'].apply(lambda x: f"{x:.5f}%")
        
        data = data.drop(columns=["genotype_probability"])

        #order of columns = rsid, gene, risk_allele, personal_genotype, risk_allele_frequency, other_allele_frequency

        column_order = ["rsid", "gene", "risk_allele", "personal_genotype", "risk_allele_frequency", "other_allele_frequency"]

        # Reorder columns
        data = data[column_order]

        # Style the DataFrame
        styled_table = (
            data.style
            .set_table_styles([
                # Table-wide styles
                {'selector': 'table', 'props': [
                    ('width', '100%'),
                    ('border-collapse', 'collapse'),
                    ('margin', '0 auto')
                ]},
                # Header styles
                {'selector': 'thead', 'props': [
                    ('background-color', '#4CAF50'),
                    ('color', 'white'),
                    ('font-weight', 'bold')
                ]},
                # Header cell styles
                {'selector': 'th', 'props': [
                    ('padding', '12px 15px'),
                    ('text-align', 'left'),
                    ('border-bottom', '2px solid #ddd')
                ]},
                # Data cell styles
                {'selector': 'td', 'props': [
                    ('padding', '10px 15px'),
                    ('border-bottom', '1px solid #ddd')
                ]},
                # Alternating row colors
                {'selector': 'tbody tr:nth-child(even)', 'props': [
                    ('background-color', '#f5f5f5')
                ]},
                # Hover effect
                {'selector': 'tbody tr:hover', 'props': [
                    ('background-color', '#f0f0f0')
                ]}
            ])
            .set_properties(**{
                'text-align': 'left',
                'font-size': '13px',
                'white-space': 'nowrap'  # Prevents text wrapping in cells
            })
        )
        
        
        # Convert to HTML and wrap in ui.HTML
        return ui.HTML(styled_table.to_html())
    
    @output
    @render.ui
    def rare_genotypes():
        if rare_user_snps.get() is None:
            return None
        
        data = rare_user_snps.get()

                # Reset the index if it's not unique
        data = data.reset_index(drop=True)

        # Ensure columns are unique
        data.columns = pd.Index(data.columns).drop_duplicates()

        data = data.reset_index(drop=True)

        data['genotype_probability'] = data['genotype_probability']*100
        #make string with %

        data['genotype_probability'] = data['genotype_probability'].apply(lambda x: f"{x:.5f}%")

        data = data.drop(columns=["risk_allele", "risk_allele_frequency", "other_allele_frequency", "Category"])

        #remove dublicate rsids
        data = data.drop_duplicates(subset=['rsid'])

        #order by genotype_probability ascending

        data = data.sort_values(by="genotype_probability", ascending=True)
        
        # Style the DataFrame
        styled_table = (
            data.style
            .set_table_styles([
                # Table-wide styles
                {'selector': 'table', 'props': [
                    ('width', '100%'),
                    ('border-collapse', 'collapse'),
                    ('margin', '0 auto')
                ]},
                # Header styles
                {'selector': 'thead', 'props': [
                    ('background-color', '#4CAF50'),
                    ('color', 'white'),
                    ('font-weight', 'bold')
                ]},
                # Header cell styles
                {'selector': 'th', 'props': [
                    ('padding', '12px 15px'),
                    ('text-align', 'left'),
                    ('border-bottom', '2px solid #ddd')
                ]},
                # Data cell styles
                {'selector': 'td', 'props': [
                    ('padding', '10px 15px'),
                    ('border-bottom', '1px solid #ddd')
                ]},
                # Alternating row colors
                {'selector': 'tbody tr:nth-child(even)', 'props': [
                    ('background-color', '#f5f5f5')
                ]},
                # Hover effect
                {'selector': 'tbody tr:hover', 'props': [
                    ('background-color', '#f0f0f0')
                ]}
            ])
            .set_properties(**{
                'text-align': 'left',
                'font-size': '13px',
                'white-space': 'nowrap'  # Prevents text wrapping in cells
            })
        )
        
        
        # Convert to HTML and wrap in ui.HTML
        return ui.HTML(styled_table.to_html())
            

        


        

        
    
    @output
    @render.text
    def trait_description():
        trait = input.traits()
        user_data = user_prs_data.get()
 
        if user_data is None or trait is None:
            return
        if trait not in user_data['trait'].values:
            return
        pmids = user_data['pmid'][user_data['trait'] == trait]
        #make pandas object string
        pmids = pmids.to_string(index=False)
        n_pmids = len(pmids.split(","))
        try:
            percentile_number = user_data['percentile'][user_data['trait'] == trait].values[0]
        except IndexError:
            percentile_number = user_data['percentile'][user_data['trait'] == trait]
        
        risk_snps = user_data['rsid'][user_data['trait'] == trait].values[0]
        n_user_risk_snps = len(risk_snps.split(","))
 

        
        report_text = (
        f"Based on your genome analysis, you are in the {percentile_number:.1f}th percentile "
        f"for developing {trait}, meaning that your risk is higher than {percentile_number:.1f}% of the population. \n\n"
        f"Analysis Details:\n"
        f"• {n_user_risk_snps} risk SNPs identified\n"
        f"• {n_pmids} research articles referenced\n\n"

    )
        return report_text
    
    @output
    @render_widget
    def blood_hypertension_plot():
        # Blood Pressure–Related Traits (approximate reference values)

        merged_data = user_prs_data.get()

        merged_data = merged_data[merged_data['or']==False]


        blood_pressure_dict = {
            "Diastolic blood pressure": 80,
            "Mean arterial pressure": 93, 
            "Pulse pressure": 40, 
            "Response to beta blocker use in hypertension (systolic blood pressure)": 120,
            "Systolic blood pressure": 120,
            "Systolic blood pressure (dietary potassium intake interaction)": 120
        }

        #subset merged_data by blood pressure traits

        merged_data = merged_data[merged_data['trait'].isin(blood_pressure_dict.keys())]
    

        #plot barplots with "final_risk" on y-axis and "trait" on x-axis
        aggregated = merged_data.groupby("trait", as_index=False)["final_risk"].sum()
        
        
        fig = px.bar(
            aggregated,
            x="trait",
            y="final_risk",
            labels={
                "trait": "Trait",
                "final_risk": "mmHg"
            },
            color = "trait",
            title="Blood Pressure–Related Traits",
        )

        #make background transparent
        fig.update_layout(
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Transparent background
            plot_bgcolor="rgba(0, 0, 0, 0)",
               legend=dict(
                    orientation="h",  # Horizontal legend
                    y=-0.2,           # Position below the plot (-ve moves it below)
                    x=0.5,            # Centered horizontally
                    xanchor="center", # Anchor the legend's horizontal center to the plot center
                    yanchor="top"     # Anchor the legend's top to the plot's bottom
                )   # Transparent plot area
        )

        fig.update_xaxes(showticklabels=False)


        return fig



        
    
    @output
    @render_widget
    def blood_lipid_plot():
        # Blood Pressure–Related Traits (approximate reference values)

        merged_data = user_prs_data.get()

        merged_data = merged_data[merged_data['or']==False]


        blood_lipids_dict = {
            "HDL cholesterol": 60,  # Often considered protective level
            "HDL cholesterol levels": 60, 
            "High density lipoprotein cholesterol levels": 60,
            "LDL cholesterol": 100,  # Near-optimal
            "LDL cholesterol levels": 100,
            "Low density lipoprotein cholesterol levels": 100,
            "Response to fenofibrate (LDL cholesterol levels)": 100,
            "Response to fenofibrate (triglyceride levels)": 150, # Upper normal for TG
            "Triglyceride levels": 150,
            "Triglycerides": 150
        }

        lipid_mapping_dict = {
            # LDL cholesterol
            "LDL cholesterol levels": "LDL cholesterol",
            "LDL cholesterol": "LDL cholesterol",
            "Low density lipoprotein cholesterol levels": "LDL cholesterol",
            
            # Triglycerides
            "Triglyceride levels": "Triglycerides",
            "Triglycerides": "Triglycerides",
            
            # HDL cholesterol
            "HDL cholesterol levels": "HDL cholesterol",
            "HDL cholesterol": "HDL cholesterol",
            "High density lipoprotein cholesterol levels": "HDL cholesterol",
            
            
            # Response to fenofibrate
            "Response to fenofibrate (LDL cholesterol levels)": "Response to fenofibrate (LDL cholesterol levels)",
            "Response to fenofibrate (triglyceride levels)": "Response to fenofibrate (LDL cholesterol levels)"
        }

        #subset merged_data by blood pressure traits
        merged_data = merged_data[merged_data['trait'].isin(blood_lipids_dict.keys())]
    
        merged_data["trait"] = merged_data["trait"].replace(lipid_mapping_dict)

        #plot barplots with "final_risk" on y-axis and "trait" on x-axis

        aggregated = merged_data.groupby("trait", as_index=False)["final_risk"].sum()

        fig = px.bar(
            aggregated,
            x="trait",
            y="final_risk",
            labels={
                "trait": "Trait",
                "final_risk": "Unit Increase"
            },
            color="trait",
            title="Blood Lipid Traits"
        )

        #make background transparent
        fig.update_layout(
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Transparent background
            plot_bgcolor="rgba(0, 0, 0, 0)",   # Transparent plot area
            legend=dict(
                    orientation="h",  # Horizontal legend
                    y=-0.2,           # Position below the plot (-ve moves it below)
                    x=0.5,            # Centered horizontally
                    xanchor="center", # Anchor the legend's horizontal center to the plot center
                    yanchor="top"     # Anchor the legend's top to the plot's bottom
                )   # Transparent plot area
        )

        fig.update_xaxes(showticklabels=False)

        #remove x axis text

   

        return fig
    

    @output
    @render_widget
    def cvd_pie_plot():

        merged_data = full_merged_data.get()
        merged_data = merged_data[merged_data['or']==False]

        blood_lipids_dict = {
            "HDL cholesterol": 60,  # Often considered protective level
            "HDL cholesterol levels": 60, 
            "High density lipoprotein cholesterol levels": 60,
            "LDL cholesterol": 100,  # Near-optimal
            "LDL cholesterol levels": 100,
            "Low density lipoprotein cholesterol levels": 100,
            "Response to fenofibrate (LDL cholesterol levels)": 100,
            "Response to fenofibrate (triglyceride levels)": 150, # Upper normal for TG
            "Triglyceride levels": 150,
            "Triglycerides": 150
        }

        blood_pressure_dict = {
            "Diastolic blood pressure": 80,
            "Mean arterial pressure": 93, 
            "Pulse pressure": 40, 
            "Response to beta blocker use in hypertension (systolic blood pressure)": 120,
            "Systolic blood pressure": 120,
            "Systolic blood pressure (dietary potassium intake interaction)": 120
        }

        all_keys = list(blood_pressure_dict.keys()) + list(blood_lipids_dict.keys())

        merged_data = merged_data[merged_data['trait'].isin(all_keys)]

        #group by risk_status_category

        good_bad_counts =(
            merged_data
            .groupby(["risk_status_category"], dropna=False)
            .size()
            .reset_index(name="count")
        )

        #make good = protective and bad = predisposition

        good_bad_counts['risk_status_category'] = good_bad_counts['risk_status_category'].replace(
            {
                "good": "Protective",
                "bad": "Predisposition"
            }
        )

     

        fig = px.pie(
            good_bad_counts,
            names="risk_status_category",
            values="count",
            title="Cardiovascular Disease Variant Distribution"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Transparent background
            plot_bgcolor="rgba(0, 0, 0, 0)",   # Transparent plot area
            legend=dict(
                    orientation="h",  # Horizontal legend
                    y=-0.2,           # Position below the plot (-ve moves it below)
                    x=0.5,            # Centered horizontally
                    xanchor="center", # Anchor the legend's horizontal center to the plot center
                    yanchor="top"     # Anchor the legend's top to the plot's bottom
                )   # Transparent plot area
        )


        return fig
    
    @output
    @render_widget
    def other_blood_markers_plot():

        merged_data = full_merged_data.get()
        merged_data = merged_data[merged_data['or']==False]


        other_blood_metrics_traits = [
            "Asymmetrical dimethylarginine levels",
            "Bilirubin levels",
            "Creatine kinase levels in statin users",
            "Creatinine clearance",
            "Estimated glomerular filtration rate",
            "Gamma glutamyl transferase levels",
            "Glycated hemoglobin levels",
            "Hematocrit",
            "Hemoglobin levels",
            "Hepcidin levels",
            "kallikrein-11 levels",
            "L-arginine levels",
            "Mean corpuscular hemoglobin",
            "Platelet count",
            "Plasma free asparagine levels",
            "Plasma free proline levels",
            "Plasma plasminogen activator levels",
            "Plasma thyroid-stimulating hormone levels",
            "Red blood cell count",
            "Resistin levels",
            "Serum alkaline phosphatase levels",
            "Serum creatinine levels",
            "Serum galactose-deficient IgA1 levels",
            "Serum IgE levels",
            "Serum parathyroid hormone levels",
            "Serum uric acid levels",
            "Symmetrical dimethylarginine levels",
            "Thyroid stimulating hormone levels",
            "Uric acid clearance",
            "White blood cell count"
        ]

        merged_data = merged_data[merged_data['trait'].isin(other_blood_metrics_traits)]
        merged_data = merged_data[merged_data['trait']!="Platelet count"]
        #plot barplots with "final_risk" on y-axis and "trait" on x-axis

        aggregated = merged_data.groupby("trait", as_index=False)["final_risk"].sum()

        

        fig = px.bar(
            aggregated,
            x="trait",
            y="final_risk",
            labels={
                "trait": "Trait",
                "final_risk": "Unit Change"
            },
            color="trait",
            title="Other Blood Metrics"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0, 0, 0, 0)",  # Transparent background
            plot_bgcolor="rgba(0, 0, 0, 0)",
            showlegend=False   # Transparent plot area
        )


        return fig



            



app = App(app_ui, server)
