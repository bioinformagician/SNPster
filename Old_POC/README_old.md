# SNPster - Genomic Analysis Platform

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Python-green.svg)](https://shiny.rstudio.com/py/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Proof of Concept (POC) for the upcoming SNPster web application**

A comprehensive web application for personal genomic analysis, providing insights into disease predispositions, inherited traits, and athletic performance potential through polygenic risk scoring (PGS) and highly predictive SNPs. This repository contains the proof-of-concept implementation that demonstrates the core functionality and user interface design for the full production application.

![SNPster Dashboard](web_app/www/gene2.webp)

## Features

- **Polygenic Risk Scoring**: Advanced statistical analysis of genetic variants
- **Disease Risk Assessment**: Comprehensive analysis of health predispositions
- **Athletic Performance Analysis**: Genetic insights into power vs. endurance traits
- **Blood Panel Insights**: Analysis of cardiovascular and metabolic markers
- **Interactive Visualizations**: Radar plots, sunburst charts, and distribution analyses
- **Privacy-Focused**: Secure, anonymous analysis with no personally identifiable information
- **Cross-Platform**: Works on Windows, macOS, and Linux

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Scientific Background](#scientific-background)
- [Data Sources](#data-sources)
- [Contributing](#contributing)
- [Contact](#contact)

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Dependencies

Install the required packages using:

```bash
pip install -r requirements.txt
```

Or install manually:
```bash
pip install shiny pandas numpy scipy plotly htmltools shinywidgets
```

### Quick Start

1. **Clone the repository:**
```bash
git clone https://github.com/bioinformagician/SNIPster.git
cd SNIPster/genomics_project/travel
```

2. **Verify setup:**
```bash
python verify_setup.py
```

3. **Run the application:**
```bash
cd web_app
python -c "from shiny_ui import app; app.run()"
```

4. **Open your browser to:** `http://localhost:8000`

## Usage

### Uploading Your Data

SNPster accepts genetic data files from popular direct-to-consumer testing companies:

- **23andMe** raw data files (.txt)
- **AncestryDNA** raw data files (.txt)
- **MyHeritage** raw data files (.csv)
- **FamilyTreeDNA** raw data files (.csv)

### Sample Analysis

Try the platform with included sample genomes:
- `genome_1.txt` - European ancestry sample
- `genome_2.txt` - East Asian ancestry sample  
- `genome_3.txt` - African ancestry sample

### Key Analysis Sections

1. **Genomic Analysis**
   - Overall disease category predisposition
   - SNP sunburst overview
   - Detailed disease risk analysis with literature references

2. **Sports Performance**
   - Power vs. endurance genetic traits
   - Athletic performance indicators
   - Muscle fiber type predictions

3. **Blood Panel**
   - Cardiovascular risk markers
   - Blood pressure predispositions
   - Lipid metabolism insights

## Project Structure

```
SNIPster/
├── README.md                          # This file
├── README_PATHS.md                    # Path configuration documentation
├── SOLUTION_SUMMARY.md               # Technical implementation summary
├── example_usage.py                  # Usage examples
├── verify_setup.py                   # Setup verification script
├── gwas_data_filtered_phased_low_filt.csv  # GWAS dataset (8.5MB)
├── background_data/
│   └── sports_snps.txt               # Athletic performance SNPs
├── genomes/                          # Sample genome files
│   ├── genome_1.txt
│   ├── genome_2.txt
│   ├── genome_3.txt
│   └── genome_4.txt
└── web_app/                          # Main application
    ├── shiny_ui.py                   # UI and server logic (2,380+ lines)
    ├── utilities.py                  # Data processing functions (357 lines)
    ├── custom_css.py                 # Styling components
    ├── old_intro_consent_card_ui.py  # Legacy UI components
    └── www/                          # Static assets
        ├── blood5.jpg
        ├── sports2.jpg
        ├── gene2.webp
        └── other_images...
```

## Scientific Background

### Polygenic Risk Scoring (PGS)

SNPster implements polygenic risk scoring based on the methodology described in PMID 32714365. The core formula used:

**Population SNP Score:**
```
Population_score_snp = frequency_snp × 2 × beta_snp
```

**Zero-centered Score:**
```
Zero_centered_score = Σ(Beta_snp × Effect_allele_count_snp - Population_score_snp)
```

**Z-score:**
```
Z_score = Zero_centered_score / Standard_deviation_population
```

### Data Processing Pipeline

1. **Variant Annotation**: Cross-reference user SNPs with GWAS catalog
2. **Population Frequency Calculation**: Use 1000 Genomes Project phase 3 data
3. **Risk Score Computation**: Calculate polygenic risk scores for each trait
4. **Statistical Analysis**: Convert scores to percentiles and z-scores
5. **Visualization**: Generate interactive plots and summaries

## Data Sources

- **GWAS Catalog**: European Bioinformatics Institute (accessed November 20, 2024)
- **1000 Genomes Project**: Phase 3 30x dataset (GRCh38)
- **SNPedia**: Community-curated genetic variant database accessed November 21, 2024
- **Athletic Performance SNPs**: Curated from peer-reviewed literature

## Technical Features

### Advanced Capabilities

- **Cross-platform path management** using `pathlib.Path`
- **Environment variable configuration** for flexible deployment
- **Reactive programming** with Shiny framework
- **Statistical computing** with scipy and numpy
- **Interactive visualizations** with Plotly
- **Base64 image encoding** for optimized web delivery
- **Modular architecture** for maintainability

### Code Quality

- **~3,100+ lines** of production Python code
- **Comprehensive error handling** with user feedback
- **Portable path system** for repository cloning
- **Documentation** and setup verification scripts
- **Separation of concerns** between UI and data processing

## Privacy & Ethics

- **Anonymous Analysis**: No personally identifiable information stored
- **Secure Processing**: Data analyzed locally with no external transmission
- **Informed Consent**: Clear privacy terms and user agreement
- **Data Retention**: Raw files stored anonymously for platform improvement
- **Contact**: snpster@gmail.com for questions or concerns

## Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Development Setup

```bash
# Clone the repository
git clone https://github.com/bioinformagician/SNIPster.git
cd SNIPster/genomics_project/travel

# Verify all paths are working
python verify_setup.py

# Run example usage
python example_usage.py
```


## Known Issues

- Some rare variants may not have sufficient population frequency data
- Athletic performance predictions are based on limited genetic markers
- Vendor bias from their custom micro arrays
- Some populations are poorly characterized leading to an inaccurate result

## Future Enhancements

- [ ] Additional ancestry groups and population-specific analysis
- [ ] Integration with SNPster user data to improve score accuracy
- [ ] Implementation of gene imputation to normalize across vendors and improve accuracy
- [ ] Mobile-responsive design improvements
- [ ] Implementation of SnakeMake or Nextflow for pipeline orchestration

## Contact

- **Email**: snpster@gmail.com
- **GitHub**: [bioinformagician](https://github.com/bioinformagician)
- **Project Repository**: [SNIPster](https://github.com/bioinformagician/SNIPster)

## Acknowledgments

- **GWAS Catalog** team at the European Bioinformatics Institute
- **1000 Genomes Project** consortium
- **SNPedia** community contributors
- **Lasse Folkersen** The main inspiration behind this project
---

**Star this repository if you find it useful!**

*SNPster is an educational and research tool. Genetic analysis results should not be used as medical advice. Consult healthcare professionals for medical decisions.*
