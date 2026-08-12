#!/usr/bin/env python3
"""
Path verification script for the genomics project.
Run this script to verify that all required files and directories are in place.
"""

import sys
from pathlib import Path

# Add the web_app directory to the path so we can import utilities
sys.path.append(str(Path(__file__).parent / "web_app"))

try:
    from utilities import (
        PROJECT_ROOT, 
        GENOMES_DIR, 
        BACKGROUND_DATA_DIR, 
        WWW_DIR, 
        gwas_data_path,
        list_available_genome_files,
        get_genome_file_path,
        get_background_data_path
    )
    print("✅ Successfully imported utilities")
except ImportError as e:
    print(f"❌ Failed to import utilities: {e}")
    sys.exit(1)

def check_path(path, description, required=True):
    """Check if a path exists and report status."""
    if path.exists():
        print(f"✅ {description}: {path}")
        return True
    else:
        status = "❌" if required else "⚠️ "
        print(f"{status} {description}: {path} (NOT FOUND)")
        return False

def main():
    print("🔬 Genomics Project - Path Verification")
    print("=" * 50)
    
    print(f"\nProject Structure:")
    print(f"PROJECT_ROOT: {PROJECT_ROOT}")
    print(f"GENOMES_DIR: {GENOMES_DIR}")
    print(f"BACKGROUND_DATA_DIR: {BACKGROUND_DATA_DIR}")
    print(f"WWW_DIR: {WWW_DIR}")
    
    print(f"\nChecking required directories...")
    all_good = True
    
    # Check directories
    all_good &= check_path(PROJECT_ROOT, "Project root directory")
    all_good &= check_path(GENOMES_DIR, "Genomes directory")
    all_good &= check_path(BACKGROUND_DATA_DIR, "Background data directory")
    all_good &= check_path(WWW_DIR, "Web assets directory")
    
    print(f"\nChecking required files...")
    
    # Check GWAS data
    all_good &= check_path(gwas_data_path, "GWAS data file")
    
    # Check background data
    sports_snps_path = get_background_data_path("sports_snps.txt")
    all_good &= check_path(sports_snps_path, "Sports SNPs file")
    
    # Check web assets
    required_images = [
        "blood5.jpg",
        "sports2.jpg", 
        "gene2.webp",
        "pexels-pixabay-531880.jpg"
    ]
    
    for image in required_images:
        image_path = WWW_DIR / image
        all_good &= check_path(image_path, f"Image: {image}")
    
    print(f"\nGenome files available:")
    genome_files = list_available_genome_files()
    if genome_files:
        for genome_file in genome_files:
            print(f"  📊 {genome_file}")
    else:
        print("  ⚠️  No genome files found")
        all_good = False
    
    print(f"\n" + "=" * 50)
    if all_good:
        print("🎉 All paths verified! Your project is ready to run.")
        print("\nTo start the application:")
        print("  cd web_app")
        print("  python shiny_ui.py")
    else:
        print("❌ Some files or directories are missing.")
        print("Please ensure all required files are in place before running the application.")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
