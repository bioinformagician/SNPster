#!/bin/bash
# Update HGDP + 1000G panel file from 5 to 10 populations
# Maps specific populations to the new 10-category system

PANEL_FILE="/srv/dependencies/imputation_runner/hgdp_1kg_panel.txt"
OUTPUT_DIR="/home/frederik/github_projects/SNPster/data pipeline/ancestry_module"
BACKUP_FILE="${OUTPUT_DIR}/hgdp_1kg_panel.5pop_backup.txt"
OUTPUT_FILE="${OUTPUT_DIR}/hgdp_1kg_panel_10pop.txt"

# Backup original
cp "$PANEL_FILE" "$BACKUP_FILE"
echo "✓ Backed up original to $BACKUP_FILE"

# Create updated panel with 10 populations
awk -F'\t' 'BEGIN {OFS="\t"}
NR==1 {
    # Print header unchanged
    print $0
    next
}
{
    sample = $1
    pop = $2
    super_pop = $3
    gender = $4
    
    # Remap populations to 10-category system
    
    # MID - Greater Middle Eastern
    if (pop == "Bedouin" || pop == "Druze" || pop == "Mozabite" || pop == "Palestinian") {
        super_pop = "MID"
    }
    # SEA - South East Asian
    else if (pop == "Cambodian" || pop == "Dai" || pop == "CDX" || pop == "KHV") {
        super_pop = "SEA"
    }
    # CAS - Central Asian
    else if (pop == "Uygur") {
        super_pop = "CAS"
    }
    # NAM - Native American (indigenous)
    else if (pop == "Karitiana" || pop == "Maya" || pop == "Pima" || pop == "Surui") {
        super_pop = "NAM"
    }
    # SSA - Sub-Saharan African (all AFR except Middle Eastern)
    else if (super_pop == "AFR") {
        super_pop = "SSA"
    }
    # AMR stays as-is (admixed Hispanic/Latin American)
    # EUR stays as-is (European)
    # EAS stays as-is (East Asian - Chinese, Japanese, Korean, Mongolian)
    # SAS stays as-is (South Asian - Indian, Pakistani, Bangladeshi)
    
    print sample, pop, super_pop, gender
}' "$PANEL_FILE" > "$OUTPUT_FILE"

echo "✓ Created updated panel: $OUTPUT_FILE"

# Show summary of changes
echo ""
echo "Population distribution in updated panel:"
cut -f3 "$OUTPUT_FILE" | tail -n +2 | sort | uniq -c | sort -rn

echo ""
echo "Specific remappings:"
echo "  Bedouin, Druze, Mozabite, Palestinian → MID (Greater Middle Eastern)"
echo "  Cambodian, Dai, CDX, KHV → SEA (South East Asian)"
echo "  Uygur → CAS (Central Asian)"
echo "  Karitiana, Maya, Pima, Surui → NAM (Native American)"
echo "  All AFR → SSA (Sub-Saharan African)"
echo "  AMR, EUR, EAS, SAS unchanged"

echo ""
echo "To activate the new panel, either:"
echo "  1. Copy to original location (requires sudo):"
echo "     sudo cp $OUTPUT_FILE $PANEL_FILE"
echo ""
echo "  2. Update docker-compose.yaml to mount this file"
echo ""
echo "  3. Set environment variable in config:"
echo "     export POPULATION_PANEL_FILE=\"$OUTPUT_FILE\""

