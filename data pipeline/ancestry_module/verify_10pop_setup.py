#!/usr/bin/env python3
"""
Verification script for 10-population ancestry setup.
Checks that all components are properly configured.
"""

import sys
import os

# Add current directory to path to import config
sys.path.insert(0, os.path.dirname(__file__))

from config import (
    K_POPULATIONS, 
    ANCESTRY_LABELS, 
    LEGACY_POPULATION_MAP,
    DB_ANCESTRY_COLUMNS,
    POPULATION_PANEL_FILE
)

def verify_configuration():
    """Verify all configuration is correct for 10-population system."""
    
    print("=" * 70)
    print("10-POPULATION ANCESTRY SYSTEM VERIFICATION")
    print("=" * 70)
    
    issues = []
    warnings = []
    
    # Check 1: K_POPULATIONS setting
    print(f"\n1. K_POPULATIONS: {K_POPULATIONS}")
    if K_POPULATIONS == 5:
        warnings.append("K_POPULATIONS is set to 5 (backward compatibility mode)")
        print("   ⚠️  Currently set to 5 for backward compatibility")
        print("   💡 Set K_POPULATIONS=10 environment variable to use full system")
    elif K_POPULATIONS == 10:
        print("   ✅ Configured for 10 populations")
    else:
        warnings.append(f"K_POPULATIONS is {K_POPULATIONS} (expected 5 or 10)")
    
    # Check 2: ANCESTRY_LABELS
    print(f"\n2. ANCESTRY_LABELS: {len(ANCESTRY_LABELS)} labels defined")
    expected_labels = ['EUR', 'EAS', 'SAS', 'SSA', 'MID', 'AMR', 'SEA', 'CAS', 'NAM', 'OCE', 'AFR']
    missing_labels = set(expected_labels) - set(ANCESTRY_LABELS.keys())
    if missing_labels:
        issues.append(f"Missing ancestry labels: {missing_labels}")
        print(f"   ❌ Missing: {missing_labels}")
    else:
        print("   ✅ All 10 populations + legacy AFR defined")
        for code, label in sorted(ANCESTRY_LABELS.items()):
            marker = "   " if code != 'AFR' else "   (legacy)"
            print(f"      {code}: {label}{' [LEGACY]' if code == 'AFR' else ''}")
    
    # Check 3: DB_ANCESTRY_COLUMNS
    print(f"\n3. DB_ANCESTRY_COLUMNS: {len(DB_ANCESTRY_COLUMNS)} columns")
    expected_cols = ['eur', 'eas', 'sas', 'ssa', 'mid', 'amr', 'sea', 'cas', 'nam', 'oce']
    if DB_ANCESTRY_COLUMNS == expected_cols:
        print("   ✅ Correctly defined for 10 populations")
        print(f"      {', '.join(DB_ANCESTRY_COLUMNS)}")
    else:
        issues.append("DB_ANCESTRY_COLUMNS doesn't match expected 10 populations")
        print(f"   ❌ Expected: {expected_cols}")
        print(f"   ❌ Got: {DB_ANCESTRY_COLUMNS}")
    
    # Check 4: LEGACY_POPULATION_MAP
    print(f"\n4. LEGACY_POPULATION_MAP: {len(LEGACY_POPULATION_MAP)} mappings")
    if 'AFR' in LEGACY_POPULATION_MAP and LEGACY_POPULATION_MAP['AFR'] == 'SSA':
        print("   ✅ AFR → SSA mapping configured")
    else:
        warnings.append("AFR → SSA legacy mapping not found")
        print("   ⚠️  AFR → SSA mapping not configured")
    
    # Check 5: Panel file
    print(f"\n5. POPULATION_PANEL_FILE: {POPULATION_PANEL_FILE}")
    if os.path.exists(POPULATION_PANEL_FILE):
        print("   ✅ Panel file exists")
        
        # Quick check of panel file content
        try:
            with open(POPULATION_PANEL_FILE, 'r') as f:
                header = f.readline().strip().split('\t')
                if 'super_pop' in header:
                    # Read first 1000 lines to check populations
                    pops = set()
                    for i, line in enumerate(f):
                        if i >= 1000:
                            break
                        cols = line.strip().split('\t')
                        if len(cols) > 2:
                            pops.add(cols[2])
                    
                    print(f"      Populations found: {sorted(pops)}")
                    
                    # Check if 10-pop system
                    if 'MID' in pops or 'SSA' in pops:
                        print("   ✅ Panel appears to use 10-population system")
                    elif 'AFR' in pops and 'MID' not in pops:
                        warnings.append("Panel file still uses 5-population system (AFR instead of SSA/MID)")
                        print("   ⚠️  Panel appears to use old 5-population system")
                        print("   💡 Run: bash update_panel_to_10pop.sh")
        except Exception as e:
            warnings.append(f"Could not parse panel file: {e}")
            print(f"   ⚠️  Could not parse: {e}")
    else:
        issues.append(f"Panel file not found: {POPULATION_PANEL_FILE}")
        print(f"   ❌ File not found")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    if issues:
        print(f"\n❌ ISSUES FOUND ({len(issues)}):")
        for i, issue in enumerate(issues, 1):
            print(f"   {i}. {issue}")
    
    if warnings:
        print(f"\n⚠️  WARNINGS ({len(warnings)}):")
        for i, warning in enumerate(warnings, 1):
            print(f"   {i}. {warning}")
    
    if not issues and not warnings:
        print("\n✅ ALL CHECKS PASSED - System fully configured for 10 populations!")
    elif not issues:
        print("\n✅ No critical issues, but review warnings above")
    else:
        print("\n❌ Please fix issues before using 10-population system")
        return False
    
    print("\n" + "=" * 70)
    print("ACTIVATION CHECKLIST")
    print("=" * 70)
    print("""
To activate the 10-population system:

1. ✓ Database schema updated (init-db.sql)
2. ✓ Panel file updated to 10 populations
3. ✓ Config files updated (config.py, ancestry_classes.py)
4. ⚠  Set K_POPULATIONS environment variable:
   
   For docker-compose.yaml, add under ancestry service:
   
   environment:
     - K_POPULATIONS=10
   
   Or for K=6 (core populations only):
   
   environment:
     - K_POPULATIONS=6

5. ⚠  Migrate existing database (if you have user data):
   
   psql -U postgres -d snpster_db < docker/migrate_ancestry_5_to_10.sql

6. Test with sample data to verify MID population detection
""")
    
    return True

if __name__ == "__main__":
    success = verify_configuration()
    sys.exit(0 if success else 1)
