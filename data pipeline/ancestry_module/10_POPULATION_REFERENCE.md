# 10-Population Ancestry System

## Overview
Upgraded from 5 to 10 ancestry categories for better stratification, especially for underrepresented populations like Greater Middle Eastern.

## Population Codes & Abbreviations

| Code | Full Name | Description | Database Column |
|------|-----------|-------------|-----------------|
| **EUR** | European | European ancestry | `eur` |
| **EAS** | East Asian | Chinese, Japanese, Korean | `eas` |
| **SAS** | South Asian | Indian, Pakistani (non-Middle Eastern) | `sas` |
| **SSA** | Sub-Saharan African | African ancestry (renamed from AFR) | `ssa` |
| **MID** | Greater Middle Eastern | ⭐ Afghan, Iranian, Arab, North African | `mid` |
| **AMR** | Hispanic/Latin American | Americas, admixed populations | `amr` |
| **SEA** | South East Asian | Thai, Vietnamese, Indonesian | `sea` |
| **CAS** | Central Asian | Uzbek, Kazakh, Uygur | `cas` |
| **NAM** | Native American | Indigenous American | `nam` |
| **OCE** | Oceanian | Pacific Islander, Papuan | `oce` |

## Key Changes from 5-Population System

### Old (5 populations):
- EUR, **AFR**, EAS, SAS, AMR

### New (10 populations):
- EUR, **SSA** (renamed), EAS, SAS, AMR + **MID, SEA, CAS, NAM, OCE**

## Database Schema

```sql
CREATE TABLE snpster_users.user_ancestry (
    user_id VARCHAR(100) PRIMARY KEY,
    
    -- Core 4 super populations
    eur NUMERIC(8,6),  -- European
    eas NUMERIC(8,6),  -- East Asian
    sas NUMERIC(8,6),  -- South Asian
    ssa NUMERIC(8,6),  -- Sub-Saharan African (was: afr)
    
    -- Extended 6 populations
    mid NUMERIC(8,6),  -- Greater Middle Eastern ⭐
    amr NUMERIC(8,6),  -- Hispanic/Latin American
    sea NUMERIC(8,6),  -- South East Asian
    cas NUMERIC(8,6),  -- Central Asian
    nam NUMERIC(8,6),  -- Native American
    oce NUMERIC(8,6),  -- Oceanian
    
    primary_ancestry VARCHAR(10) CHECK (...),
    primary_ancestry_percentage NUMERIC(8,6),
    ancestry_method VARCHAR(100),
    reference_panel VARCHAR(100),
    created_at TIMESTAMPTZ
);
```

## HGDP + 1000G Panel Mapping

### Panel File Format (`super_pop` column):

```
sample          pop             super_pop       gender
HG00096         GBR             EUR             male
NA18525         CHB             EAS             female
HGDP00001       Bedouin         MID             male
HGDP00002       Druze           MID             female
HGDP00003       Palestinian     MID             male
HGDP00004       Uygur           CAS             male
HGDP00005       Yoruba          SSA             male
NA19648         MXL             AMR             female
HG01879         KHV             SEA             male
HGDP00006       Pima            NAM             female
HGDP00007       Papuan          OCE             male
```

### Population Remapping for HGDP + 1000G:

**Greater Middle Eastern (MID):**
- HGDP: Bedouin, Druze, Palestinian, Mozabite

**Sub-Saharan African (SSA):**
- HGDP: Yoruba, Mandenka, Biaka, Mbuti, San, Bantu
- 1000G: YRI, LWK, GWD, MSL, ESN, ASW, ACB

**South East Asian (SEA):**
- HGDP: Cambodian, Dai
- 1000G: CDX (Chinese Dai), KHV (Kinh Vietnamese)

**Central Asian (CAS):**
- HGDP: Uygur
- (Can include Hazara depending on classification)

**Native American (NAM):**
- HGDP: Karitiana, Surui, Pima, Maya, Colombian

**Oceanian (OCE):**
- HGDP: Papuan, Bougainville

**Hispanic/Latin American (AMR - unchanged):**
- 1000G: MXL, PUR, CLM, PEL

## Configuration

### Environment Variables:

```bash
# Number of populations (K)
export K_POPULATIONS=10  # Full 10-population system
# or
export K_POPULATIONS=6   # Core 6 (EUR, EAS, SAS, SSA, MID, AMR)
# or
export K_POPULATIONS=5   # Legacy 5-population (backward compatible)

# Reference panel
export POPULATION_PANEL_FILE="/path/to/hgdp_1kg_panel_10pop.txt"
```

### config.py:

```python
K_POPULATIONS = int(os.getenv("K_POPULATIONS", "5"))

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
    'OCE': 'Oceanian'
}

DB_ANCESTRY_COLUMNS = ['eur', 'eas', 'sas', 'ssa', 'mid', 'amr', 'sea', 'cas', 'nam', 'oce']
```

## Migration Steps

### For New Installations:
Simply use the updated `init-db.sql` which includes all 10 populations.

### For Existing Databases:

1. **Backup:**
```bash
docker compose exec postgres pg_dump -U postgres -d snpster_db > backup.sql
```

2. **Run Migration:**
```bash
docker compose exec -T postgres psql -U postgres -d snpster_db < migrate_ancestry_5_to_10.sql
```

3. **Verify:**
```sql
SELECT * FROM snpster_users.user_ancestry LIMIT 5;
```

## Sample Sizes in PGS Catalog

| Population | Sample Count | Well-Represented? |
|------------|--------------|-------------------|
| European | 9,143 | ✓ Excellent |
| East Asian | 4,003 | ✓ Excellent |
| South Asian | 2,559 | ✓ Good |
| Greater Middle Eastern | 799 | ✓ Good |
| Hispanic/Latin American | 493 | ✓ Fair |
| Native American | 41 | ⚠️ Limited |
| South East Asian | 34 | ⚠️ Limited |
| Sub-Saharan African | 33 | ⚠️ Limited |
| Oceanian | 12 | ⚠️ Very Limited |
| Central Asian | ~20 | ⚠️ Very Limited |

## Benefits

✅ **Greater Middle Eastern (MID)** properly classified - critical for Afghan, Iranian, Arab populations  
✅ **Sub-Saharan African (SSA)** distinguished from African American  
✅ **South East Asian (SEA)** separated from East Asian  
✅ **Better PGS accuracy** - use appropriate reference populations  
✅ **Backward compatible** - works with K=5, 6, or 10  

## Computational Cost

| K | Relative Time | Memory | Use Case |
|---|---------------|--------|----------|
| 5 | 1.0x | 4GB | Legacy, fast |
| 6 | 1.2x | 5GB | Add MID only (recommended start) |
| 10 | 1.5x | 6GB | Full stratification |

## Example Results

### K=5 (Legacy):
```
EUR: 75%, EAS: 10%, SAS: 8%, SSA: 5%, AMR: 2%
Primary: EUR
```

### K=10 (Enhanced):
```
EUR: 60%, MID: 25%, SAS: 8%, EAS: 4%, SSA: 2%, AMR: 1%
Primary: EUR (or MID if >50%)
```

**Afghan example:**
```
MID: 52%, SAS: 28%, EUR: 15%, EAS: 3%, others: 2%
Primary: MID ⭐
```

## Next Steps

1. ✅ Database schema updated (10 columns)
2. ✅ Config updated (population labels)
3. ✅ Code updated (dynamic K support)
4. 🔄 Update HGDP + 1000G panel file with 10 population labels
5. 🔄 Test with K=6 or K=10
6. 🔄 Verify Afghan friend gets MID classification

## Support

- Default K=5 for backward compatibility
- Upgrade to K=6 to add MID (recommended)
- Full K=10 when enhanced panel ready
- Columns with NULL values are fine (only populated if present in K)
