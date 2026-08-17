-- Migration script to upgrade user_ancestry table from 5 to 10 populations
-- Upgrades: 5 populations (EUR, AFR, EAS, SAS, AMR) → 10 populations (EUR, EAS, SAS, SSA, MID, AMR, SEA, CAS, NAM, OCE)
-- This script should be run on existing databases BEFORE applying new init-db.sql

-- ===================================
-- Step 1: Backup existing data
-- ===================================

CREATE TABLE IF NOT EXISTS snpster_users.user_ancestry_backup AS 
SELECT * FROM snpster_users.user_ancestry;

-- ===================================
-- Step 2: Update constraints FIRST (before renaming)
-- ===================================

-- Drop old primary_ancestry constraint
ALTER TABLE snpster_users.user_ancestry 
    DROP CONSTRAINT IF EXISTS user_ancestry_primary_ancestry_check;

-- Add new constraint with all 10 populations (+ AFR for temporary backward compatibility during migration)
ALTER TABLE snpster_users.user_ancestry
    ADD CONSTRAINT user_ancestry_primary_ancestry_check 
    CHECK (primary_ancestry IN (
        'EUR', 'EAS', 'SAS', 'SSA', 'MID', 'AMR', 'SEA', 'CAS', 'NAM', 'OCE', 'AFR'
    ));

-- ===================================
-- Step 3: Rename AFR to SSA (Sub-Saharan African)
-- ===================================

ALTER TABLE snpster_users.user_ancestry
    RENAME COLUMN afr TO ssa;

-- ===================================
-- Step 4: Add new ancestry columns
-- ===================================

ALTER TABLE snpster_users.user_ancestry
    ADD COLUMN IF NOT EXISTS mid NUMERIC(8,6) CHECK (mid >= 0 AND mid <= 1),  -- Greater Middle Eastern
    ADD COLUMN IF NOT EXISTS sea NUMERIC(8,6) CHECK (sea >= 0 AND sea <= 1),  -- South East Asian
    ADD COLUMN IF NOT EXISTS cas NUMERIC(8,6) CHECK (cas >= 0 AND cas <= 1),  -- Central Asian
    ADD COLUMN IF NOT EXISTS nam NUMERIC(8,6) CHECK (nam >= 0 AND nam <= 1),  -- Native American
    ADD COLUMN IF NOT EXISTS oce NUMERIC(8,6) CHECK (oce >= 0 AND oce <= 1);  -- Oceanian

-- ===================================
-- Step 5: Update primary_ancestry values
-- ===================================

-- Update AFR to SSA
UPDATE snpster_users.user_ancestry
SET primary_ancestry = 'SSA'
WHERE primary_ancestry = 'AFR';

-- Mark legacy data for potential re-analysis
UPDATE snpster_users.user_ancestry
SET ancestry_method = COALESCE(ancestry_method, 'ADMIXTURE') || '_5POP_LEGACY'
WHERE ancestry_method IS NULL 
   OR (ancestry_method = 'ADMIXTURE' AND mid IS NULL AND sea IS NULL);

-- ===================================
-- Step 6: Verify migration
-- ===================================

DO $$
DECLARE
    original_count INTEGER;
    current_count INTEGER;
    ssa_count INTEGER;
    afr_count INTEGER;
BEGIN
    SELECT COUNT(*) INTO original_count FROM snpster_users.user_ancestry_backup;
    SELECT COUNT(*) INTO current_count FROM snpster_users.user_ancestry;
    SELECT COUNT(*) INTO ssa_count FROM snpster_users.user_ancestry WHERE ssa IS NOT NULL;
    SELECT COUNT(*) INTO afr_count FROM snpster_users.user_ancestry_backup WHERE afr IS NOT NULL;
    
    IF original_count != current_count THEN
        RAISE EXCEPTION 'Migration verification failed: record count mismatch (original: %, current: %)', 
            original_count, current_count;
    END IF;
    
    IF ssa_count != afr_count THEN
        RAISE WARNING 'AFR→SSA migration: original AFR=%, new SSA=%', afr_count, ssa_count;
    END IF;
    
    RAISE NOTICE 'Migration successful: % records migrated', current_count;
    RAISE NOTICE 'AFR column renamed to SSA: % records', ssa_count;
    RAISE NOTICE 'New columns added: MID, SEA, CAS, NAM, OCE (all NULL initially)';
END $$;

-- Show sample of migrated data
SELECT 
    user_id,
    eur, eas, sas, ssa,
    mid, amr, sea, cas, nam, oce,
    primary_ancestry,
    primary_ancestry_percentage,
    ancestry_method
FROM snpster_users.user_ancestry
LIMIT 5;

-- ===================================
-- Summary
-- ===================================

SELECT 
    COUNT(*) as total_users,
    COUNT(eur) as has_eur,
    COUNT(ssa) as has_ssa,
    COUNT(mid) as has_mid,
    STRING_AGG(DISTINCT primary_ancestry, ', ' ORDER BY primary_ancestry) as distinct_primary_ancestries
FROM snpster_users.user_ancestry;

-- ===================================
-- Rollback script (in case of issues)
-- ===================================
-- DROP TABLE IF EXISTS snpster_users.user_ancestry;
-- ALTER TABLE snpster_users.user_ancestry_backup RENAME TO user_ancestry;
