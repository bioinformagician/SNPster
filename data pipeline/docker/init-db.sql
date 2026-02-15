-- Initialize database with custom schema
-- This file will be executed when the PostgreSQL container starts

-- Create your custom schema
CREATE SCHEMA IF NOT EXISTS genomics;

-- Set search path to include your schema
ALTER DATABASE snpster_db SET search_path = genomics, public;

-- ===================================
-- TABLES
-- ===================================

-- SNP (Single Nucleotide Polymorphism) data table
CREATE TABLE genomics.snp_data (
    id SERIAL PRIMARY KEY,
    rsid VARCHAR(20) NOT NULL UNIQUE,
    chromosome INTEGER NOT NULL CHECK (chromosome BETWEEN 1 AND 23),
    position BIGINT NOT NULL,
    ref_allele CHAR(1) NOT NULL,
    alt_allele CHAR(1) NOT NULL,
    gene_symbol VARCHAR(50),
    consequence VARCHAR(100),
    maf DECIMAL(5,4), 
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);


CREATE TABLE genomics.samples (
    id SERIAL PRIMARY KEY,
    sample_id VARCHAR(50) UNIQUE NOT NULL,
    population VARCHAR(10),
    sex CHAR(1) CHECK (sex IN ('M', 'F')),
    age INTEGER CHECK (age > 0 AND age < 150),
    phenotype_status VARCHAR(20),
    batch_id VARCHAR(50),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);


CREATE TABLE genomics.genotypes (
    id SERIAL PRIMARY KEY,
    sample_id INTEGER REFERENCES genomics.samples(id) ON DELETE CASCADE,
    snp_id INTEGER REFERENCES genomics.snp_data(id) ON DELETE CASCADE,
    genotype VARCHAR(3) NOT NULL, -- e.g., "AA", "AT", "TT"
    dosage DECIMAL(3,2) CHECK (dosage >= 0 AND dosage <= 2), -- 0, 1, or 2
    quality_score DECIMAL(5,2),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(sample_id, snp_id)
);


CREATE TABLE genomics.pgs_scores (
    id SERIAL PRIMARY KEY,
    sample_id INTEGER REFERENCES genomics.samples(id) ON DELETE CASCADE,
    trait_name VARCHAR(100) NOT NULL,
    pgs_score DECIMAL(10,6) NOT NULL,
    percentile DECIMAL(5,2) CHECK (percentile >= 0 AND percentile <= 100),
    risk_category VARCHAR(20),
    calculated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);


CREATE TABLE genomics.gwas_results (
    id SERIAL PRIMARY KEY,
    snp_id INTEGER REFERENCES genomics.snp_data(id) ON DELETE CASCADE,
    trait VARCHAR(100) NOT NULL,
    beta DECIMAL(10,6),
    se DECIMAL(10,6), -- Standard error
    p_value SCIENTIFIC NOT NULL,
    odds_ratio DECIMAL(10,6),
    ci_lower DECIMAL(10,6),
    ci_upper DECIMAL(10,6),
    study VARCHAR(200),
    sample_size INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ===================================
-- INDEXES for Performance
-- ===================================

-- SNP data indexes
CREATE INDEX idx_snp_rsid ON genomics.snp_data(rsid);
CREATE INDEX idx_snp_chr_pos ON genomics.snp_data(chromosome, position);
CREATE INDEX idx_snp_gene ON genomics.snp_data(gene_symbol);

-- Sample indexes
CREATE INDEX idx_sample_id ON genomics.samples(sample_id);
CREATE INDEX idx_sample_population ON genomics.samples(population);

-- Genotype indexes (most important for performance)
CREATE INDEX idx_genotype_sample ON genomics.genotypes(sample_id);
CREATE INDEX idx_genotype_snp ON genomics.genotypes(snp_id);
CREATE INDEX idx_genotype_composite ON genomics.genotypes(sample_id, snp_id);

-- PGS indexes
CREATE INDEX idx_pgs_sample ON genomics.pgs_scores(sample_id);
CREATE INDEX idx_pgs_trait ON genomics.pgs_scores(trait_name);

-- GWAS indexes
CREATE INDEX idx_gwas_snp ON genomics.gwas_results(snp_id);
CREATE INDEX idx_gwas_trait ON genomics.gwas_results(trait);
CREATE INDEX idx_gwas_pvalue ON genomics.gwas_results(p_value);

-- ===================================
-- FUNCTIONS & TRIGGERS
-- ===================================

-- Function to update the updated_at timestamp
CREATE OR REPLACE FUNCTION genomics.update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ language 'plpgsql';

-- Trigger for SNP data
CREATE TRIGGER update_snp_data_updated_at 
    BEFORE UPDATE ON genomics.snp_data 
    FOR EACH ROW EXECUTE FUNCTION genomics.update_updated_at_column();

-- ===================================
-- VIEWS for Common Queries
-- ===================================

-- View for complete genotype information
CREATE VIEW genomics.complete_genotypes AS
SELECT 
    s.sample_id,
    s.population,
    sd.rsid,
    sd.chromosome,
    sd.position,
    sd.gene_symbol,
    g.genotype,
    g.dosage,
    g.quality_score
FROM genomics.genotypes g
JOIN genomics.samples s ON g.sample_id = s.id
JOIN genomics.snp_data sd ON g.snp_id = sd.id;

-- View for high-impact variants
CREATE VIEW genomics.high_impact_variants AS
SELECT 
    rsid,
    chromosome,
    position,
    gene_symbol,
    consequence,
    maf
FROM genomics.snp_data
WHERE consequence IN ('stop_gained', 'frameshift_variant', 'start_lost')
   OR maf < 0.01;

-- ===================================
-- PERMISSIONS
-- ===================================

-- Grant permissions to postgres user
GRANT ALL PRIVILEGES ON SCHEMA genomics TO postgres;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA genomics TO postgres;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA genomics TO postgres;
GRANT EXECUTE ON ALL FUNCTIONS IN SCHEMA genomics TO postgres;

-- ===================================
-- SAMPLE DATA (Optional - for testing)
-- ===================================

-- Insert some test SNPs
INSERT INTO genomics.snp_data (rsid, chromosome, position, ref_allele, alt_allele, gene_symbol, maf) VALUES
('rs123456', 1, 1000000, 'A', 'G', 'GENE1', 0.25),
('rs789012', 2, 2000000, 'C', 'T', 'GENE2', 0.15),
('rs345678', 3, 3000000, 'G', 'A', 'GENE3', 0.35);

-- Insert test samples
INSERT INTO genomics.samples (sample_id, population, sex, age) VALUES
('SAMPLE001', 'EUR', 'M', 35),
('SAMPLE002', 'AFR', 'F', 42),
('SAMPLE003', 'ASN', 'M', 28);