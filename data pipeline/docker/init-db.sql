-- Initialize database with custom schema
-- This file will be executed when the PostgreSQL container starts

CREATE SCHEMA IF NOT EXISTS snpster_users;
CREATE SCHEMA IF NOT EXISTS data_libraries;

ALTER DATABASE snpster_db SET search_path = snpster_users, public;

CREATE SEQUENCE IF NOT EXISTS snpster_users.imputation_jobs_seq START 1;
CREATE SEQUENCE IF NOT EXISTS snpster_users.prsc_jobs_seq START 1;
CREATE SEQUENCE IF NOT EXISTS snpster_users.report_jobs_seq START 1;

-- ===================================
-- User information table(s) to store user details
-- ===================================

CREATE TABLE snpster_users.user_information (
    user_id varchar(100) PRIMARY KEY NOT NULL,
    email VARCHAR(100) NOT NULL UNIQUE,
    password_hash VARCHAR(255) NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    genefile_location TEXT NOT NULL
);

-- ===================================
-- Tables for storing information about the PGS scores
-- ===================================

CREATE TABLE data_libraries.pgscatalog_data (
    pgs_id varchar(100) PRIMARY KEY,
    pgs_name VARCHAR(255),
    reported_trait VARCHAR(255),
    mapped_trait_efo_label VARCHAR(255),
    efo_id VARCHAR(255),
    pgs_development_method VARCHAR(255),
    pgs_development_details TEXT,
    original_genome_build VARCHAR(20),
    number_of_variants INTEGER,
    number_of_interaction_terms INTEGER,
    type_of_variant_weight TEXT,
    pgp_id varchar(100),
    publication_pmid int,
    publication_doi VARCHAR(255),
    score_and_results_match_original_publication BOOLEAN,
    ancestry_distribution_source_of_variant_associations_gwas VARCHAR(255),
    ancestry_distribution_score_development_training VARCHAR(255),
    ancestry_distribution_pgs_evaluation VARCHAR(255),
    ftp_link VARCHAR(255),
    release_date DATE,
    license_terms_of_use TEXT
);


CREATE TABLE data_libraries.pgs_publications (
    pgp_id varchar(100) PRIMARY KEY,
    first_author VARCHAR(255),
    title VARCHAR(255),
    journal_name VARCHAR(255),
    publication_date DATE,
    release_date DATE,
    authors TEXT,
    digital_object_identifier_doi VARCHAR(255),
    pubmed_id_pmid VARCHAR(20)
);

CREATE TABLE data_libraries.ontology_mappings (
    ontology_id varchar(100) PRIMARY KEY,
    ontology_label VARCHAR(255),
    ontology_description TEXT,
    ontology_url VARCHAR(255)
);

CREATE TABLE data_libraries.pgs_performance (
    ppm_id varchar(100) PRIMARY KEY,
    pgs_id varchar(100) REFERENCES data_libraries.pgscatalog_data(pgs_id) ON DELETE CASCADE,
    pss_id varchar(100),
    pgp_id varchar(100) REFERENCES data_libraries.pgs_publications(pgp_id) ON DELETE CASCADE,
    reported_trait VARCHAR(255),
    covariates_included_in_model TEXT,
    pgs_performance_other_relevant_info TEXT,
    publication_pmid int,
    publication_doi VARCHAR(255),
    hazard_ratio DECIMAL(10,6),
    odds_ratio DECIMAL(10,6),
    beta DECIMAL(10,6),
    auroc DECIMAL(5,2),
    concordance_statistic DECIMAL(5,2),
    other_metric TEXT
);

CREATE TABLE data_libraries.score_development_samples ( 
    pgs_id VARCHAR(50), 
    stage_of_pgs_development VARCHAR(100), 
    individuals_development INTEGER, 
    cases_development INTEGER, 
    controls_development INTEGER, 
    percent_male_development FLOAT
); 

CREATE TABLE data_libraries.evaluation_sample_sets (
    pss_id VARCHAR(50), 
    individuals_evaluation INTEGER, 
    cases_evaluation INTEGER, 
    controls_evaluation INTEGER, 
    percent_male_evaluation FLOAT
);

CREATE TABLE data_libraries.pipeline_dependencies (
    module varchar(100),
    dependency_name VARCHAR(100),
    file_path TEXT,
    PRIMARY KEY (module, dependency_name, file_path)
);

-- ===================================
-- Tables for handling everything from user uploaded data to imputation
-- ===================================

CREATE TABLE snpster_users.imputation_jobs (
    imputation_id integer PRIMARY KEY DEFAULT nextval('snpster_users.imputation_jobs_seq') NOT NULL,
    user_id varchar(100) NOT NULL REFERENCES snpster_users.user_information(user_id) ON DELETE CASCADE,
    imputation_status VARCHAR(20) CHECK (imputation_status IN ('queued', 'running', 'completed', 'failed')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ
);

CREATE TABLE snpster_users.imputed_data (
    imputation_id integer NOT NULL REFERENCES snpster_users.imputation_jobs(imputation_id) ON DELETE CASCADE,
    file_type VARCHAR(20) CHECK (file_type IN ('imputed VCF', 'samplesheet')),
    file_path TEXT,
    PRIMARY KEY (imputation_id, file_type, file_path)
);

-- ===================================
-- Reporting tables for post-imputation jobs
-- ===================================

CREATE TABLE snpster_users.prsc_jobs (
    imputation_id integer NOT NULL REFERENCES snpster_users.imputation_jobs(imputation_id) ON DELETE CASCADE,
    prsc_id integer PRIMARY KEY DEFAULT nextval('snpster_users.prsc_jobs_seq') NOT NULL,
    prsc_status VARCHAR(20) CHECK (prsc_status IN ('queued', 'running', 'completed', 'failed')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ
);

CREATE TABLE snpster_users.prsc_job_parameters (
    prsc_id integer NOT NULL REFERENCES snpster_users.prsc_jobs(prsc_id) ON DELETE CASCADE,
    pgs_id varchar(100) NOT NULL REFERENCES data_libraries.pgscatalog_data(pgs_id) ON DELETE CASCADE
);

CREATE TABLE snpster_users.prsc_job_results (
    prsc_id integer NOT NULL REFERENCES snpster_users.prsc_jobs(prsc_id) ON DELETE CASCADE,
    pgs_id varchar(100) NOT NULL REFERENCES data_libraries.pgscatalog_data(pgs_id) ON DELETE CASCADE,
    percentile DECIMAL(5,2) CHECK (percentile >= 0 AND percentile <= 100),
    z_most_similar_pop DECIMAL(5,2)
);

CREATE TABLE snpster_users.pgs_reports_shop(
    pgs_id varchar(100) NOT NULL REFERENCES data_libraries.pgscatalog_data(pgs_id),
    report_name VARCHAR(255),
    scoring_file_path TEXT
);

CREATE TABLE snpster_users.report_jobs (
    prsc_id integer NOT NULL REFERENCES snpster_users.prsc_jobs(prsc_id) ON DELETE CASCADE,
    report_id integer PRIMARY KEY DEFAULT nextval('snpster_users.report_jobs_seq') NOT NULL,
    report_status VARCHAR(20) CHECK (report_status IN ('queued', 'running', 'completed', 'failed')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ
);

-- ===================================
-- Trigger: create imputation job when a new user is created
-- ===================================

CREATE OR REPLACE FUNCTION snpster_users.create_imputation_job_for_new_user()
RETURNS TRIGGER
LANGUAGE plpgsql
AS $$
BEGIN
    INSERT INTO snpster_users.imputation_jobs (
        user_id,
        imputation_status
    )
    VALUES (
        NEW.user_id,
        'queued'
    );

    RETURN NEW;
END;
$$;

DROP TRIGGER IF EXISTS trg_create_imputation_job_on_user_insert ON snpster_users.user_information;

CREATE TRIGGER trg_create_imputation_job_on_user_insert
AFTER INSERT ON snpster_users.user_information
FOR EACH ROW
EXECUTE FUNCTION snpster_users.create_imputation_job_for_new_user();