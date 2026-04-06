-- Initialize database with custom schema
-- This file will be executed when the PostgreSQL container starts

CREATE SCHEMA IF NOT EXISTS snpster_users;
CREATE SCHEMA IF NOT EXISTS data_libraries;

ALTER DATABASE snpster_db SET search_path = snpster_users, public;

-- ===================================
-- User information table(s) to store user details
-- ===================================

CREATE TABLE snpster_users.user_information (
    user_id varchar(100) PRIMARY KEY,
    email VARCHAR(100) NOT NULL UNIQUE,
    phone_number VARCHAR(20),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    genefile_storage_backend VARCHAR(20) NOT NULL,
    genefile_location TEXT -- stored on linux server in folder tbd
);

-- ===================================
-- Tables for storing information about the PGS scores, mostly containing the information from metadata files
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

CREATE TABLE data_libraries.scoring_files (
    pgs_id varchar(100) PRIMARY KEY REFERENCES data_libraries.pgscatalog_data(pgs_id),
    storage_backend VARCHAR(20) NOT NULL,
    file_path TEXT
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

-- Ontology Trait ID | Ontology Trait Label | Ontology Trait Description | Ontology URL
CREATE TABLE data_libraries.pgs_performance (
    ppm_id varchar(100) PRIMARY KEY,
    pgs_id varchar(100) REFERENCES data_libraries.pgscatalog_data(pgs_id),
    pss_id varchar(100),
    pgp_id varchar(100) REFERENCES data_libraries.pgs_publications(pgp_id),
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
    percent_male_development FLOAT ); 

CREATE TABLE data_libraries.evaluation_sample_sets (
     pss_id VARCHAR(50), 
     individuals_evaluation INTEGER, 
     cases_evaluation INTEGER, 
     controls_evaluation INTEGER, 
     percent_male_evaluation FLOAT );



CREATE TABLE data_libraries.pipeline_dependencies (
    module varchar(100), -- e.g imputer, harmonizer, standardizer dependencies
    dependency_name VARCHAR(100),
    storage_backend VARCHAR(20) NOT NULL,
    file_path TEXT,
    PRIMARY KEY (module, dependency_name, file_path)
);



CREATE TABLE data_libraries.pipeline_dependencies (
    module varchar(100), -- e.g imputer, harmonizer, standardizer dependencies
    dependency_name VARCHAR(100),
    storage_backend VARCHAR(20) NOT NULL,
    file_path TEXT,
    PRIMARY KEY (module, dependency_name, file_path)
);

-- ===================================
-- Reporting tables for post-imputation jobs (PRS calculations, reports)
-- ===================================

CREATE TABLE snpster_users.job_board (
    job_id varchar(100) PRIMARY KEY,
    user_id varchar(100) REFERENCES snpster_users.user_information(user_id),
    job_status VARCHAR(20) NOT NULL CHECK (job_status IN ('queued', 'running', 'completed', 'failed', 'cancelled')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMPTZ,
    finished_at TIMESTAMPTZ
);

CREATE TABLE snpster_users.prsc_job_parameters (
    job_id varchar(100) REFERENCES snpster_users.job_board(job_id) ON DELETE CASCADE,
    pgs_id varchar(100) REFERENCES data_libraries.pgscatalog_data(pgs_id),
    PRIMARY KEY (job_id, pgs_id)
);

CREATE TABLE snpster_users.prsc_job_results (
    job_id varchar(100) PRIMARY KEY REFERENCES snpster_users.job_board(job_id) ON DELETE CASCADE,
    pgs_id varchar(100) REFERENCES data_libraries.pgscatalog_data(pgs_id),
    percentile DECIMAL(5,2) CHECK (percentile >= 0 AND percentile <= 100)
);

CREATE TABLE snpster_users.reports (
    report_id SERIAL PRIMARY KEY,
    job_id varchar(100) REFERENCES snpster_users.job_board(job_id) ON DELETE CASCADE,
    report_status VARCHAR(20) NOT NULL CHECK (report_status IN ('queued', 'running', 'completed', 'failed')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMPTZ,
    finished_at TIMESTAMPTZ,
    report_storage_backend VARCHAR(20) NOT NULL DEFAULT 'local_fs' CHECK (report_storage_backend IN ('local_fs', 's3')),
    report_location TEXT
);

-- ===================================
-- Tables for handling everything from user uploaded data to imputation (including imputation)
-- ===================================

CREATE TABLE snpster_users.imputation_jobs (
    job_id varchar(100) PRIMARY KEY REFERENCES snpster_users.job_board(job_id) ON DELETE CASCADE,
    imputation_status VARCHAR(20) NOT NULL CHECK (imputation_status IN ('queued', 'running', 'completed', 'failed')),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    imputation_started_at TIMESTAMPTZ,
    imputation_finished_at TIMESTAMPTZ,
    imputed_genotype_storage_backend VARCHAR(20) NOT NULL DEFAULT 'local_fs' CHECK (imputed_genotype_storage_backend IN ('local_fs', 's3')),
    imputed_genotype_location TEXT
);

CREATE TABLE snpster_users.imputed_data (
    job_id varchar(100) REFERENCES snpster_users.imputation_jobs(job_id) ON DELETE CASCADE,
    file_type VARCHAR(20), -- e.g, imputed VCF, samplesheet
    storage_backend VARCHAR(20) NOT NULL DEFAULT 'local_fs' CHECK (storage_backend IN ('local_fs', 's3')),
    file_path TEXT,
    PRIMARY KEY (job_id, file_type, file_path)
);
