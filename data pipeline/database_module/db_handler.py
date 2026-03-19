from config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH
import psycopg2
from psycopg2 import sql
import pandas as pd

class DbHandler:
    
    def __init__(self, port, db_url, user, password, host, connection=None, cursor=None):
        self.db_url = db_url
        self.user = user
        self.password = password
        self.port = port
        self.host = host
        self.connection = connection
        self.cursor = cursor

    def connect(self):
        # Code to establish a connection to the database using the provided parameters
        try:
            self.connection = psycopg2.connect(
                dbname=DATABASE_NAME,
                user=USERNAME,
                password=PASSWORD,
                host=HOST,
                port=PORT
            )
            
            self.cursor = self.connection.cursor()
            print("Database connection established successfully.")
            return True
        except Exception as e:
            print(f"Error connecting to the database: {e}")
            self.connection = None
            self.cursor = None
            return False
        

    def execute_query(self, query, params=None):
        # Code to execute a given SQL query using the established connection
        try:
            self.cursor.execute(query, params)
            self.connection.commit()
            print("Query executed successfully.")
        except Exception as e:
            print(f"Error executing query: {e}")
            print(f"Query: {query}")
            self.connection.rollback()

    def close(self):
        # Code to close the database connection
        if self.cursor:
            self.cursor.close()
        if self.connection:
            self.connection.close()
        print("Database connection closed.")
        
        
class DbUtils:
    
    def __init__(self, db_handler):
        self.db_handler = db_handler
        
    
    def read_pgs_metadata_excel(self, excel_filepath:str) -> pd.DataFrame: #C:\Users\frezz\Downloads\snpster\data pipeline\reporting_module\data\pgs_all_metadata.xlsx
        
        pgscatalog_data = pd.read_excel(excel_filepath, sheet_name="Scores")
        pgs_publications = pd.read_excel(excel_filepath, sheet_name="Publications")
        ontology_mappings = pd.read_excel(excel_filepath, sheet_name="EFO Traits")
        pgs_performance = pd.read_excel(excel_filepath, sheet_name="Performance Metrics")
        
        pgscatalog_data.rename(columns={
            'Polygenic Score (PGS) ID': 'pgs_id',
            'PGS Name': 'pgs_name',
            'Reported Trait': 'reported_trait',
            'Mapped Trait(s) (EFO label)': 'mapped_trait_efo_label',
            'Mapped Trait(s) (EFO ID)': 'mapped_trait_efo_id',
            'PGS Development Method': 'pgs_development_method',
            'PGS Development Details/Relevant Parameters': 'pgs_development_details',
            'Original Genome Build': 'original_genome_build',
            'Number of Variants': 'number_of_variants',
            'Number of Interaction Terms': 'number_of_interaction_terms',
            'Type of Variant Weight': 'type_of_variant_weight',
            'PGS Publication (PGP) ID': 'pgp_id',
            'Publication (PMID)': 'publication_pmid',
            'Publication (doi)': 'publication_doi',
            'Score and results match the original publication': 'score_and_results_match_original_publication',
            'Ancestry Distribution (%) - Source of Variant Associations (GWAS)': 'ancestry_distribution_source_of_variant_associations_gwas',
            'Ancestry Distribution (%) - Score Development/Training': 'ancestry_distribution_score_development_training',
            'Ancestry Distribution (%) - PGS Evaluation': 'ancestry_distribution_pgs_evaluation',
            'FTP link': 'ftp_link',
            'Release Date': 'release_date',
            'License/Terms of Use': 'license_terms_of_use'
                }, inplace=True)

        pgs_publications.rename(columns={
            'PGS Publication/Study (PGP) ID': 'pgp_id',
            'First Author': 'first_author',
            'Title': 'title',
            'Journal Name': 'journal_name',
            'Publication Date': 'publication_date',
            'Release Date': 'release_date',
            'Authors': 'authors',
            'digital object identifier (doi)': 'digital_object_identifier_doi',
            'PubMed ID (PMID)': 'pubmed_id_pmid'
            }, inplace=True)

        ontology_mappings.rename(columns={
            'Ontology Trait ID': 'ontology_id',
            'Ontology Trait Label': 'ontology_label',
            'Ontology Trait Description': 'ontology_description',
            'Ontology URL': 'ontology_url'
        }, inplace=True)

        pgs_performance.rename(columns={
            'PGS Performance Metric (PPM) ID': 'performance_id',
            'Evaluated Score': 'pgs_id',
            'PGS Sample Set (PSS)': 'pss_id',
            'PGS Publication (PGP) ID': 'pgp_id',
            'Reported Trait': 'reported_trait',
            'Covariates Included in the Model': 'covariates_included_in_model',
            'PGS Performance: Other Relevant Information': 'pgs_performance_other_relevant_info',
            'Publication (PMID)': 'publication_pmid',
            'Publication (doi)': 'publication_doi',
            'Hazard Ratio (HR)': 'hazard_ratio',
            'Odds Ratio (OR)': 'odds_ratio',
            'Beta': 'beta',
            'Area Under the Receiver-Operating Characteristic Curve (AUROC)': 'auroc',
            'Concordance Statistic (C-index)': 'concordance_statistic',
            'Other Metric(s)': 'other_metric'
        }, inplace=True)
        
        return {"pgscatalog_data": pgscatalog_data,
                "pgs_publications": pgs_publications,
                "ontology_mappings": ontology_mappings,
                "pgs_performance": pgs_performance}
    
    def truncate_all_tables(self):
        # Code to truncate all relevant tables in the database before inserting new data
        tables = ["pgscatalog_data", "pgs_publications", "ontology_mappings", "pgs_performance"]
        for table in tables:
            query = sql.SQL("TRUNCATE TABLE {schema}.{table} RESTART IDENTITY CASCADE").format(
                schema=sql.Identifier("data_libraries"),
                table=sql.Identifier(table)
            )
            self.db_handler.execute_query(query)
            print(f"Table {table} truncated successfully.")

        
    
    def insert_dataframe_to_db(self, dataframe: pd.DataFrame, table_name: str):
        # Code to insert the given dataframe into the specified table in the database
        for index, row in dataframe.iterrows():
            columns = list(row.index)
            values = [None if pd.isna(value) else value for value in row.values]

            query = sql.SQL("INSERT INTO {schema}.{table} ({fields}) VALUES ({placeholders})").format(
                schema=sql.Identifier("data_libraries"),
                table=sql.Identifier(table_name),
                fields=sql.SQL(", ").join(sql.Identifier(column) for column in columns),
                placeholders=sql.SQL(", ").join(sql.Placeholder() for _ in columns),
            )
            self.db_handler.execute_query(query, values)
        
        


if __name__ == "__main__":
    db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
    if not db_handler.connect():
        print("Exiting: database connection could not be established.")
        raise SystemExit(1)
    
    
    db_utils = DbUtils(db_handler)
    
    db_utils.truncate_all_tables()
    
    excel_filepath = PGS_EXCEL_FILEPATH
    dataframes = db_utils.read_pgs_metadata_excel(excel_filepath)
    
    # Example of how to access the dataframes
    pgscatalog_data_df = dataframes["pgscatalog_data"]
    pgs_publications_df = dataframes["pgs_publications"]
    ontology_mappings_df = dataframes["ontology_mappings"]
    pgs_performance_df = dataframes["pgs_performance"]
    
    # Print the columns of each dataframe to verify
    print("pgscatalog_data columns:")
    print(pgscatalog_data_df.columns)
    
    print("pgs_publications columns:")
    print(pgs_publications_df.columns)
    
    print("ontology_mappings columns:")
    print(ontology_mappings_df.columns)
    
    print("pgs_performance columns:")
    print(pgs_performance_df.columns)
    
    
    for key, df in dataframes.items():
        print(f"{key} dataframe shape: {df.shape}")
        db_utils.insert_dataframe_to_db(df, table_name=key)
    
    db_handler.close()
    
    