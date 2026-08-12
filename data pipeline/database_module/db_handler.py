from pathlib import Path
from db_config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT, PGS_EXCEL_FILEPATH
import psycopg2
from psycopg2 import sql
from psycopg2.extras import execute_values
import pandas as pd
import time

class DbHandler:
    
    def __init__(self, port:int, user:str, password:str, host:str, connection=None, cursor=None, db_url:str=None):
        self.db_url = db_url
        self.user = user
        self.password = password
        self.port = port
        self.host = host
        self.connection = connection
        self.cursor = cursor

    def __enter__(self): #enable usage of with
        if not self.connect():
            raise RuntimeError("Failed to establish database connection.")
        return self

    def __exit__(self, exc_type, exc, tb):#enable usage of with
        self.close()
        return False

    def connect(self, retries = 10, wait_time = 60) -> bool:
        # Code to establish a connection to the database using the provided parameters
        for attempt in range(retries):
            try:
                self.connection = psycopg2.connect(
                    dbname=self.db_url or DATABASE_NAME,
                    user=self.user,
                    password=self.password,
                    host=self.host,
                port=self.port
            )
            
                self.cursor = self.connection.cursor()
                print("Database connection established successfully.")
                return True
            except Exception as e:
                print(f"Error connecting to the database: {e}")
                self.connection = None
                self.cursor = None
                if attempt < retries - 1:
                    time.sleep(wait_time)  # Wait before retrying
        return False
        

    def close(self) -> None:
        # Code to close the database connection
        if self.cursor:
            self.cursor.close()
        if self.connection:
            self.connection.close()
        self.cursor = None
        self.connection = None
        print("Database connection closed.")

    def _ensure_connection(self) -> bool:
        connection_closed = self.connection is None or self.connection.closed != 0
        cursor_closed = self.cursor is None or self.cursor.closed

        if connection_closed or cursor_closed:
            return self.connect(retries=1, wait_time=0)

        return True
        
    def execute_query(self, query, params=None, retries = 10, wait_time = 60) -> list:
        # Code to execute a given SQL query using the established connection
        
        for attempt in range(retries):

            try:
                if not self._ensure_connection():
                    raise psycopg2.InterfaceError("Database connection is not available.")

                self.cursor.execute(query, params)
                self.connection.commit()
                print("Query executed successfully.")
                
                return self.cursor.fetchall() if self.cursor.description else None
                
            except Exception as e:
                print(f"Error executing query: {e}")
                print(f"Query: {query}")
                if self.connection is not None:
                    try:
                        if self.connection.closed == 0:
                            self.connection.rollback()
                    except Exception:
                        pass
                if isinstance(e, (psycopg2.InterfaceError, psycopg2.OperationalError)):
                    self.close()
                if attempt < retries - 1:
                    time.sleep(wait_time)  # Wait before retrying
        return None
        
        
class DbUtils:
    
    def __init__(self, db_handler:DbHandler):
        self.db_handler = db_handler
        
    
    def get_pd_dataframe_from_query(self, query:str) -> pd.DataFrame:
        # Code to execute a SELECT query and return the results as a pandas DataFrame
        query_result = self.db_handler.execute_query(query)
        if query_result is not None:
            column_names = [desc[0] for desc in self.db_handler.cursor.description]
            df = pd.DataFrame(query_result, columns=column_names)
            return df
        else:
            print("No results returned from the query.")
            return pd.DataFrame()  # Return empty DataFrame if no results
        
        

    
    def read_pgs_metadata_excel(self, excel_filepath:str) -> pd.DataFrame: #C:\Users\frezz\Downloads\snpster\data pipeline\reporting_module\data\pgs_all_metadata.xlsx
        
        pgscatalog_data = pd.read_excel(excel_filepath, sheet_name="Scores")
        pgs_publications = pd.read_excel(excel_filepath, sheet_name="Publications")
        ontology_mappings = pd.read_excel(excel_filepath, sheet_name="EFO Traits")
        pgs_performance = pd.read_excel(excel_filepath, sheet_name="Performance Metrics")
        score_development_samples = pd.read_excel(excel_filepath, sheet_name="Score Development Samples")
        evaluation_sample_sets = pd.read_excel(excel_filepath, sheet_name="Evaluation Sample Sets")
        
        print(score_development_samples.columns)
        print(evaluation_sample_sets.columns)
        
        def _pick_existing_column(df: pd.DataFrame, candidates: list[str], default=None) -> pd.Series:
            cols = list(df.columns)
            normalized = {" ".join(str(c).strip().lower().split()): c for c in cols}

            for candidate in candidates:
                if candidate in cols:
                    return df[candidate]

                key = " ".join(str(candidate).strip().lower().split())
                if key in normalized:
                    return df[normalized[key]]

            return pd.Series([default] * len(df), index=df.index)

        score_development_samples = score_development_samples[["Polygenic Score (PGS) ID", "Stage of PGS Development",
                                                              "Number of Individuals", "Number of Cases",
                                                              "Number of Controls", "Percent of Participants Who are Male"]]

        evaluation_sample_sets = pd.DataFrame({
            "PGS Sample Set (PSS)": _pick_existing_column(evaluation_sample_sets, ["PGS Sample Set (PSS)"]),
            "Polygenic Score (PGS) ID": _pick_existing_column(evaluation_sample_sets, ["Polygenic Score (PGS) ID", "Evaluated Score", "PGS ID"]),
            "Number of Individuals": _pick_existing_column(evaluation_sample_sets, ["Number of Individuals"]),
            "Number of Cases": _pick_existing_column(evaluation_sample_sets, ["Number of Cases"]),
            "Number of Controls": _pick_existing_column(evaluation_sample_sets, ["Number of Controls"]),
            "Percent of Participants Who are Male": _pick_existing_column(evaluation_sample_sets, ["Percent of Participants Who are Male"]),
            "Broad Ancestry Category": _pick_existing_column(evaluation_sample_sets, ["Broad Ancestry Category", "Broad Ancestry Category(ies)", "Broad ancestry category"]),
            "Country of Recruitment": _pick_existing_column(evaluation_sample_sets, ["Country of Recruitment", "Country/Countries of Recruitment", "Countries of recruitment"]),
            "Cohort": _pick_existing_column(evaluation_sample_sets, ["Cohort", "Cohort(s)", "Cohorts", "Cohort(s) Included", "Cohort Name"]),
        })
        
        score_development_samples.rename(columns={
            'Polygenic Score (PGS) ID': 'pgs_id',
            'Stage of PGS Development': 'stage_of_pgs_development',
            'Number of Individuals': 'individuals_development',
            'Number of Cases': 'cases_development',
            'Number of Controls': 'controls_development',
            'Percent of Participants Who are Male': 'percent_male_development'
        }, inplace=True)
        
        evaluation_sample_sets.rename(columns={
            'PGS Sample Set (PSS)': 'pss_id',
            'Polygenic Score (PGS) ID': 'pgs_id',
            'Number of Individuals': 'individuals_evaluation',
            'Number of Cases': 'cases_evaluation',
            'Number of Controls': 'controls_evaluation',
            'Percent of Participants Who are Male': 'percent_male_evaluation',
            'Broad Ancestry Category': 'broad_ancestry_category',
            'Country of Recruitment': 'country_of_recruitment',
            'Cohort': 'cohort',
        }, inplace=True)
        
        pgscatalog_data.rename(columns={
            'Polygenic Score (PGS) ID': 'pgs_id',
            'PGS Name': 'pgs_name',
            'Reported Trait': 'reported_trait',
            'Mapped Trait(s) (EFO label)': 'mapped_trait_efo_label',
            'Mapped Trait(s) (EFO ID)': 'efo_id',
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
            'PGS Performance Metric (PPM) ID': 'ppm_id',
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
        
        # Clean numeric columns: extract first number before confidence intervals, store NULL if no number found
        numeric_cols = ['hazard_ratio', 'odds_ratio', 'beta', 'auroc', 'concordance_statistic']
        for col in numeric_cols:
            if col in pgs_performance.columns:
                # Extract first signed number (e.g., "-0.7" from "-0.7 (0.15)").
                pgs_performance[col] = pgs_performance[col].astype(str).str.extract(
                    r'^\s*([+-]?(?:\d+(?:\.\d+)?|\.\d+))',
                    expand=False,
                )
                # Convert to numeric, NaN for invalid values
                pgs_performance[col] = pd.to_numeric(pgs_performance[col], errors='coerce')
        
        return {"pgscatalog_data": pgscatalog_data,
                "pgs_publications": pgs_publications,
                "ontology_mappings": ontology_mappings,
                "pgs_performance": pgs_performance,
                "score_development_samples": score_development_samples,
                "evaluation_sample_sets": evaluation_sample_sets
                }
    
    def truncate_all_tables(self) -> None:
        # Code to truncate all relevant tables in the database before inserting new data
        tables = ["pgscatalog_data", "pgs_publications", "ontology_mappings", "pgs_performance", "score_development_samples", "evaluation_sample_sets"]
        for table in tables:
            query = sql.SQL("TRUNCATE TABLE {schema}.{table} RESTART IDENTITY CASCADE").format(
                schema=sql.Identifier("data_libraries"),
                table=sql.Identifier(table)
            )
            self.db_handler.execute_query(query)
            print(f"Table {table} truncated successfully.")

        
    
    def insert_dataframe_to_db(self, dataframe: pd.DataFrame, table_name: str, schema: str = "data_libraries") -> None:
        # Bulk insert all dataframe rows in a single efficient batch operation
        if dataframe.empty:
            print(f"No rows to insert for {table_name}")
            return
        
        columns = list(dataframe.columns)
        # Convert dataframe rows to tuples, handling NaN/NaT as None
        values = [
            tuple(None if pd.isna(val) else val for val in row)
            for row in dataframe.itertuples(index=False, name=None)
        ]
        
        # Build INSERT statement with execute_values format
        query = sql.SQL(
            "INSERT INTO {schema}.{table} ({fields}) VALUES %s"
        ).format(
            schema=sql.Identifier(schema),
            table=sql.Identifier(table_name),
            fields=sql.SQL(",").join(sql.Identifier(col) for col in columns),
        )
        
        try:
            execute_values(
                self.db_handler.cursor,
                query.as_string(self.db_handler.connection),
                values,
                page_size=1000,
            )
            self.db_handler.connection.commit()
            print(f"Inserted {len(values)} rows into {schema}.{table_name}")
        except Exception as e:
            print(f"Error inserting into {schema}.{table_name}: {e}")
            self.db_handler.connection.rollback()
            raise

    def upsert_dataframe_to_db(self, dataframe: pd.DataFrame, table_name: str, schema: str, conflict_columns: list) -> None:
        """
        Bulk upsert dataframe rows into database table.
        Updates existing rows on conflict, inserts new ones otherwise.
        
        Args:
            dataframe: DataFrame to upsert
            table_name: Target table name
            schema: Target schema (e.g., "snpster_users")
            conflict_columns: List of columns to use for conflict detection (e.g., ['user_id'])
        """
        if dataframe.empty:
            print(f"No rows to upsert for {table_name}")
            return
        
        columns = list(dataframe.columns)
        # Convert dataframe rows to tuples, handling NaN/NaT as None
        values = [
            tuple(None if pd.isna(val) else val for val in row)
            for row in dataframe.itertuples(index=False, name=None)
        ]
        
        # Build UPSERT query (INSERT ... ON CONFLICT DO UPDATE)
        update_cols = [col for col in columns if col not in conflict_columns]
        query = sql.SQL(
            "INSERT INTO {schema}.{table} ({fields}) VALUES %s "
            "ON CONFLICT ({conflict_fields}) DO UPDATE SET {updates}"
        ).format(
            schema=sql.Identifier(schema),
            table=sql.Identifier(table_name),
            fields=sql.SQL(",").join(sql.Identifier(col) for col in columns),
            conflict_fields=sql.SQL(",").join(sql.Identifier(col) for col in conflict_columns),
            updates=sql.SQL(",").join(
                sql.SQL("{col} = EXCLUDED.{col}").format(col=sql.Identifier(col))
                for col in update_cols
            ),
        )
        
        try:
            execute_values(
                self.db_handler.cursor,
                query.as_string(self.db_handler.connection),
                values,
                page_size=1000,
            )
            self.db_handler.connection.commit()
            print(f"Upserted {len(values)} rows into {schema}.{table_name}")
        except Exception as e:
            print(f"Error upserting into {schema}.{table_name}: {e}")
            self.db_handler.connection.rollback()
            raise
            self.db_handler.connection.rollback()
    
    def get_user_id_from_imputation_id(self, imputation_id:list) -> pd.DataFrame:
        """Input a list of imputation_ids and return a pandas df with the columns imputation_id and user_id. One imputation_id can have multiple user_ids"""
        
        query = sql.SQL("SELECT user_id, imputation_id FROM snpster_users.imputation_jobs WHERE imputation_id = ANY(%s)")
        result = self.db_handler.execute_query(query, (imputation_id,))
        if result:
            return pd.DataFrame(result, columns=['user_id', 'imputation_id'])
        else:
            print(f"No user_id found for imputation_id: {imputation_id}")
            return None
    
        
        


if __name__ == "__main__":
    db_handler = DbHandler(port=PORT, db_url=None, user=USERNAME, password=PASSWORD, host=HOST)
    if not db_handler.connect():
        print("Exiting: database connection could not be established.")
        raise SystemExit(1)
    
    
    db_utils = DbUtils(db_handler)
    
    db_utils.truncate_all_tables()
    
    excel_filepath = PGS_EXCEL_FILEPATH
    dataframes = db_utils.read_pgs_metadata_excel(excel_filepath)
    
    
    for key, df in dataframes.items():
        print(f"{key} dataframe shape: {df.shape}")
        db_utils.insert_dataframe_to_db(df, table_name=key)
    

    #insert 
    
    db_handler.close()
    
    