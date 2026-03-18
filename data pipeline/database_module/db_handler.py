from config import USERNAME, PASSWORD, DATABASE_NAME, HOST, PORT
import psycopg2

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
            print("Database connection established successfully.")
        except Exception as e:
            print(f"Error connecting to the database: {e}")
        
        self.cursor = self.connection.cursor()
        

    def execute_query(self, query):
        # Code to execute a given SQL query using the established connection
        try:
            self.cursor.execute(query)
            self.connection.commit()
            print("Query executed successfully.")
        except Exception as e:
            print(f"Error executing query: {e}")
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
        
    
    @staticmethod
    def initialize_pgs_catalog_tables(self, excel_file_path): #C:\Users\frezz\Downloads\snpster\data pipeline\reporting_module\data\pgs_all_metadata.xlsx
        pass
        