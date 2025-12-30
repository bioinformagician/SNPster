from data_models import DataContainer, FileHandler, pd, os
from pyliftover import LiftOver

class EnvironmentHandler:
    
    def __init__(self,
                 chain_file_dict: dict[str, str],
                 grch_to_hg_identifier_dict: dict[str, str],
                 user_file:str,
                 output_dir: str
                 ) -> None:
        
        self.output_dir = output_dir
        self.chain_file_dict = chain_file_dict
        self.user_file = user_file
        self.grch_to_hg_identifier_dict = grch_to_hg_identifier_dict
    
    def validate_environment(self) -> None:
        
        paths = [self.output_dir] + list(self.chain_file_dict.values()) + [self.user_file]
        #remove None values (for GRCh38 chain file)
        paths = [path for path in paths if path is not None]
        for path in paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")
    

class WorkflowOrchestrator:
    def __init__(self,
                    data_container: DataContainer,
                    file_handler: FileHandler,
                    environment_handler: EnvironmentHandler
                    ):
        
        self.data_container = data_container
        self.file_handler = file_handler
        self.environment_handler = environment_handler

    def check_dictionary_coherence(self) -> None:
        """Check that the keys and values in the build extraction tool matches the liftover files dicts"""
        for key in self.file_handler.genome_build_dict.values():
            if key not in self.environment_handler.chain_file_dict.keys():
                raise ValueError(f"Genome build '{key}' in genome_build_dict does not have a corresponding chain file in environment_handler.chain_file_dict")

        for key in self.environment_handler.grch_to_hg_identifier_dict.keys():
            if key not in self.file_handler.genome_build_dict.values():
                raise ValueError(f"Genome build '{key}' in grch_to_hg_identifier_dict does not have a corresponding entry in file_handler.genome_build_dict")
        
    def set_vendor(self) -> None:
        self.data_container.vendor = self.file_handler.identify_vendor()
        print(f"Identified vendor: {self.data_container.vendor}")
    
    def set_genome_build(self) -> None:
        self.data_container.genome_build = self.file_handler.identify_genome_build()
        print(f"Identified genome build: {self.data_container.genome_build}")
        
        if self.data_container.genome_build != 'GRCh38':
            self.data_container.lift_over = True
            
    
    def set_microarray_data(self) -> pd.DataFrame:
        self.data_container.microarray_data = self.file_handler.normalize_file()
    
    def evaluate_liftover(self) -> None:
        
        if self.data_container.lift_over:
            print(f"Performing liftover from {self.data_container.genome_build} to GRCh38...")
            lo = LiftOver(self.environment_handler.grch_to_hg_identifier_dict[self.data_container.genome_build], 
                          self.environment_handler.grch_to_hg_identifier_dict['GRCh38'])
            self.data_container.lift_over_data(lo)
    
    def evaluate_zipping(self) -> None:
        if self.file_handler.is_zipped_file():
            unzipped_file = self.file_handler.unzip_file()
            self.file_handler.user_file = unzipped_file
            self.environment_handler.user_file = unzipped_file
    
    def write_parquet_output(self) -> None:
        output_path = os.path.join(self.environment_handler.output_dir, "standardized_microarray_data.parquet")
        self.data_container.microarray_data.to_parquet(output_path)
        print(f"Standardized data written to {output_path}")
            
        
            
        


