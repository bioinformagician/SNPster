from data_models import DataContainer, FileHandler, pd, os
from pyliftover import LiftOver
import pyarrow as pa
import pyarrow.parquet as pq

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
        
        for vendor in self.file_handler.forward_strand_vendors:
            if vendor not in self.file_handler.accepted_vendors_dict.keys():
                raise ValueError(f"Vendor '{vendor}' in forward_strand_vendors is not in accepted_vendors_dict")
        
    def set_vendor(self) -> None:
        self.data_container.vendor = self.file_handler.identify_vendor()
        print(f"Identified vendor: {self.data_container.vendor}")
    
    def set_genome_build(self) -> None:
        self.data_container.genome_build = self.file_handler.identify_genome_build()
        print(f"Identified genome build: {self.data_container.genome_build}")
        
        if self.data_container.genome_build != 'GRCh38':
            self.data_container.lift_over = True
    
    def set_strand_direction(self) -> None:
        if self.data_container.vendor is None:
            raise ValueError("Vendor must be set before determining strand direction.")
        self.data_container.is_forward_strand = self.data_container.vendor in self.file_handler.forward_strand_vendors
            
    
    def set_microarray_data(self) -> pd.DataFrame:
        self.data_container.microarray_data = self.file_handler.normalize_file()
    
    def evaluate_liftover(self) -> None:
        
        if self.data_container.lift_over:
            print(f"Performing liftover from {self.data_container.genome_build} to GRCh38...")
            lo = LiftOver(self.environment_handler.grch_to_hg_identifier_dict[self.data_container.genome_build], 
                          self.environment_handler.grch_to_hg_identifier_dict['GRCh38'])
            self.data_container.lift_over_data(lo)
            # Update genome_build metadata to reflect the liftover
            self.data_container.genome_build = 'GRCh38'
            print(f"Liftover complete. Genome build updated to GRCh38")
    
    def evaluate_zipping(self) -> None:
        if self.file_handler.is_zipped_file():
            unzipped_file = self.file_handler.unzip_file()
            self.file_handler.user_file = unzipped_file
            self.environment_handler.user_file = unzipped_file
    
    def write_parquet_output(self) -> None:
        output_path = os.path.join(self.environment_handler.output_dir, f"standardized_microarray_data_user_{self.data_container.identifier}.parquet")
        
        # Convert DataFrame to PyArrow Table
        table = pa.Table.from_pandas(self.data_container.microarray_data)
        
        # Add custom metadata
        metadata = {
            b'vendor': self.data_container.vendor.encode('utf-8'),
            b'genome_build': self.data_container.genome_build.encode('utf-8'),
            b'is_forward_strand': str(self.data_container.is_forward_strand).encode('utf-8'),
            b'identifier': str(self.data_container.identifier).encode('utf-8')
        }
        
        # Merge with existing schema metadata
        existing_metadata = table.schema.metadata or {}
        combined_metadata = {**existing_metadata, **metadata}
        
        # Create new schema with metadata
        new_schema = table.schema.with_metadata(combined_metadata)
        table = table.cast(new_schema)
        
        # Write parquet with metadata
        pq.write_table(table, output_path)
        print(f"Standardized data written to {output_path}")
        print(f"Metadata: vendor={self.data_container.vendor}, genome_build={self.data_container.genome_build}, is_forward_strand={self.data_container.is_forward_strand}, identifier={self.data_container.identifier}")
    
    
    def write_meta_data_output(self) -> None:
        meta_data_path = os.path.join(self.environment_handler.output_dir, f"metadata_user_{self.data_container.identifier}.txt")
        with open(meta_data_path, 'w') as f:
            f.write(f"Vendor: {self.data_container.vendor}\n")
            f.write(f"Genome Build: {self.data_container.genome_build}\n")
            f.write(f"Is Forward Strand: {self.data_container.is_forward_strand}\n")
            f.write(f"Identifier: {self.data_container.identifier}\n")
        print(f"Metadata written to {meta_data_path}")
        
            
        


