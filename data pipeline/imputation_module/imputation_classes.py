@dataclass
class DataContainer:
    microarray_data: Optional[pd.DataFrame] = None
    reference_data: Optional[pd.DataFrame] = None
    harmonized_data: Optional[pd.DataFrame] = None
    harmonization_stats: Optional[pd.DataFrame] = None
    
    

class EnvironmentHandler:
    def __init__(self,
                 working_dir: str,
                 user_upload_file: str,
                 plink_1_9_path: str,
                 plink_2_0_path: str,
                 plink_map_file: str,
                 pvar_ref_file: str,
                 PLINK_PREFIX: str,
                 plink_reference_fasta: str,
                 beagle_references: str,
                 chromosome_split_files: dict = None,
                 vcf_plink_reference_mapping: pd.DataFrame = None,
                 user_snp_list_path: str = None,
                 reference_data_path: str = None,
                 split_harmonized_file_paths: dict[str, str] = None,
                 bed_file_paths: dict[str, str] = None,
                 vcf_file_paths: dict[str, str] = None
                 ):
        
        self.working_dir = working_dir
        self.user_upload_file = user_upload_file
        self.plink_1_9_path = plink_1_9_path
        self.plink_2_0_path = plink_2_0_path
        self.plink_map_file = plink_map_file
        self.PLINK_PREFIX = PLINK_PREFIX
        self.plink_reference_fasta = plink_reference_fasta
        self.chromosome_split_files = chromosome_split_files
        self.vcf_plink_reference_mapping = vcf_plink_reference_mapping
        self.user_snp_list_path = user_snp_list_path
        self.pvar_ref_file = pvar_ref_file
        self.reference_data_path = reference_data_path
        self.split_harmonized_file_paths = split_harmonized_file_paths
        self.bed_file_paths = bed_file_paths
        self.vcf_file_paths = vcf_file_paths
        self.beagle_references = beagle_references
        self.validate_paths()
    
    
    def validate_paths(self) -> None:
        
        for path in [
            self.working_dir,
            self.user_upload_file,
            self.plink_1_9_path,
            self.plink_2_0_path,
            self.plink_map_file,
            self.pvar_ref_file,
            self.plink_reference_fasta
        ]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required path does not exist: {path}")