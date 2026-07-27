from typing import Optional
from dataclasses import dataclass
import pandas as pd


@dataclass
class DataContainer:
    imputation_id: Optional[str] = None
    file_id: Optional[str] = None
    vendor: Optional[str] = None
    genome_build: Optional[str] = None
    is_forward_strand: Optional[bool] = True
    harmonized_data: Optional[pd.DataFrame] = None
    


