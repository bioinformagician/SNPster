from pathlib import Path
from combiner_classes import Combiner
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=False, default=os.getcwd())
parser.add_argument(
    '--parquet_files',
    type=str,
    nargs='+',
    required=True,
    help='One or more standardized parquet files to combine.'
)
parser.add_argument(
    '--imputation_id',
    type=str,
    required=False,
    default=None,
    help='Optional imputation id. If omitted, read from parquet metadata.',
)

args = parser.parse_args()

combiner = Combiner(
    file_paths=[Path(file_path) for file_path in args.parquet_files],
    output_dir=args.output_dir,
    imputation_id=args.imputation_id,
)

combiner.combine_microarray_data()
combiner.validate_combined_data()
combiner.write_parquet_output()




