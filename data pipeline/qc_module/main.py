import argparse
import os
import pandas as pd
from qc_classes import *
from config import *

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=False, default=os.getcwd())
parser.add_argument('--input_file', type=str, required=False, default = os.getenv("VCF_FILES_TEST_DIR"), help='Path to the user microarray data file (e.g., 23andMe, Ancestry, Myheritage...).')


args = parser.parse_args()


qc_thresholds = QCThresholds(GP_MIN, DS_TOL, SNPS_ONLY, BIALLELIC_ONLY)

imputed_data_container = ImputedDataContainer(qc_thresholds = qc_thresholds,
                                              file_path=args.input_file
                                              )




if ENGINE == "pandas":
    imputed_data_container.load_vcf_to_df_pandas()
else:
    imputed_data_container.load_vcf_to_df_polars()

print(f"Running QC on {imputed_data_container.file_path} imputed data...")
imputed_data_container.qc_imputed_data()

vcf_out_path = imputed_data_container.write_pandas_to_vcf(output_dir = args.output_dir)

print(f"Completed chromosome {imputed_data_container.file_path} and written to {args.output_dir}. Memory freed.")

imputed_data_container.zip_vcf(vcf_path = vcf_out_path)

