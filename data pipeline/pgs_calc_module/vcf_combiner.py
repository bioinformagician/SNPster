import argparse
import os


import pandas as pd

from .vcf_combiner_classes import VCFEnvironmentHandler, VCFHandler


#add vcf_sheet_path to environment handler and set it in main.py when creating combined samplesheet, then use it in vcf handler to set vcf file dict

args = argparse.ArgumentParser(description="Combine VCF files for PGS calculation.")
args.add_argument("--output_dir", required=True, help="Directory for combined VCF output.")
args.add_argument("--output_dir", required=True, help="Directory for output of merged vcf files")



