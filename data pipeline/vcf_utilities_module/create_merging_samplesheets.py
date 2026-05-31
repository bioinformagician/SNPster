import argparse
from vcf_classes import VCFUtilities


#add vcf_sheet_path to environment handler and set it in main.py when creating combined samplesheet, then use it in vcf handler to set vcf file dict

args = argparse.ArgumentParser(description="Combine VCF files for PGS calculation.")
args.add_argument("--vcf_file_dir", default="/home/frederik/github_projects/SNPster/data pipeline/vcf_merge_sheets_export/vcf_merge_sheet_chr22.csv", help="Path to the CSV file containing VCF file information for merging.")
args.add_argument("--output_dir", default="/home/frederik/github_projects/SNPster/data pipeline/vcf_merge_sheets_export", help="Directory for combined VCF output.")



args = args.parse_args()


vcf_utilities = VCFUtilities()


vcf_utilities.make_vcf_merging_samplesheet(args.vcf_file_dir, args.output_dir)



