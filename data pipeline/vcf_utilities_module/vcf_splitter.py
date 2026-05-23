import argparse
from vcf_classes import VCFEnvironmentHandler, VCFHandler


#add vcf_sheet_path to environment handler and set it in main.py when creating combined samplesheet, then use it in vcf handler to set vcf file dict

args = argparse.ArgumentParser(description="Combine VCF files for PGS calculation.")
args.add_argument("--vcf_samplesheet_path", default="/home/frederik/github_projects/SNPster/data pipeline/vcf_merge_sheets_export/vcf_merge_sheet_chr22.csv", help="Path to the CSV file containing VCF file information for merging.")
args.add_argument("--output_dir", default="/home/frederik/github_projects/SNPster/data pipeline/vcf_merge_sheets_export", help="Directory for combined VCF output.")



args = args.parse_args()

vcf_environment_handler = VCFEnvironmentHandler(
    output_dir=args.output_dir,
    vcf_samplesheet_path=args.vcf_samplesheet_path
    )

vcf_handler = VCFHandler(vcf_environment_handler=vcf_environment_handler)


vcf_handler.split_vcf_files()


