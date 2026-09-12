import argparse
import os
import re
import sys
from ancestry_classes import AncestryEnvironmentHandler, AncestryInference

parser = argparse.ArgumentParser(description='Infer genetic ancestry using 1000 Genomes reference')
parser.add_argument('--vcf_files', type=str, nargs='+', required=True,
                    help='One or more per-chromosome VCF files (e.g. chr1.merged.vcf.gz '
                         'chr2.merged.vcf.gz ...). The chromosome is inferred from each filename.')
parser.add_argument('--output_dir', type=str, default='.', help='Output directory')
parser.add_argument('--reference_vcf_dir', type=str, 
                    default=os.getenv('REFERENCE_VCF_DIR', '/data/references'),
                    help='Directory containing reference VCFs')
parser.add_argument('--population_panel', type=str,
                    default=os.getenv('POPULATION_PANEL_FILE', '/data/hgdp_1kg_panel.txt'),
                    help='Population panel file')

args = parser.parse_args()

CHROM_FILENAME_PATTERN = re.compile(r'chr([0-9]{1,2}|X|Y|MT)\.', re.IGNORECASE)


def chrom_from_filename(path: str) -> str:
    match = CHROM_FILENAME_PATTERN.search(os.path.basename(path))
    if not match:
        raise ValueError(f"Could not infer chromosome from VCF filename: {path}")
    return match.group(1)


vcf_files_by_chrom = {chrom_from_filename(path): path for path in args.vcf_files}

try:
    # Initialize environment
    env_handler = AncestryEnvironmentHandler(
        vcf_files=vcf_files_by_chrom,
        output_dir=args.output_dir,
        reference_vcf_dir=args.reference_vcf_dir,
        population_panel_file=args.population_panel,
    )

    # Run ancestry inference
    print("Starting ancestry inference...")
    print(f"Sample VCFs: {vcf_files_by_chrom}")
    print(f"Reference directory: {args.reference_vcf_dir}")
    print(f"Using chromosomes: {', '.join(sorted(vcf_files_by_chrom, key=str))}")

    ancestry_inference = AncestryInference(env_handler)
    ancestry_inference.run_ancestry_inference()

    print("\n✓ Ancestry inference complete!")

    ancestry_inference.upload_ancestry_results()
    
except Exception as e:
    print(f"\n✗ Ancestry inference failed: {e}", file=sys.stderr)
    sys.exit(1)