import argparse
import os
import sys
from ancestry_classes import AncestryEnvironmentHandler, AncestryInference

parser = argparse.ArgumentParser(description='Infer genetic ancestry using 1000 Genomes reference')
parser.add_argument('--vcf_file', type=str, required=True, help='Input VCF file from harmonization')
parser.add_argument('--output_dir', type=str, default='.', help='Output directory')
parser.add_argument('--chromosome', type=str, default='22', help='Chromosome to use (default: 22 for speed)')
parser.add_argument('--reference_vcf_dir', type=str, 
                    default='/data',
                    help='Directory containing 1000G reference VCFs')
parser.add_argument('--population_panel', type=str,
                    default='/data/1000G_panel.txt',
                    help='1000 Genomes population panel file')

args = parser.parse_args()

try:
    # Initialize environment
    env_handler = AncestryEnvironmentHandler(
        vcf_file=args.vcf_file,
        output_dir=args.output_dir,
        reference_vcf_dir=args.reference_vcf_dir,
        population_panel_file=args.population_panel,
        use_chromosome=args.chromosome
    )

    # Run ancestry inference
    print("Starting ancestry inference...")
    print(f"Sample VCF: {args.vcf_file}")
    print(f"Reference directory: {args.reference_vcf_dir}")
    print(f"Using chromosome: {args.chromosome}")

    ancestry_inference = AncestryInference(env_handler)
    ancestry_inference.run_ancestry_inference()

    print("\n✓ Ancestry inference complete!")

    ancestry_inference.upload_ancestry_results()
    
except Exception as e:
    print(f"\n✗ Ancestry inference failed: {e}", file=sys.stderr)
    sys.exit(1)