import argparse
import os
import sys
from ancestry_classes import AncestryEnvironmentHandler, AncestryInference

parser = argparse.ArgumentParser(description='Infer genetic ancestry using 1000 Genomes reference')
parser.add_argument('--vcf_file', type=str, required=True, help='Input VCF file from harmonization')
parser.add_argument('--output_dir', type=str, default='.', help='Output directory')
parser.add_argument('--chromosomes', type=str, default='1,2,21,22', 
                    help='Comma-separated list of chromosomes to use (default: 1,2,21,22)')
parser.add_argument('--reference_vcf_dir', type=str, 
                    default=os.getenv('REFERENCE_VCF_DIR', '/data/references'),
                    help='Directory containing reference VCFs')
parser.add_argument('--population_panel', type=str,
                    default=os.getenv('POPULATION_PANEL_FILE', '/data/hgdp_1kg_panel.txt'),
                    help='Population panel file')

args = parser.parse_args()

# Parse chromosome list
chromosomes = [c.strip() for c in args.chromosomes.split(',')]

try:
    # Initialize environment
    env_handler = AncestryEnvironmentHandler(
        vcf_file=args.vcf_file,
        output_dir=args.output_dir,
        reference_vcf_dir=args.reference_vcf_dir,
        population_panel_file=args.population_panel,
        use_chromosomes=chromosomes
    )

    # Run ancestry inference
    print("Starting ancestry inference...")
    print(f"Sample VCF: {args.vcf_file}")
    print(f"Reference directory: {args.reference_vcf_dir}")
    print(f"Using chromosomes: {', '.join(chromosomes)}")

    ancestry_inference = AncestryInference(env_handler)
    ancestry_inference.run_ancestry_inference()

    print("\n✓ Ancestry inference complete!")

    ancestry_inference.upload_ancestry_results()
    
except Exception as e:
    print(f"\n✗ Ancestry inference failed: {e}", file=sys.stderr)
    sys.exit(1)