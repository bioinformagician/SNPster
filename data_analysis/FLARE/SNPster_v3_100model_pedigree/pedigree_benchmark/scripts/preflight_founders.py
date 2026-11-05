#!/usr/bin/env python3
"""Check founder eligibility and exact plan capacity before expensive model jobs."""
import argparse
import json
import tempfile
from pathlib import Path

from make_plans import generate
from prepare_founders import check_existing_output, prepare, read_reference_panel, sample_ids


def run(params, outdir, families, batch_size, seed, reuse, design):
    check_existing_output(outdir, design)
    vcf = Path(params['pedigree_vcf_dir']) / (
        'hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz')
    with tempfile.TemporaryDirectory(prefix='snpster_founder_preflight_') as directory:
        base = Path(directory)
        donors = base / 'eligible.tsv'
        summary = prepare(
            params['pedigree_validation_meta'], params['pedigree_training_meta'],
            sample_ids(vcf), donors,
            reference_samples=sample_ids(params['flare_ref']),
            additional_metadata=params.get('pedigree_additional_meta'),
            label_metadata=params['pedigree_label_meta'],
            reference_labels=read_reference_panel(params['flare_panel']),
            mosaic_metadata=params['pedigree_mosaic_meta'],
        )
        generate(donors, [Path(params['pedigree_training_meta'])], base / 'plans',
                 families, batch_size, seed, reuse, design)
        manifest = json.loads((base / 'plans/manifest.json').read_text())
    report = {'founder_audit': summary, 'plan': manifest}
    print(json.dumps(report, indent=2, sort_keys=True))
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--families', type=int, required=True)
    parser.add_argument('--batch-size', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--allow-reuse', choices=('true', 'false'), required=True)
    parser.add_argument('--design', choices=('full', 'reduced'), required=True)
    args = parser.parse_args()
    if args.families < 9:
        parser.error('Unified selection requires at least nine families')
    run(json.loads(Path(args.params).read_text()), args.outdir, args.families,
        args.batch_size, args.seed, args.allow_reuse == 'true', args.design)


if __name__ == '__main__':
    main()
