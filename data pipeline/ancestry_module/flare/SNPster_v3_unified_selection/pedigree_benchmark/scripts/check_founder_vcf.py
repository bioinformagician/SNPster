#!/usr/bin/env python3
"""Reject unphased or missing founder genotypes before Ped-sim consumes a VCF.
The VCF should already be sample-subsetted, biallelic SNP-only and compressed.
"""
import argparse
import gzip
from pathlib import Path

def run(vcf, sample_file, expected_chromosome, map_path):
    expected=set(Path(sample_file).read_text().split())
    chrom=str(expected_chromosome).removeprefix('chr')
    sample_names=[]; n=0; vcf_chrom=None
    with gzip.open(vcf,'rt') if str(vcf).endswith('.gz') else open(vcf) as fh:
        for raw in fh:
            if raw.startswith('##'): continue
            if raw.startswith('#CHROM'):
                sample_names=raw.rstrip('\n').split('\t')[9:]
                if set(sample_names)!=expected: raise ValueError('Founder VCF sample IDs do not match plan')
                continue
            f=raw.rstrip('\n').split('\t')
            if not sample_names: raise ValueError('VCF header missing')
            if f[0].removeprefix('chr')!=chrom: raise ValueError(f'Wrong chromosome: {f[0]} expected {chrom}')
            vcf_chrom=f[0]
            if len(f[3])!=1 or len(f[4])!=1 or f[3] not in 'ACGT' or f[4] not in 'ACGT':
                raise ValueError(f'Expected only biallelic ACGT SNVs at {f[0]}:{f[1]}')
            gt_idx=f[8].split(':').index('GT')
            for sample, rec in zip(sample_names,f[9:]):
                gt=rec.split(':')[gt_idx]
                if gt not in ('0|0','0|1','1|0','1|1'):
                    raise ValueError(f'Missing / unphased genotype {gt} at {f[0]}:{f[1]} sample={sample}')
            n+=1
    if n==0: raise ValueError('No usable genotype sites')
    nmap=0;prev_pos=None;lastcm=None
    with open(map_path) as f:
        for raw in f:
            if not raw.strip() or raw.lstrip().startswith('#'): continue
            cols=raw.split()
            if len(cols) not in (3,4): raise ValueError('Ped-sim .simmap expects chrom bp cM [female_cM], not a PLINK .map')
            ch,pos=cols[0].removeprefix('chr'),int(cols[1]); cm=[float(x) for x in cols[2:]]
            if ch!=chrom or cols[0]!=vcf_chrom: raise ValueError(f'Ped-sim requires exact VCF/map chromosome label match: map={cols[0]} VCF={vcf_chrom}')
            if prev_pos is not None and (pos<=prev_pos or any(x<y for x,y in zip(cm,lastcm))):
                raise ValueError('Genetic map must be strictly increasing physically and nondecreasing in cM')
            nmap+=1;prev_pos=pos;lastcm=cm
    if nmap<2: raise ValueError('Need at least two positions in Ped-sim genetic map')
    print(f'Preflight passed: chr{chrom}, {len(sample_names)} donors, {n:,} SNVs, {nmap:,} map positions')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--vcf',required=True)
    p.add_argument('--samples',required=True)
    p.add_argument('--chromosome',required=True)
    p.add_argument('--map',required=True)
    a=p.parse_args();run(a.vcf,a.samples,a.chromosome,a.map)
if __name__=='__main__':main()
