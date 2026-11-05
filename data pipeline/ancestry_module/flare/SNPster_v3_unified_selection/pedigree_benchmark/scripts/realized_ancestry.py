#!/usr/bin/env python3
"""Infer exact source ancestry by tracing Ped-sim's --bp_all founder haplotypes.

Uses physical base pairs rather than observed variant counts. The BP records
contain [chrom]|[start] and [founder_haplotype_id]:[inclusive_end] tokens.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

SOURCES = ('EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA')


def load_bp(bp_path, target_chr):
    """Return person -> hap -> list[(start,end,founder_hap_id)]."""
    target_chr=str(target_chr).removeprefix('chr')
    parsed=defaultdict(dict)
    with open(bp_path,encoding='utf-8') as fh:
        for raw in fh:
            fields=raw.strip().split()
            if len(fields)<5: continue
            person, sex, hap = fields[:3]
            if hap not in ('h0','h1'): raise ValueError(f'Unexpected BP hap: {hap}')
            current_chr=None; start=None; parts=[]
            for tok in fields[3:]:
                if '|' in tok:
                    chrom, start_s = tok.split('|',1)
                    current_chr=chrom.removeprefix('chr')
                    start=int(start_s)
                elif ':' in tok:
                    donor_hap,end_s=tok.split(':',1)
                    end=int(end_s)
                    if end<start: raise ValueError(f'Invalid BP segment: {raw}')
                    if current_chr==target_chr:
                        parts.append((start,end,int(donor_hap)))
                    start=end+1
                else: raise ValueError(f'Unknown BP item {tok}')
            if parts: parsed[person][hap]=parts
    return parsed


def ancestry_for_batch(plan_path,bp_path,chrom,output):
    batch=json.loads(Path(plan_path).read_text())
    bp=load_bp(bp_path,chrom)
    founder_map=batch['founders']
    hap_ancestry={}
    for person,info in founder_map.items():
        if person not in bp or set(bp[person])!={'h0','h1'}:
            raise ValueError(f'Founder BP rows absent for {person}; require Ped-sim --bp_all')
        for hap in ('h0','h1'):
            segments=bp[person][hap]
            unique=set(hid for _,_,hid in segments)
            if len(unique)!=1:
                raise ValueError(f'Founder {person} {hap} unexpectedly has multiple founder hap IDs: {unique}')
            hid=next(iter(unique))
            if hid in hap_ancestry: raise ValueError(f'Founder hap ID reused: {hid}')
            hap_ancestry[hid]=info['source']
    if len(hap_ancestry)!=2*len(founder_map):
        raise ValueError('Incomplete founder haplotype mapping')
    fieldnames=['sample_id','family_id','class','generation','chromosome','bp_weight',
                'is_founder','father_id','mother_id']+[f'expected_{s}' for s in SOURCES]+[f'realized_{s}' for s in SOURCES]
    with open(output,'w',newline='',encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=fieldnames,delimiter='\t');writer.writeheader()
        for node in batch['nodes']:
            sample=node['sample_id']
            if sample not in bp or set(bp[sample])!={'h0','h1'}:
                raise ValueError(f'No complete BP rows for {sample}')
            ancestries=defaultdict(int);denom=0
            for hap in ('h0','h1'):
                segs=bp[sample][hap]
                # Every haplotype must cover the same physical span. Contiguous BP intervals.
                for i, (start,end,hid) in enumerate(segs):
                    if i and start!=segs[i-1][1]+1:
                        raise ValueError(f'BP gap or overlap for {sample} {hap} {chrom}')
                    if hid not in hap_ancestry:
                        raise ValueError(f'Unassigned founder hap {hid} for {sample}')
                    n=end-start+1
                    denom+=n
                    ancestries[hap_ancestry[hid]]+=n
            if denom<1: raise ValueError(f'Empty genomic span for {sample}')
            row={'sample_id':sample,'family_id':node['family_id'],'class':node['class'],
                 'generation':node['generation'],'chromosome':str(chrom),
                 'bp_weight':denom,'is_founder':int(node['is_founder']),
                 'father_id':node['father_id'],'mother_id':node['mother_id']}
            for s in SOURCES:
                row[f'expected_{s}']=node['expected'].get(s,0.)
                row[f'realized_{s}']=ancestries[s]/denom
            if abs(sum(float(row[f'realized_{s}']) for s in SOURCES)-1)>1e-8:
                raise ValueError(f'Ancestry sum not one: {sample}')
            writer.writerow(row)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--plan',required=True)
    p.add_argument('--bp',required=True)
    p.add_argument('--chr',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    ancestry_for_batch(a.plan,a.bp,a.chr,a.output)
if __name__=='__main__':main()
