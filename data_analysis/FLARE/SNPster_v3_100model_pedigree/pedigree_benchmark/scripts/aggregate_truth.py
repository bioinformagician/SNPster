#!/usr/bin/env python3
"""Combine per-chromosome realized truth with physical-span weighting and checks."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
SOURCES=('EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA')

def run(files, out, expected_chromosomes):
    records={}
    seen=defaultdict(set)
    for fname in files:
        with open(fname, newline='') as fh:
            for r in csv.DictReader(fh,delimiter='\t'):
                key=r['sample_id'];chrom=r['chromosome']
                if chrom in seen[key]: raise ValueError(f'Duplicate chromosome {chrom} for {key}')
                seen[key].add(chrom)
                if key not in records:
                    records[key]={'family_id':r['family_id'],'class':r['class'],
                                  'generation':int(r['generation']),'is_founder':int(r['is_founder']),
                                  'father_id':r['father_id'],'mother_id':r['mother_id'],
                                  'expected':[float(r[f'expected_{s}']) for s in SOURCES],
                                  'realized_bases':[0.]*len(SOURCES),'bases':0}
                v=records[key]
                if any(abs(float(r[f'expected_{s}'])-v['expected'][i])>1e-10 for i,s in enumerate(SOURCES)):
                    raise ValueError(f'Inconsistent expectations for {key}')
                w=int(r['bp_weight']);v['bases']+=w
                for i,s in enumerate(SOURCES): v['realized_bases'][i]+=float(r[f'realized_{s}'])*w
    expected_set={str(c).removeprefix('chr') for c in expected_chromosomes}
    for sample, observed in seen.items():
        if observed!=expected_set: raise ValueError(f'{sample}: chromosomes expected {expected_set} but got {observed}')
    fields=['sample_id','family_id','class','generation','is_founder','father_id','mother_id','bp_weight',
            'n_sources','minimum_nonzero_ancestry','tv_expected_vs_realized'] + [f'expected_{s}' for s in SOURCES] + [f'realized_{s}' for s in SOURCES]
    with open(out,'w',newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');wr.writeheader()
        for sample,v in sorted(records.items()):
            realized=[x/v['bases'] for x in v['realized_bases']]
            nonzero=[x for x in realized if x>1e-10]
            row={'sample_id':sample,'family_id':v['family_id'],'class':v['class'],
                 'generation':v['generation'],'is_founder':v['is_founder'],
                 'father_id':v['father_id'],'mother_id':v['mother_id'],
                 'bp_weight':v['bases'],'n_sources':len(nonzero),
                 'minimum_nonzero_ancestry':min(nonzero) if nonzero else 0,
                 'tv_expected_vs_realized':0.5*sum(abs(a-b) for a,b in zip(realized,v['expected']))}
            row.update({f'expected_{s}':v['expected'][i] for i,s in enumerate(SOURCES)})
            row.update({f'realized_{s}':realized[i] for i,s in enumerate(SOURCES)})
            wr.writerow(row)
    return len(records)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--inputs',nargs='+',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--chromosomes',required=True,help='Comma-separated list (1,2,3)')
    a=p.parse_args()
    run([Path(x) for x in a.inputs],a.output,a.chromosomes.split(','))
if __name__=='__main__':main()
