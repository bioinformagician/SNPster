#!/usr/bin/env python3
"""Rank frozen finalists on an independent simulated-genome benchmark."""
import argparse, json, csv
from pathlib import Path

def summarize(paths,out):
    rows=[]
    for p in paths:
        m=json.loads(Path(p).read_text())
        x=m['aggregate']
        rows.append({'candidate_id':m['candidate_id'], 'representation':m['representation'],
                     'n':x['n'],'mean_tv':x['mean_tv'],'p95_tv':x['p95_tv'],
                     'dominant_accuracy':x['dominant_accuracy'],'off_target_mass':x['off_target_mass'],
                     'assay':m['assay'],'production_domain_validated':m['production_domain_validated']})
    if not rows: raise ValueError('No Ped-sim finalist metrics available')
    rows.sort(key=lambda x: (x['mean_tv'],x['p95_tv']))
    with open(out,'w',newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metrics',nargs='+',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();summarize(a.metrics,a.output)
