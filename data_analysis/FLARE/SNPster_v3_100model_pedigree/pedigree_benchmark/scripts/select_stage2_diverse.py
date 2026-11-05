#!/usr/bin/env python3
"""Broaden stage3 shortlist so every surviving FLARE representation gets pedigree validation.

Stage2 is still an economical reference-based decoder pre-screen; its score is
NOT the final model-selection objective. Force at least one candidate per
representation and preserve frozen controls before filling remaining slots.
"""
import argparse, json, shutil
from pathlib import Path
import pandas as pd


def select(root, search_space, output_dir, output_table, manifest=None):
    search=json.loads(Path(search_space).read_text())
    n=int(search['stage2']['retain_top_candidates'])
    fixed={c['candidate_id'] for c in search['stage2'].get('fixed_candidates',[])}
    if bool(root)==bool(manifest): raise ValueError('Exactly one of --root and --manifest required')
    if manifest:
        paths=[]
        for lineno,line in enumerate(Path(manifest).read_text().splitlines(),1):
            folder=Path(line.strip())
            if not line.strip(): raise ValueError(f'Empty manifest line {lineno}')
            for filename in ('metrics.json','candidate.json'):
                f=folder/filename
                if not f.is_file() or f.stat().st_size == 0:
                    raise ValueError(f'Manifest line {lineno}: missing/empty {f}')
            paths.append(folder/'metrics.json')
        if len(paths)!=len(set(paths)): raise ValueError('Duplicate directories in result manifest')
    else:
        paths=list(Path(root).glob('candidate_*/metrics.json'))
    rows=[]
    for p in paths:
        m=json.loads(p.read_text());rows.append({'candidate_id':m['candidate_id'],
            'representation':m['representation'], 'pass':bool(m.get('pass',False)),
            'real_hard_pass':bool(m.get('real_hard_pass',False)),
            'score':float(m.get('selection_score',999)),'dir':str(p.parent)})
    if not rows: raise ValueError('No stage2 candidate results found')
    df=pd.DataFrame(rows).drop_duplicates('candidate_id')
    df=df.sort_values(['pass','real_hard_pass','score','candidate_id'],ascending=[False,False,True,True])
    ids=[]
    # Ensure a chance for high-admixture-performing representations to survive.
    for rep in sorted(df.representation.unique()):
        cid=df[df.representation.eq(rep)].iloc[0].candidate_id
        if cid not in ids: ids.append(cid)
    for cid in sorted(fixed):
        if cid in set(df.candidate_id) and cid not in ids: ids.append(cid)
    if n and len(ids)>n:
        raise ValueError(f'Stage2 retain_top_candidates={n} < number of representations and fixed controls ({len(ids)}). Increase the stage2 shortlist.')
    for cid in df.candidate_id:
        if n and len(ids)>=n: break
        if cid not in ids: ids.append(cid)
    selected=df[df.candidate_id.isin(ids)].copy()
    selected['selection_priority']=selected.candidate_id.map({cid:i for i,cid in enumerate(ids)})
    selected=selected.sort_values('selection_priority')
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    for _,r in selected.iterrows(): shutil.copy2(Path(r.dir)/'candidate.json',out/f'{r.candidate_id}.json')
    selected.drop(columns=['dir']).to_csv(output_table,sep='\t',index=False)
    return selected

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root');p.add_argument('--manifest');p.add_argument('--search-space',required=True)
    p.add_argument('--outdir',required=True);p.add_argument('--table',required=True);a=p.parse_args()
    d=select(a.root,a.search_space,a.outdir,a.table,a.manifest)
    print('Selected',len(d),'stage3 candidates from',d.representation.nunique(),'representations')
