#!/usr/bin/env python3
import argparse,json,shutil
from pathlib import Path
import pandas as pd
ap=argparse.ArgumentParser(); ap.add_argument('--root'); ap.add_argument('--manifest'); ap.add_argument('--search-space',required=True); ap.add_argument('--outdir',required=True); ap.add_argument('--table',required=True); a=ap.parse_args()
s=json.loads(Path(a.search_space).read_text()); n=int(s['stage2']['retain_top_candidates']); fixed_ids={x['candidate_id'] for x in s['stage2'].get('fixed_candidates',[])}
if bool(a.root)==bool(a.manifest): raise SystemExit('Specify exactly one of --root and --manifest')
if a.manifest:
 paths=[]
 for lineno,line in enumerate(Path(a.manifest).read_text().splitlines(),1):
  folder=Path(line.strip())
  if not line.strip(): raise SystemExit(f'Empty manifest line {lineno}')
  for filename in ('metrics.json','candidate.json'):
   f=folder/filename
   if not f.is_file() or f.stat().st_size==0: raise SystemExit(f'Manifest line {lineno}: missing/empty {f}')
  paths.append(folder/'metrics.json')
 if len(paths)!=len(set(paths)): raise SystemExit('Duplicate stage2 result directory in manifest')
else:
 paths=list(Path(a.root).glob('candidate_*/metrics.json'))
rows=[]
for p in paths:
 m=json.loads(p.read_text()); rows.append({'candidate_id':m['candidate_id'],'pass':bool(m.get('pass',False)),'real_hard_pass':bool(m.get('real_hard_pass',False)),'score':float(m.get('selection_score',999)),'dir':str(p.parent)})
if not rows: raise SystemExit('No stage2 results found')
df=pd.DataFrame(rows).sort_values(['pass','real_hard_pass','score'],ascending=[False,False,True]); keep=list(df.head(n if n>0 else len(df)).candidate_id)
for x in fixed_ids:
 if x in set(df.candidate_id) and x not in keep: keep.append(x)
outdir=Path(a.outdir); outdir.mkdir(parents=True,exist_ok=True); selected=df[df.candidate_id.isin(keep)].copy()
for _,r in selected.iterrows(): shutil.copy2(Path(r.dir)/'candidate.json',outdir/f"{r.candidate_id}.json")
selected.drop(columns=['dir']).to_csv(a.table,sep='\t',index=False); print(selected.drop(columns=['dir']).to_string(index=False))
