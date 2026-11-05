#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--root',required=True)
ap.add_argument('--search-space',required=True)
ap.add_argument('--fixed-models',required=True)
ap.add_argument('--model-config-dir',required=True)
ap.add_argument('--output',required=True)
a=ap.parse_args()

search=json.loads(Path(a.search_space).read_text())
stage1=search['stage1']
n=int(stage1['retain_top_representations'])
per_strategy=int(stage1.get('retain_best_per_strategy',0))
model_dir=Path(a.model_config_dir)

rows=[]
for p in Path(a.root).glob('screen_*/metrics.json'):
    m=json.loads(p.read_text())
    rep=m['representation']
    mp=model_dir/f'{rep}.json'
    if mp.exists():
        rc=json.loads(mp.read_text())
        strategy=rc.get('strategy_name','unknown')
        family=rc.get('family','unknown')
        k=rc.get('k')
    else:
        strategy='fixed_or_external'; family='fixed'; k=None
    rows.append({
        'rep_id':rep,
        'strategy':strategy,
        'family':family,
        'k':k,
        'pass':bool(m.get('pass',False)),
        'score':float(m.get('selection_score',999))
    })

if not rows:
    raise SystemExit('No stage1 screen metrics found')

df=(pd.DataFrame(rows)
      .drop_duplicates('rep_id')
      .sort_values(['pass','score'],ascending=[False,True]))

keep=[]
# Preserve strategy diversity before overall ranking.
if per_strategy > 0:
    for strategy in sorted(df['strategy'].dropna().unique()):
        g=df[df.strategy.eq(strategy)].sort_values(['pass','score'],ascending=[False,True])
        for rep in g.head(per_strategy)['rep_id']:
            if rep not in keep:
                keep.append(rep)

# Fill remaining capacity by global rank. If strategy preservation exceeds n,
# keep them all: the point of this stage is to avoid prematurely killing an
# entire model-building strategy.
for rep in df['rep_id']:
    if len(keep) >= n and len(keep) >= per_strategy * df['strategy'].nunique():
        break
    if rep not in keep:
        keep.append(rep)

fixed=pd.read_csv(a.fixed_models,sep='\t',dtype=str).fillna('')
for _,r in fixed.iterrows():
    if r.get('enabled','true').lower()!='false' and r.get('force_keep','false').lower()=='true' and r.rep_id not in keep:
        keep.append(r.rep_id)

out=df[df.rep_id.isin(keep)].copy()
missing=[x for x in keep if x not in set(out.rep_id)]
for x in missing:
    out=pd.concat([out,pd.DataFrame([{
      'rep_id':x,'strategy':'fixed_or_external','family':'fixed','k':None,'pass':False,'score':999.0
    }])],ignore_index=True)

# Preserve desired keep ordering.
rank={rep:i+1 for i,rep in enumerate(keep)}
out['selected_rank']=out.rep_id.map(rank)
out=out.sort_values('selected_rank')
out.to_csv(a.output,sep='\t',index=False)
print(out.to_string(index=False))
