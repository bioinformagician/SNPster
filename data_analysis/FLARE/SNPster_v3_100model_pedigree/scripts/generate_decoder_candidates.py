#!/usr/bin/env python3
import argparse,hashlib,itertools,json
from pathlib import Path
import pandas as pd

def combos(grid):
    if not grid: return [{}]
    keys=list(grid); return [dict(zip(keys,v)) for v in itertools.product(*(grid[k] for k in keys))]

ap=argparse.ArgumentParser(); ap.add_argument('--selected-reps',required=True); ap.add_argument('--search-space',required=True); ap.add_argument('--model-config-dir',required=True); ap.add_argument('--outdir',required=True); ap.add_argument('--max-candidates',type=int,default=0); a=ap.parse_args()
s=json.loads(Path(a.search_space).read_text()); st=s['stage2']; reps=pd.read_csv(a.selected_reps,sep='\t')['rep_id'].astype(str).tolist(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); model_dir=Path(a.model_config_dir)
fixed_rep_meta={c.get('representation'):c.get('representation_config',{}) for c in st.get('fixed_candidates',[])}

def rep_config(rep):
    p=model_dir/f'{rep}.json'
    if p.exists():
        return json.loads(p.read_text())
    if rep in fixed_rep_meta and fixed_rep_meta[rep]:
        return fixed_rep_meta[rep]
    return {'strategy_name':'fixed_or_external','family':'fixed','rep_id':rep}
cands=[]
for rep in reps:
  for dec in st['decoders']:
    for dp in combos(dec.get('grid',{})):
      for cal in st['calibrators']:
        for cp in combos(cal.get('grid',{})):
          for wt in st['training_weights']:
            c={'representation':rep,'representation_config':rep_config(rep),'seed':s.get('seed',20261005),'decoder':{'type':dec['type'],'params':dp},'calibrator':{'type':cal['type'],'params':cp},'training_weights':wt,'calibration_folds':st.get('calibration_folds',3)}
            raw=json.dumps(c,sort_keys=True,separators=(',',':')); h=hashlib.sha1(raw.encode()).hexdigest()[:10]
            c['candidate_id']=f"{rep}__{dec['type']}__{cal['type']}__{h}"; c['config_sha256']=hashlib.sha256(raw.encode()).hexdigest(); cands.append(c)
# Exact fixed candidates are always included as regression baselines.
for c in st.get('fixed_candidates',[]):
    if c['representation'] in reps:
        c=dict(c)
        c.setdefault('representation_config', rep_config(c['representation']))
        cands.append(c)
# Deterministic order but always preserve fixed candidates when limiting.
fixed_ids={x['candidate_id'] for x in st.get('fixed_candidates',[])}
regular=[x for x in cands if x['candidate_id'] not in fixed_ids]; fixed=[x for x in cands if x['candidate_id'] in fixed_ids]
# Deterministic hash ordering within each representation, then round-robin across
# representations so a smoke cap exercises several representations and decoders.
def stable_key(c):
    raw=json.dumps(c,sort_keys=True,separators=(',',':'))
    return hashlib.sha1(raw.encode()).hexdigest()
by_rep={}
for c in regular:
    by_rep.setdefault(c['representation'], []).append(c)
for rep in by_rep:
    by_rep[rep]=sorted(by_rep[rep], key=stable_key)
rep_order=sorted(by_rep)
regular_rr=[]
depth=0
while True:
    added=False
    for rep in rep_order:
        if depth < len(by_rep[rep]):
            regular_rr.append(by_rep[rep][depth])
            added=True
    if not added:
        break
    depth += 1
regular=regular_rr
if a.max_candidates and a.max_candidates>0:
    regular=regular[:a.max_candidates]
cands=regular+fixed
for c in cands: (out/f"{c['candidate_id']}.json").write_text(json.dumps(c,indent=2,sort_keys=True)+'\n')
print(f"Generated {len(cands)} stage2 candidates across {len(reps)} representations")
