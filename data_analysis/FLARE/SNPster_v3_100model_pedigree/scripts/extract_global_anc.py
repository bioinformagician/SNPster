#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
df=pd.read_csv(a.input,sep=r'\s+',compression='infer'); anc=[c for c in df.columns if c.startswith('anc_')]
if not anc: raise SystemExit('No anc_* columns found')
scol='SAMPLE' if 'SAMPLE' in df.columns else df.columns[0]
X=df[anc].to_numpy(float); rs=X.sum(axis=1,keepdims=True)
if np.any(rs<=0): raise SystemExit('Zero-sum ancestry row')
X/=rs; out=pd.DataFrame(X,columns=anc); out.insert(0,'sample_id',df[scol].astype(str)); out.to_csv(a.output,sep='\t',index=False)
