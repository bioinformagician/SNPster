#!/usr/bin/env python3
"""Evaluate a frozen v3 finalist on Ped-sim-derived FLARE-clean genotypes.

The input is projected directly from high-density genomes to the FLARE marker
panel. No consumer array masking, harmonization, or Beagle are performed here.
The predicted ancestry is compared with REALIZED founder-label ancestry.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.optimize import nnls

SOURCES = ('EUR', 'EAS', 'SAS', 'SEA', 'GME', 'NAM', 'OCE', 'SSA')


def simplex_rows(x):
    x = np.atleast_2d(np.asarray(x, dtype=float))
    rows=[]
    for v in x:
        if not np.isfinite(v).all(): raise ValueError('Invalid non-finite model prediction')
        u=np.sort(v)[::-1]; cssv=np.cumsum(u)-1; j=np.arange(1,len(v)+1)
        keep=(u-cssv/j)>0
        if not np.any(keep): raise ValueError('Failed simplex projection')
        k=int(j[keep][-1]); theta=cssv[k-1]/k
        clipped=np.maximum(v-theta, 0.)
        rows.append(clipped/clipped.sum())
    return np.asarray(rows)


def decode(decoder, features):
    if isinstance(decoder, dict) and decoder.get('kind') == 'centroid_nnls':
        A=decoder['centroids'].T
        result=[]
        for x in features:
            coef,_=nnls(A,x)
            result.append(coef / coef.sum() if coef.sum()>0 else np.ones(A.shape[1])/A.shape[1])
        return np.asarray(result)
    return simplex_rows(decoder.predict(features))


def metrics(df):
    if df.empty: raise ValueError('Empty pedigree evaluation subset')
    tv=df['TV_distance'].to_numpy(float)
    return {'n':len(df), 'mean_tv':float(tv.mean()), 'median_tv':float(np.median(tv)),
            'p95_tv':float(np.quantile(tv,.95)),
            'dominant_accuracy':float(df['dominant_correct'].mean()),
            'off_target_mass':float(df['off_target_mass'].mean()),
            'source_mae':{s:float(np.abs(df[f'pred_{s}']-df[f'realized_{s}']).mean()) for s in SOURCES}}


def run(model_path, truth_path, feature_paths, out_predictions, out_metrics, candidate, representation, chromosomes,
        founder_audit_path, plan_manifest_path):
    audit = json.loads(Path(founder_audit_path).read_text())
    if (audit.get('strict_reference_disjoint') is not True
            or audit.get('labels_validated') is not True
            or audit.get('mosaic_training_exclusion_checked') is not True
            or audit.get('flare_reference_overlap_count') != 0):
        raise ValueError('Pedigree evaluation requires validated, training/reference-disjoint founders')
    manifest = json.loads(Path(plan_manifest_path).read_text())
    if manifest.get('design') not in ('full', 'reduced'):
        raise ValueError('Pedigree evaluation requires an explicitly labeled full/reduced design')
    saved=joblib.load(model_path)
    cols=saved['feature_columns']
    if tuple(saved['sources']) != SOURCES: raise ValueError('Frozen decoder source ordering differs from SNPster 8-source ontology')
    if saved.get('candidate',{}).get('candidate_id') != candidate:
        raise ValueError('Loaded model.joblib is not the requested finalist')
    if saved.get('candidate',{}).get('representation') != representation:
        raise ValueError('Loaded model.joblib representation differs from chosen FLARE model')
    if not feature_paths: raise ValueError('No FLARE features provided')
    fdfs=[]
    for path in feature_paths:
        f=pd.read_csv(path, sep='\t')
        if 'sample_id' not in f or not set(cols).issubset(f):
            raise ValueError(f'{path}: features do not match saved model.feature_columns')
        if f.sample_id.duplicated().any(): raise ValueError(f'{path}: duplicate sample_id')
        fdfs.append(f[['sample_id',*cols]])
    features=pd.concat(fdfs,ignore_index=True)
    if features.sample_id.duplicated().any(): raise ValueError('Duplicate offspring across FLARE feature batches')
    truth=pd.read_csv(truth_path, sep='\t')
    required=['sample_id','family_id','class','generation','is_founder','minimum_nonzero_ancestry'] + [f'realized_{s}' for s in SOURCES]
    if not set(required).issubset(truth): raise ValueError('Missing essential columns in pedigree truth')
    truth=truth[truth.is_founder.astype(str).str.lower().isin(['0','false'])].copy()
    if truth.sample_id.duplicated().any(): raise ValueError('Duplicate truth sample_id')
    merged=features.merge(truth,how='outer',on='sample_id',validate='one_to_one',indicator=True)
    if (merged['_merge']!='both').any():
        raise ValueError('Feature/truth IDs differ; mismatched: '+str(merged.loc[merged['_merge']!='both','sample_id'].head(10).tolist()))
    merged=merged.drop(columns=['_merge'])
    if merged.empty: raise ValueError('No pedigree descendants to evaluate')
    x=merged[cols].to_numpy(float)
    if not np.isfinite(x).all(): raise ValueError('Non-finite ancestry features')
    raw=decode(saved['decoder'],x)
    pred=simplex_rows(saved['calibrator'].predict(raw)) if saved['calibrator'] is not None else raw
    y=merged[[f'realized_{s}' for s in SOURCES]].to_numpy(float)
    if not np.allclose(y.sum(axis=1),1.,atol=1e-5): raise ValueError('Realized ancestry does not sum to 1')
    if not np.allclose(pred.sum(axis=1),1.,atol=1e-5): raise ValueError('Predicted ancestry does not sum to 1')
    out=merged[['sample_id','family_id','class','generation','minimum_nonzero_ancestry']+[f'expected_{s}' for s in SOURCES]+[f'realized_{s}' for s in SOURCES]].copy()
    out['pedigree_design'] = manifest['design']
    out['reference_disjoint'] = True
    out['founder_labels_validated'] = True
    for j,s in enumerate(SOURCES): out[f'pred_{s}']=pred[:,j]
    out['TV_distance']=np.abs(y-pred).sum(axis=1)/2
    out['dominant_correct']=np.argmax(y,axis=1)==np.argmax(pred,axis=1)
    out['off_target_mass']=(pred * (y<.01)).sum(axis=1)
    out['n_sources']=(y>=.01).sum(axis=1)
    out['minor_ancestry_bin']=pd.cut(out['minimum_nonzero_ancestry'],
                              bins=[-1,0.02,.05,.10,.25,.50,1.01],
                              labels=['<2%','2–5%','5–10%','10–25%','25–50%','>50%']).astype(str)
    Path(out_predictions).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(out_predictions,sep='\t',index=False)
    report={'candidate_id':candidate,'representation':representation,'assay':'high_density_direct_FLARE_panel',
            'production_domain_validated':False, 'truth':'bp_weighted_realized_founder_labels',
            'pedigree_design':manifest['design'], 'reference_disjoint':True,
            'founder_labels_validated':True,
            'founder_audit':audit, 'plan_manifest':manifest,
            'chromosomes':chromosomes, 'aggregate':metrics(out),
            'by_class':{str(k):metrics(g) for k,g in out.groupby('class')},
            'by_generation':{str(k):metrics(g) for k,g in out.groupby('generation')},
            'by_minor_fraction':{str(k):metrics(g) for k,g in out.groupby('minor_ancestry_bin')}}
    Path(out_metrics).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',required=True);ap.add_argument('--truth',required=True)
    ap.add_argument('--feature-files',nargs='+',required=True)
    ap.add_argument('--candidate',required=True);ap.add_argument('--representation',required=True)
    ap.add_argument('--chromosomes',required=True)
    ap.add_argument('--founder-audit',required=True);ap.add_argument('--plan-manifest',required=True)
    ap.add_argument('--predictions',required=True);ap.add_argument('--metrics',required=True)
    a=ap.parse_args()
    data=run(a.model,a.truth,a.feature_files,a.predictions,a.metrics,a.candidate,a.representation,a.chromosomes,
             a.founder_audit,a.plan_manifest)
    print(json.dumps({'candidate_id':a.candidate,'n':data['aggregate']['n'],'mean_tv':data['aggregate']['mean_tv']}))
if __name__=='__main__': main()
