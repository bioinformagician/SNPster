#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--root',required=True)
ap.add_argument('--outdir',required=True)
ap.add_argument('--title',default='SNPster ancestry strategy benchmark v3')
a=ap.parse_args()

rows=[]; raw=[]
for p in sorted(Path(a.root).glob('candidate_*/metrics.json')):
    m=json.loads(p.read_text()); raw.append(m)
    cand=m['candidate']; rc=cand.get('representation_config',{}) or {}
    sel=m.get('role_metrics',{}).get('validation',{})
    lonni=next((r for r in m.get('real_world_benchmarks',[]) if r.get('benchmark_id','').startswith('lonni')),None)
    dec=cand.get('decoder',{}); dp=dec.get('params',{}) or {}
    cal=cand.get('calibrator',{}); cp=cal.get('params',{}) or {}
    wt=cand.get('training_weights',{}) or {}
    row={
      'candidate_id':m['candidate_id'],
      'representation':m.get('representation'),
      'strategy':rc.get('strategy_name'),
      'model_family':rc.get('family'),
      'K_requested':rc.get('k'),
      'K_effective':len(m.get('feature_columns',[]) or []),
      'covariance':rc.get('covariance_type'),
      'pca_components':rc.get('pca_components',0),
      'refinement_mode':rc.get('refinement_mode'),
      'refine_em':rc.get('refine_em'),
      'weight_prior_type':rc.get('weight_concentration_prior_type'),
      'weight_prior':rc.get('weight_concentration_prior'),
      'model_seed':rc.get('model_seed'),
      'panel_cm':rc.get('panel_cm'),
      'decoder':dec.get('type'),
      'trees':dp.get('n_estimators'),
      'leaf':dp.get('min_samples_leaf'),
      'max_features':dp.get('max_features'),
      'decoder_alpha':dp.get('alpha'),
      'calibrator':cal.get('type'),
      'calibrator_alpha':cp.get('alpha'),
      'pure_weight':wt.get('pure'),
      'mosaic_weight':wt.get('mosaic'),
      'pass':m.get('pass',False),
      'real_hard_pass':m.get('real_hard_pass',False),
      'selection_score':m.get('selection_score'),
      'mean_tv':sel.get('mean_tv'),
      'p95_tv':sel.get('p95_tv'),
      'off_target_mass':sel.get('off_target_mass'),
      'dominant_accuracy':sel.get('dominant_accuracy')
    }
    if lonni:
        row['lonni_pass']=lonni.get('pass')
        pred=lonni.get('prediction',{})
        row.update({f'lonni_{k}':v for k,v in pred.items()})
        row['lonni_EAS_SEA']=(pred.get('EAS',0) or 0)+(pred.get('SEA',0) or 0)
    rows.append(row)

if not rows: raise SystemExit('No candidate results found')
df=pd.DataFrame(rows).sort_values(['pass','real_hard_pass','selection_score'],ascending=[False,False,True],na_position='last')
if 'lonni_pass' not in df.columns:
    df['lonni_pass']=False
out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
df.to_csv(out/'leaderboard.tsv',sep='\t',index=False)
(out/'all_metrics.json').write_text(json.dumps(raw,indent=2)+'\n')

# Strategy-level summaries make it easy to compare model-building approaches,
# not just individual hyperparameter combinations.
summary=(df.groupby(['strategy','model_family'],dropna=False)
           .agg(
             n_candidates=('candidate_id','size'),
             n_pass=('pass','sum'),
             n_lonni_pass=('lonni_pass','sum'),
             best_selection_score=('selection_score','min'),
             median_selection_score=('selection_score','median'),
             best_mean_tv=('mean_tv','min'),
             median_mean_tv=('mean_tv','median'),
             best_p95_tv=('p95_tv','min'),
             median_p95_tv=('p95_tv','median'),
             best_off_target=('off_target_mass','min'),
             best_dominant_accuracy=('dominant_accuracy','max')
           )
           .reset_index()
           .sort_values(['n_pass','best_selection_score'],ascending=[False,True]))
summary.to_csv(out/'strategy_summary.tsv',sep='\t',index=False)

family=(df.groupby('model_family',dropna=False)
          .agg(
            n_candidates=('candidate_id','size'),
            n_pass=('pass','sum'),
            n_lonni_pass=('lonni_pass','sum'),
            best_selection_score=('selection_score','min'),
            median_mean_tv=('mean_tv','median'),
            best_mean_tv=('mean_tv','min')
          )
          .reset_index()
          .sort_values(['n_pass','best_selection_score'],ascending=[False,True]))
family.to_csv(out/'model_family_summary.tsv',sep='\t',index=False)
view=df.head(1000)
html=view.to_html(index=False,float_format=lambda x:f'{x:.5f}')
(out/'leaderboard.html').write_text(
    f'<!doctype html><meta charset="utf-8"><title>{a.title}</title>'
    '<style>body{font-family:system-ui;margin:24px}table{border-collapse:collapse;font-size:12px}'
    'th,td{border:1px solid #ddd;padding:5px 7px}th{position:sticky;top:0;background:#eee}</style>'
    f'<h1>{a.title}</h1>{html}'
)
print(df.head(40).to_string(index=False))
