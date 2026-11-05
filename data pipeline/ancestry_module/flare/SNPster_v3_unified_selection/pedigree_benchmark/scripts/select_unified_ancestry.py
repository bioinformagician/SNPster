#!/usr/bin/env python3
"""Select frozen v3 decoders using a single ancestry-proportion regression loss.

Reference labels and Ped-sim realized founder-label truth use the SAME 8-component
TV error. Balance by scenario, ancestry combination, and founder family so neither
sample counts nor many siblings dominate. Does not refit any models.
"""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd

SOURCES = ('EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA')

def reference_frame(path):
    df = pd.read_csv(path,sep='\t',dtype={'sample_id':str})
    df = df.loc[df['role'].eq('validation')].copy()
    if df.empty: raise ValueError(f'{path}: missing reference-validation predictions')
    for s in SOURCES:
        if s not in df or f'truth_{s}' not in df: raise ValueError(f'{path}: missing {s} or truth_{s}')
    out = df[['sample_id']].copy()
    for s in SOURCES:
        out[f'pred_{s}'] = df[s].astype(float)
        out[f'truth_{s}'] = df[f'truth_{s}'].astype(float)
    out['scenario'] = 'reference_source'
    out['family_id'] = 'ref:' + out['sample_id'].astype(str)
    return out

def pedigree_frame(path):
    df = pd.read_csv(path,sep='\t',dtype={'sample_id':str})
    for s in SOURCES:
        if f'pred_{s}' not in df or f'realized_{s}' not in df: raise ValueError(f'{path}: missing {s}')
    for c in ['family_id','class']:
        if c not in df: raise ValueError(f'{path}: missing {c}')
    out = df[['sample_id','family_id']].copy()
    out['scenario'] = df['class'].astype(str)
    for s in SOURCES:
        out[f'pred_{s}'] = df[f'pred_{s}'].astype(float)
        out[f'truth_{s}'] = df[f'realized_{s}'].astype(float)
    return out

def validate_sample_vectors(df, label):
    pred = df[[f'pred_{s}' for s in SOURCES]].to_numpy(float)
    truth = df[[f'truth_{s}' for s in SOURCES]].to_numpy(float)
    if not len(df) or df[['sample_id']].duplicated().any():
        raise ValueError(f'{label}: empty or duplicated sample ids')
    if not (np.isfinite(pred).all() and np.isfinite(truth).all()):
        raise ValueError(f'{label}: non-finite fractions')
    if np.any(pred < -1e-7) or np.any(truth < -1e-7) or not np.allclose(pred.sum(axis=1),1,atol=2e-4) or not np.allclose(truth.sum(axis=1),1,atol=2e-4):
        raise ValueError(f'{label}: predicted/truth vectors must sum to 1 and be nonnegative')

def weighted_quantile(values, weights, q):
    """Upper empirical quantile under nonnegative per-row weights."""
    values=np.asarray(values,float);weights=np.asarray(weights,float)
    if not len(values) or not np.isfinite(weights).all() or np.any(weights < 0) or weights.sum()<=0:
        raise ValueError('Invalid weighted quantile input')
    order=np.argsort(values,kind='mergesort')
    return float(values[order][np.searchsorted(np.cumsum(weights[order]),q*weights.sum(),side='left')])

def calculate(df, cfg):
    validate_sample_vectors(df, 'combined')
    d=df.copy()
    p=d[[f'pred_{s}' for s in SOURCES]].to_numpy(float)
    y=d[[f'truth_{s}' for s in SOURCES]].to_numpy(float)
    d['tv']=.5*np.abs(p-y).sum(axis=1)
    d['off_target']=(p*(y < float(cfg.get('absent_ancestry_max',.01)))).sum(axis=1)
    d['dominant_correct']=(p.argmax(axis=1)==y.argmax(axis=1)).astype(float)
    # Broad category signature only for balancing (not part of the truth definition).
    d['combination']=['+'.join(SOURCES[i] for i,v in enumerate(row) if v>=.01) for row in y]
    assert not d['combination'].eq('').any()
    scenarios=[]
    for scenario,g in d.groupby('scenario',sort=True):
        combos=[]
        row_weights=pd.Series(0.,index=g.index)
        ncomb=g['combination'].nunique()
        for combination,group in g.groupby('combination'):
            # Equal weight per independent founder family within each ancestry combination.
            # Each combination and each founder family receives equal total mass.
            nfamilies=group['family_id'].nunique()
            for _, family_rows in group.groupby('family_id'):
                row_weights.loc[family_rows.index] = 1.0/(ncomb*nfamilies*len(family_rows))
            fam=group.groupby('family_id').agg(tv=('tv','mean'), off=('off_target','mean'),
                                                dom=('dominant_correct','mean'))
            combos.append({'combination':combination,'n':len(group),'families':len(fam),
                           'mean_tv':float(fam.tv.mean()),'off_target':float(fam.off.mean()),
                           'dominant_accuracy':float(fam.dom.mean())})
        scenarios.append({'scenario':scenario,'n':len(g),'combinations':len(combos),
                          'mean_tv':float(np.mean([r['mean_tv'] for r in combos])),
                          'p95_tv':weighted_quantile(g.tv.to_numpy(),row_weights.loc[g.index].to_numpy(),.95),
                          'off_target':float(np.mean([r['off_target'] for r in combos])),
                          'dominant_accuracy':float(np.mean([r['dominant_accuracy'] for r in combos])),
                          'by_combination':combos})
    # Every benchmark class is one scenario; neither 159 pure people nor 5k siblings
    # can dominate the aggregate just because there are more rows.
    mean_tv=float(np.mean([r['mean_tv'] for r in scenarios]))
    p95=float(np.mean([r['p95_tv'] for r in scenarios]))
    off=float(np.mean([r['off_target'] for r in scenarios]))
    dom=float(np.mean([r['dominant_accuracy'] for r in scenarios]))
    minor_mask=(y>=float(cfg.get('minor_min',.02))) & (y<=float(cfg.get('minor_max',.05)))
    minor_n=int(minor_mask.sum())
    minor_detect=(p[minor_mask]>=float(cfg.get('minor_detection_threshold',.01))) if minor_n else np.array([])
    false_mask=y<float(cfg.get('absent_ancestry_max',.01))
    return {
        'n':len(d),'n_reference':int(d.scenario.eq('reference_source').sum()),
        'n_pedigree':int(d.scenario.ne('reference_source').sum()),
        'scenario_count':len(scenarios),'balanced_mean_tv':mean_tv,
        'balanced_p95_tv':p95,'balanced_off_target_mass':off,
        'balanced_dominant_accuracy':dom,'unstratified_mean_tv':float(d.tv.mean()),
        'minor_2_5_n':minor_n,
        'minor_2_5_recall':float(minor_detect.mean()) if minor_n else None,
        'false_component_calls_5pct':int((p[false_mask]>=.05).sum()),
        'scenarios':scenarios,
        'per_source_mae':{s:float(np.mean(np.abs(p[:,i]-y[:,i]))) for i,s in enumerate(SOURCES)}
    }

def choose(candidate_dirs,ped_paths,config,outdir):
    cfg=json.loads(Path(config).read_text())
    dirs={}
    for path in candidate_dirs:
        path=Path(path)
        m=json.loads((path/'metrics.json').read_text())
        cid=str(m['candidate_id'])
        if cid in dirs: raise ValueError(f'Duplicate stage3 candidate ID: {cid}')
        dirs[cid]=(path,m)
    peds={}
    for path in ped_paths:
        path=Path(path)
        cid=path.name.removesuffix('.pedigree_predictions.tsv')
        if cid in peds: raise ValueError(f'Duplicate pedigree predictions: {cid}')
        peds[cid]=path
    if not dirs or set(dirs)!=set(peds):
        raise ValueError(f'Candidate/pedigree ID mismatch. Missing pedigree: {sorted(set(dirs)-set(peds))}, extra pedigree: {sorted(set(peds)-set(dirs))}')
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    rows=[]; reports=[]
    for cid,(directory,legacy) in sorted(dirs.items()):
        pure=reference_frame(directory/'predictions.tsv')
        admixed=pedigree_frame(peds[cid])
        if set(pure.sample_id) & set(admixed.sample_id): raise ValueError(f'{cid}: overlapping reference and pedigree sample IDs')
        frame=pd.concat([pure,admixed],ignore_index=True)
        m=calculate(frame,cfg)
        if m['n_reference'] < cfg.get('min_reference_samples',1) or m['n_pedigree'] < cfg.get('min_pedigree_samples',1):
            raise ValueError(f'{cid}: insufficient samples for unified evaluation')
        required=set(cfg.get('required_scenarios',[]))
        have={r['scenario'] for r in m['scenarios']}
        if not required.issubset(have): raise ValueError(f'{cid}: missing required simulation scenarios {sorted(required-have)}')
        gate=[]
        for metric,maximum in cfg.get('max_gates',{}).items():
            gate.append((metric,float(m[metric])<=float(maximum)))
        for metric,minimum in cfg.get('min_gates',{}).items():
            gate.append((metric,float(m[metric])>=float(minimum)))
        n_support=int(cfg.get('minor_gate_min_support',20))
        if m['minor_2_5_n'] >= n_support:
            gate.append(('minor_2_5_recall',m['minor_2_5_recall']>=float(cfg.get('minor_recall_min',.40))))
        else:
            # Missing support is not treated as a pass; score includes coverage flag.
            m['minor_gate_not_evaluable']=True
        hard=bool(legacy.get('real_hard_pass',False))
        if cfg.get('require_real_hard_gate',True): gate.append(('real_world_regression',hard))
        passed=all(status for _,status in gate)
        # Same objective for all individuals: TV. Tail and rare metrics are safety gates,
        # not independently weighted classification objectives.
        row={'candidate_id':cid,'representation':legacy.get('representation',''),
             'pass':passed,'real_hard_pass':hard,'balanced_mean_tv':m['balanced_mean_tv'],
             'balanced_p95_tv':m['balanced_p95_tv'],
             'balanced_off_target_mass':m['balanced_off_target_mass'],
             'balanced_dominant_accuracy':m['balanced_dominant_accuracy'],
             'minor_2_5_recall':m['minor_2_5_recall'],
             'minor_2_5_n':m['minor_2_5_n'],
             'n_reference':m['n_reference'],'n_pedigree':m['n_pedigree'],
             'scenario_count':m['scenario_count'],'legacy_selection_score':legacy.get('selection_score'),
             'gate_failures':','.join(name for name,status in gate if not status)}
        rows.append(row);reports.append({'candidate_id':cid,'metrics':m,'gate_results':dict(gate),'pass':passed})
    ranks=pd.DataFrame(rows).sort_values(['pass','balanced_mean_tv','balanced_p95_tv','balanced_off_target_mass','candidate_id'],ascending=[False,True,True,True,True])
    ranks.insert(0,'unified_rank',range(1,len(ranks)+1))
    ranks.to_csv(out/'unified_leaderboard.tsv',sep='\t',index=False)
    (out/'unified_diagnostics.json').write_text(json.dumps(reports,indent=2)+'\n')
    n=int(cfg.get('retain_finalists',5))
    eligible=ranks[ranks['pass']]
    if eligible.empty and cfg.get('require_any_pass',True):
        raise ValueError('No candidate passed unified error/rare-ancestry/regression gates; review unified_leaderboard.tsv rather than promoting a failing model')
    chosen=eligible.head(n) if not eligible.empty else ranks.head(n)
    (out/'finalists').mkdir(exist_ok=True)
    for cid in chosen['candidate_id']:
        shutil.copy2(dirs[cid][0]/'candidate.json',out/'finalists'/f'{cid}.json')
    chosen.to_csv(out/'finalists.tsv',sep='\t',index=False)
    return ranks,chosen

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--candidate-dirs',nargs='+',required=True)
    p.add_argument('--ped-predictions',nargs='+',required=True)
    p.add_argument('--config',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args()
    ranks,chosen=choose(a.candidate_dirs,a.ped_predictions,a.config,a.outdir)
    print(f'Unified selection ranked {len(ranks)} candidates; selected {len(chosen)} finalists')
if __name__=='__main__': main()
