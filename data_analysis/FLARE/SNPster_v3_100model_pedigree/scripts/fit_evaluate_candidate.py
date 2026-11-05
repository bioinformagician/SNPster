#!/usr/bin/env python3
import argparse, inspect, json, math, os, re
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold, GroupKFold

SOURCES = ['EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA']
EPS = 1e-12


def simplex_row(v):
    v = np.asarray(v, float)
    u = np.sort(v)[::-1]
    cssv = np.cumsum(u) - 1
    ind = np.arange(1, len(v)+1)
    cond = u - cssv / ind > 0
    if not np.any(cond):
        return np.ones_like(v) / len(v)
    rho = ind[cond][-1]
    theta = cssv[rho-1] / rho
    w = np.maximum(v-theta, 0)
    s = w.sum()
    return w/s if s > 0 else np.ones_like(v)/len(v)


def simplex(X):
    X = np.atleast_2d(np.asarray(X, float))
    return np.vstack([simplex_row(x) for x in X])


def detect_sample_col(df):
    for c in ['sample_id','sample','SAMPLE','iid','IID']:
        if c in df.columns:
            return c
    return df.columns[0]


def load_truth(path, mode='auto'):
    df = pd.read_csv(path, sep='\t')
    sample_col = detect_sample_col(df)
    out = pd.DataFrame({'sample_id': df[sample_col].astype(str)})

    modes = [mode] if mode != 'auto' else ['exact','true','truth','prop','p','source_group']
    found = False
    for m in modes:
        if m == 'source_group' and 'source_group' in df.columns:
            for s in SOURCES:
                out[s] = (df['source_group'].astype(str) == s).astype(float)
            found = True; break
        prefix = {'exact':'', 'true':'true_', 'truth':'truth_', 'prop':'prop_', 'p':'p_'}.get(m)
        if prefix is not None and all(prefix+s in df.columns for s in SOURCES):
            for s in SOURCES: out[s] = pd.to_numeric(df[prefix+s], errors='raise')
            found = True; break
    if not found:
        if 'source_group' in df.columns:
            for s in SOURCES:
                out[s] = (df['source_group'].astype(str) == s).astype(float)
            found = True
    if not found:
        raise ValueError(f'Could not infer truth columns in {path}; expected {SOURCES}, true_*, truth_*, prop_*, p_* or source_group')

    Y = out[SOURCES].to_numpy(float)
    rs = Y.sum(axis=1, keepdims=True)
    if np.any(rs <= 0): raise ValueError(f'Zero-sum truth row in {path}')
    out[SOURCES] = Y / rs

    # Preserve useful stratification/group columns if present.
    extras = ['scenario','n_sources','template','platform','fine','source_group','family_id','cv_group','donor_group','mosaic_family']
    for c in extras:
        if c in df.columns and c not in out.columns:
            out[c] = df[c].values
    return out


def load_features(path):
    df = pd.read_csv(path, sep='\t')
    sample_col = detect_sample_col(df)
    df = df.rename(columns={sample_col:'sample_id'})
    df['sample_id'] = df['sample_id'].astype(str)
    feats = [c for c in df.columns if c.startswith('anc_') or c.startswith('f_') or re.fullmatch(r'f\d+', c or '')]
    if not feats:
        feats = [c for c in df.columns if c != 'sample_id' and pd.api.types.is_numeric_dtype(df[c])]
    if not feats: raise ValueError(f'No numeric feature columns found in {path}')
    return df[['sample_id']+feats], feats


def make_estimator(cfg, seed):
    t = cfg['type']
    p = dict(cfg.get('params', {}))
    n_jobs = int(os.environ.get('SKLEARN_N_JOBS', '1'))
    if t == 'extratrees':
        return ExtraTreesRegressor(random_state=seed, n_jobs=n_jobs, **p)
    if t == 'randomforest':
        return RandomForestRegressor(random_state=seed, n_jobs=n_jobs, **p)
    if t == 'ridge':
        return Ridge(**p)
    if t == 'histgb':
        base = HistGradientBoostingRegressor(random_state=seed, **p)
        return MultiOutputRegressor(base, n_jobs=n_jobs)
    if t == 'mlp':
        if 'hidden_layer_sizes' in p and isinstance(p['hidden_layer_sizes'], list):
            p['hidden_layer_sizes'] = tuple(p['hidden_layer_sizes'])
        p.setdefault('max_iter', 1200)
        p.setdefault('early_stopping', True)
        return MLPRegressor(random_state=seed, **p)
    if t == 'centroid_nnls':
        return {'kind':'centroid_nnls'}
    if t == 'frozen_joblib':
        return joblib.load(p['path'])
    raise ValueError(f'Unsupported decoder: {t}')


def fit_with_optional_weight(model, X, Y, w=None):
    if isinstance(model, dict):
        # Weighted source centroids using truth fractions as responsibilities.
        centroids = []
        for j in range(Y.shape[1]):
            r = Y[:,j].copy()
            if w is not None: r *= w
            if r.sum() <= 0: raise ValueError(f'No training mass for source index {j}')
            centroids.append(np.average(X, axis=0, weights=r))
        model['centroids'] = np.vstack(centroids)
        return model
    if w is not None:
        if 'sample_weight' not in inspect.signature(model.fit).parameters:
            raise ValueError(f'{type(model).__name__} does not support required sample weights')
        model.fit(X, Y, sample_weight=w)
    else:
        model.fit(X, Y)
    return model


def predict_model(model, X):
    if isinstance(model, dict) and model.get('kind') == 'centroid_nnls':
        A = model['centroids'].T
        rows = []
        for x in X:
            coef, _ = nnls(A, x)
            if coef.sum() <= 0: coef = np.ones(A.shape[1])
            rows.append(coef/coef.sum())
        return np.vstack(rows)
    return simplex(model.predict(X))


def group_vector(expr, row):
    return sum(float(row[s]) for s in expr.split('+'))


def evaluate_constraints(pred_row, constraints):
    details = {}
    passed = True
    for expr, rule in constraints.items():
        val = group_vector(expr, pred_row)
        ok = True
        if 'min' in rule: ok &= val >= float(rule['min'])
        if 'max' in rule: ok &= val <= float(rule['max'])
        if 'target' in rule and 'tol' in rule: ok &= abs(val-float(rule['target'])) <= float(rule['tol'])
        details[expr] = {'value':val, 'rule':rule, 'pass':bool(ok)}
        passed &= bool(ok)
    return passed, details


def metrics_from_df(df):
    Y = df[[f'truth_{s}' for s in SOURCES]].to_numpy(float)
    P = df[SOURCES].to_numpy(float)
    tv = 0.5*np.abs(P-Y).sum(axis=1)
    off = np.where(Y <= 1e-12, P, 0).sum(axis=1)
    dom = (P.argmax(axis=1) == Y.argmax(axis=1)).astype(float)
    out = {
        'n': int(len(df)),
        'mean_tv': float(tv.mean()),
        'median_tv': float(np.median(tv)),
        'p95_tv': float(np.quantile(tv, .95)),
        'max_tv': float(tv.max()),
        'dominant_accuracy': float(dom.mean()),
        'off_target_mass': float(off.mean()),
        'source_mae': {s: float(np.mean(np.abs(P[:,i]-Y[:,i]))) for i,s in enumerate(SOURCES)},
        'source_bias': {s: float(np.mean(P[:,i]-Y[:,i])) for i,s in enumerate(SOURCES)},
    }
    return out



def stratified_metrics(df):
    out = {}
    for col in ['scenario','n_sources','template','platform','source_group','fine']:
        if col not in df.columns:
            continue
        groups = {}
        for key, g in df.groupby(col, dropna=False):
            groups[str(key)] = metrics_from_df(g)
        out[col] = groups
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate', required=True)
    ap.add_argument('--datasets', required=True)
    ap.add_argument('--features-dir', required=True)
    ap.add_argument('--real-benchmarks', required=True)
    ap.add_argument('--success-criteria', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--save-model', action='store_true', help='Persist model.joblib (recommended only for stage 3/finalists)')
    args = ap.parse_args()

    cfg = json.loads(Path(args.candidate).read_text())
    rep = cfg['representation']; seed = int(cfg.get('seed', 20261005))
    ds = pd.read_csv(args.datasets, sep='\t', dtype=str).fillna('')
    ds = ds[ds['enabled'].str.lower().ne('false')].copy()
    feat_dir = Path(args.features_dir)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    data = {}
    feature_cols = None
    for _, r in ds.iterrows():
        did = r['dataset_id']
        fp = feat_dir / f'{rep}__{did}.features.tsv'
        if not fp.exists():
            raise FileNotFoundError(f'{rep}: missing features for enabled dataset {did}: {fp}')
        fdf, fcols = load_features(fp)
        if feature_cols is None: feature_cols = fcols
        if fcols != feature_cols:
            raise ValueError(f'Feature columns differ for {did}: {fcols} vs {feature_cols}')
        truth = None
        if r.get('truth_path','') and r['truth_path'] != '.':
            truth = load_truth(r['truth_path'], r.get('truth_mode','auto') or 'auto')
            merged = fdf.merge(truth, on='sample_id', how='inner', validate='one_to_one')
            if len(merged) != len(fdf):
                missing = set(fdf.sample_id)-set(merged.sample_id)
                raise ValueError(f'{did}: {len(missing)} feature samples missing truth, e.g. {list(missing)[:5]}')
        else:
            merged = fdf.copy()
        data[did] = (r, merged)

    train_parts=[]; y_parts=[]; w_parts=[]; group_parts=[]
    weight_map = {str(k):float(v) for k,v in cfg.get('training_weights',{}).items()}
    for did,(r,d) in data.items():
        if r['role'] != 'train_decoder': continue
        train_parts.append(d[feature_cols].to_numpy(float))
        y_parts.append(d[SOURCES].to_numpy(float))
        tag = r.get('weight_tag','') or 'default'
        w_parts.append(np.full(len(d), weight_map.get(tag, weight_map.get('default',1.0))))
        gcol = next((c for c in ['cv_group','family_id','mosaic_family','donor_group'] if c in d.columns), None)
        if gcol: group_parts.extend([f'{did}:{x}' for x in d[gcol].astype(str)])
        else: group_parts.extend([f'{did}:{x}' for x in d['sample_id'].astype(str)])
    if not train_parts: raise ValueError('No role=train_decoder datasets available')
    X=np.vstack(train_parts); Y=np.vstack(y_parts); W=np.concatenate(w_parts); groups=np.asarray(group_parts)

    # OOF raw predictions for calibration; exact frozen artifacts can bypass retraining.
    cal_cfg = cfg.get('calibrator', {'type':'none','params':{}})
    calibrator = None
    frozen_decoder = cfg['decoder']['type'] == 'frozen_joblib'

    if frozen_decoder:
        decoder = make_estimator(cfg['decoder'], seed)
        if cal_cfg['type'] == 'frozen_joblib':
            calibrator = joblib.load(cal_cfg['params']['path'])
        elif cal_cfg['type'] != 'none':
            raise ValueError('Frozen decoder candidates must use calibrator=none or frozen_joblib')
    else:
        if cal_cfg['type'] != 'none':
            folds = int(cfg.get('calibration_folds',3))
            unique_groups = np.unique(groups)
            if len(unique_groups) >= folds:
                splitter = GroupKFold(n_splits=folds).split(X, Y, groups)
            else:
                splitter = KFold(n_splits=folds, shuffle=True, random_state=seed).split(X)
            oof = np.zeros_like(Y, dtype=float)
            for fold,(tr,va) in enumerate(splitter):
                m = make_estimator(cfg['decoder'], seed+fold+1)
                m = fit_with_optional_weight(m, X[tr], Y[tr], W[tr])
                oof[va] = predict_model(m, X[va])
            if cal_cfg['type'] == 'ridge':
                calibrator = Ridge(**cal_cfg.get('params',{})).fit(oof, Y, sample_weight=W)
            else:
                raise ValueError(f"Unsupported calibrator: {cal_cfg['type']}")

        decoder = make_estimator(cfg['decoder'], seed)
        decoder = fit_with_optional_weight(decoder, X, Y, W)

    prediction_rows=[]; dataset_metrics={}; role_frames={}
    for did,(r,d) in data.items():
        if r['role'] == 'train_decoder':
            continue
        raw = predict_model(decoder, d[feature_cols].to_numpy(float))
        pred = simplex(calibrator.predict(raw)) if calibrator is not None else raw
        pdf = pd.DataFrame(pred, columns=SOURCES)
        pdf.insert(0,'sample_id', d['sample_id'].astype(str).values)
        pdf.insert(0,'role', r['role'])
        pdf.insert(0,'dataset_id', did)
        for i,s in enumerate(SOURCES): pdf[f'raw_{s}'] = raw[:,i]
        if all(s in d.columns for s in SOURCES):
            for s in SOURCES: pdf[f'truth_{s}'] = d[s].to_numpy(float)
            # Preserve grouping columns for downstream stratified metrics.
            for c in ['scenario','n_sources','template','platform','fine','source_group']:
                if c in d.columns: pdf[c] = d[c].values
            dataset_metrics[did] = metrics_from_df(pdf)
            dataset_metrics[did]['stratified'] = stratified_metrics(pdf)
            role_frames.setdefault(r['role'], []).append(pdf)
        prediction_rows.append(pdf)

    role_metrics={}
    for role,frames in role_frames.items():
        role_metrics[role] = metrics_from_df(pd.concat(frames, ignore_index=True))

    predictions = pd.concat(prediction_rows, ignore_index=True) if prediction_rows else pd.DataFrame()

    # Real-world regression benchmarks (Lonni is intentionally a regression benchmark, not an unbiased final test).
    real_results=[]; real_hard_pass=True
    rb = pd.read_csv(args.real_benchmarks, sep='\t', dtype=str).fillna('')
    for _,b in rb.iterrows():
        did=b['dataset_id']
        sub=predictions[predictions['dataset_id']==did]
        if sub.empty:
            real_results.append({'benchmark_id':b['benchmark_id'],'dataset_id':did,'pass':False,'error':'no prediction'})
            if b.get('hard_gate','').lower()=='true': real_hard_pass=False
            continue
        row=sub.iloc[0]
        constraints=json.loads(b['constraints_json'])
        ok,details=evaluate_constraints(row,constraints)
        rec={'benchmark_id':b['benchmark_id'],'dataset_id':did,'pass':bool(ok),'details':details,
             'prediction':{s:float(row[s]) for s in SOURCES}, 'hard_gate':b.get('hard_gate','').lower()=='true'}
        real_results.append(rec)
        if rec['hard_gate'] and not ok: real_hard_pass=False

    criteria=json.loads(Path(args.success_criteria).read_text())
    sel_role=criteria.get('selection_role','validation')
    sel=role_metrics.get(sel_role,{})
    if not sel:
        raise ValueError(f'{rep}: no truth data for required selection role={sel_role}')
    gates=criteria.get('gates',{})
    gate_details={}; overall_pass=True
    def gate(name, val, op, thr):
        nonlocal overall_pass
        ok = (val <= thr) if op=='max' else (val >= thr)
        gate_details[name]={'value':float(val),'threshold':float(thr),'pass':bool(ok)}
        overall_pass &= bool(ok)
    if sel:
        if 'mean_tv_max' in gates: gate('mean_tv',sel['mean_tv'],'max',gates['mean_tv_max'])
        if 'p95_tv_max' in gates: gate('p95_tv',sel['p95_tv'],'max',gates['p95_tv_max'])
        if 'off_target_mass_max' in gates: gate('off_target_mass',sel['off_target_mass'],'max',gates['off_target_mass_max'])
        if 'dominant_accuracy_min' in gates: gate('dominant_accuracy',sel['dominant_accuracy'],'min',gates['dominant_accuracy_min'])
    if criteria.get('require_real_hard_gates', True): overall_pass &= real_hard_pass

    sw=criteria.get('score_weights',{'mean_tv':0.5,'p95_tv':0.3,'off_target_mass':0.2})
    score=0.0
    if sel:
        score += float(sw.get('mean_tv',0))*sel['mean_tv']
        score += float(sw.get('p95_tv',0))*sel['p95_tv']
        score += float(sw.get('off_target_mass',0))*sel['off_target_mass']
        score += float(sw.get('dominant_error',0))*(1-sel['dominant_accuracy'])
    if not real_hard_pass: score += float(criteria.get('real_gate_failure_penalty',1.0))

    metrics={
        'candidate_id':cfg['candidate_id'], 'config_sha256':cfg.get('config_sha256'), 'representation':rep,
        'candidate':cfg, 'dataset_metrics':dataset_metrics, 'role_metrics':role_metrics,
        'real_world_benchmarks':real_results, 'gate_details':gate_details,
        'real_hard_pass':bool(real_hard_pass), 'pass':bool(overall_pass), 'selection_score':float(score),
        'feature_columns':feature_cols,
    }
    (outdir/'metrics.json').write_text(json.dumps(metrics, indent=2, sort_keys=True)+'\n')
    predictions.to_csv(outdir/'predictions.tsv', sep='\t', index=False)
    if args.save_model:
        joblib.dump({'decoder':decoder,'calibrator':calibrator,'feature_columns':feature_cols,'sources':SOURCES,'candidate':cfg}, outdir/'model.joblib')
    (outdir/'candidate.json').write_text(json.dumps(cfg, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'candidate_id':cfg['candidate_id'],'pass':overall_pass,'score':score}))

if __name__ == '__main__':
    main()
