#!/usr/bin/env python3
import argparse,csv,json,sys
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument('--params',required=True); ap.add_argument('--search-space',required=True); a=ap.parse_args(); params=json.loads(Path(a.params).read_text()); errors=[]
for key in ['flare_jar','flare_ref','flare_panel','source_train_vcf','genetic_map']:
 if not Path(params[key]).exists(): errors.append(f'Missing {key}: {params[key]}')
for table_key in ['fixed_models','dev_datasets','diagnostic_datasets','full_datasets','real_benchmarks']:
 if not Path(params[table_key]).exists(): errors.append(f'Missing config table {table_key}: {params[table_key]}')
if Path(params['fixed_models']).exists():
 with open(params['fixed_models']) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r.get('enabled','true').lower()!='false' and not Path(r['model_path']).exists(): errors.append(f"Missing fixed model {r['rep_id']}: {r['model_path']}")
for table_key in ['dev_datasets','diagnostic_datasets']:
 p=Path(params[table_key])
 if p.exists():
  with open(p) as f:
   for r in csv.DictReader(f,delimiter='\t'):
    if r.get('enabled','true').lower()=='false': continue
    if not Path(r['input_path']).exists(): errors.append(f"Missing dataset input {r['dataset_id']}: {r['input_path']}")
    tp=r.get('truth_path','').strip()
    if tp and tp!='.' and not Path(tp).exists(): errors.append(f"Missing dataset truth {r['dataset_id']}: {tp}")
if errors:
 print('PREFLIGHT FAILED',file=sys.stderr); [print(' - '+e,file=sys.stderr) for e in errors]; sys.exit(2)
print('Preflight OK')
