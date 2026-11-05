#!/usr/bin/env python3
import argparse,json
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument('--rep-id',required=True); ap.add_argument('--search-space',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
s=json.loads(Path(a.search_space).read_text()); st=s['stage1']
c={
 'candidate_id':f"screen__{a.rep_id}", 'representation':a.rep_id, 'seed':s.get('seed',20261005),
 'decoder':st['screen_decoder'],'calibrator':st['screen_calibrator'],
 'training_weights':st['training_weights'],'calibration_folds':st.get('calibration_folds',3)
}
Path(a.output).write_text(json.dumps(c,indent=2,sort_keys=True)+'\n')
