import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.dummy import DummyRegressor

ROOT=Path(__file__).resolve().parents[1]

def mod(filename):
    p=ROOT/'pedigree_benchmark'/'scripts'/filename
    spec=importlib.util.spec_from_file_location(filename[:-3],p)
    m=importlib.util.module_from_spec(spec)
    sys.modules[filename[:-3]]=m
    spec.loader.exec_module(m)
    return m

s=mod('evaluate_finalist.py')
summ=mod('summarize_metrics.py')

class TestScoring(unittest.TestCase):
    def test_saved_frozen_candidate_and_realized_truth(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td)
            y=[.5,0,0,.5,0,0,0,0]
            estimator=DummyRegressor(strategy='constant',constant=y)
            estimator.fit(np.ones((3,2)),np.zeros((3,8)))
            joblib.dump({'decoder':estimator,'calibrator':None,'feature_columns':['anc_1','anc_2'],
                         'sources':s.SOURCES,'candidate':{'candidate_id':'c1','representation':'r1'}},d/'model.joblib')
            f=pd.DataFrame({'sample_id':['p1','p2'],'anc_1':[.2,.6],'anc_2':[.8,.4]})
            f.to_csv(d/'features.tsv',sep='\t',index=False)
            t=pd.DataFrame([{'sample_id':k,'family_id':'fam1','class':'f1','generation':2,
                             'is_founder':0,'minimum_nonzero_ancestry':.5,
                             **{f'expected_{src}':frac for src,frac in zip(s.SOURCES,y)},
                             **{f'realized_{src}':frac for src,frac in zip(s.SOURCES,y)}} for k in ['p1','p2']])
            t.to_csv(d/'truth.tsv',sep='\t',index=False)
            audit={'strict_reference_disjoint':True,'labels_validated':True,
                   'mosaic_training_exclusion_checked':True,'flare_reference_overlap_count':0}
            (d/'audit.json').write_text(json.dumps(audit))
            (d/'manifest.json').write_text(json.dumps({'design':'reduced'}))
            r=s.run(d/'model.joblib',d/'truth.tsv',[d/'features.tsv'],d/'pred.tsv',d/'m.json','c1','r1','1,2',
                    d/'audit.json',d/'manifest.json')
            self.assertEqual(r['aggregate']['n'],2)
            self.assertAlmostEqual(r['aggregate']['mean_tv'],0.0)
            self.assertEqual(r['production_domain_validated'],False)
            self.assertEqual(r['pedigree_design'],'reduced')
            self.assertTrue(r['reference_disjoint'])
            ranks=summ.summarize([d/'m.json'],d/'rank.tsv')
            self.assertEqual(ranks[0]['candidate_id'],'c1')
            self.assertIn('high_density_direct', (d/'rank.tsv').read_text())
            audit['flare_reference_overlap_count']=1
            (d/'audit.json').write_text(json.dumps(audit))
            with self.assertRaisesRegex(ValueError,'reference-disjoint'):
                s.run(d/'model.joblib',d/'truth.tsv',[d/'features.tsv'],d/'pred.tsv',d/'m.json','c1','r1','1,2',
                      d/'audit.json',d/'manifest.json')

if __name__=='__main__':unittest.main()
