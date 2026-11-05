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
            r=s.run(d/'model.joblib',d/'truth.tsv',[d/'features.tsv'],d/'pred.tsv',d/'m.json','c1','r1','1,2')
            self.assertEqual(r['aggregate']['n'],2)
            self.assertAlmostEqual(r['aggregate']['mean_tv'],0.0)
            self.assertEqual(r['production_domain_validated'],False)
            ranks=summ.summarize([d/'m.json'],d/'rank.tsv')
            self.assertEqual(ranks[0]['candidate_id'],'c1')
            self.assertIn('high_density_direct', (d/'rank.tsv').read_text())

if __name__=='__main__':unittest.main()
