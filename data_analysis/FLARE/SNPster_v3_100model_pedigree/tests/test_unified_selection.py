import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pedigree_benchmark'/'scripts'))
from select_unified_ancestry import SOURCES, calculate, choose, pedigree_frame

class UnifiedAncestrySelection(unittest.TestCase):
    def row(self,sample,scenario,family,truth,pred):
        d={'sample_id':sample,'scenario':scenario,'family_id':family}
        for s in SOURCES:
            d['truth_'+s]=truth.get(s,0.)
            d['pred_'+s]=pred.get(s,0.)
        return d
    def test_same_tv_for_pure_and_admixed(self):
        rows=[self.row('pure','reference_source','rp1',{'EUR':1},{'EUR':.9,'SSA':.1}),
              self.row('f1','f1','f1',{'EUR':.5,'SSA':.5},{'EUR':.6,'SSA':.4})]
        m=calculate(pd.DataFrame(rows),{})
        self.assertAlmostEqual(m['balanced_mean_tv'],.1)
        self.assertEqual(m['scenario_count'],2)
    def test_more_siblings_do_not_overweight_pedigree(self):
        a=self.row('a','f1','fam_a',{'EUR':.5,'SSA':.5},{'EUR':.7,'SSA':.3})
        b=self.row('b','f1','fam_b',{'EUR':.5,'SSA':.5},{'EUR':.5,'SSA':.5})
        rows=[a,b]+[dict(a,sample_id=f'x{i}') for i in range(100)]
        m=calculate(pd.DataFrame(rows),{})
        self.assertAlmostEqual(m['balanced_mean_tv'],.1)
        self.assertGreater(m['unstratified_mean_tv'],.19)
    def test_rare_ancestry_recall_separate_from_tv(self):
        rows=[self.row('minor','dilution_3pct','fam',{'EUR':.96875,'SEA':.03125},{'EUR':1})]
        m=calculate(pd.DataFrame(rows),{})
        self.assertAlmostEqual(m['balanced_mean_tv'],.03125)
        self.assertEqual(m['minor_2_5_n'],1)
        self.assertEqual(m['minor_2_5_recall'],0)
    def make_candidate(self,root,name,eur_prediction,real_pass=True):
        folder=root/f'candidate_{name}';folder.mkdir()
        (folder/'metrics.json').write_text(json.dumps({'candidate_id':name,'representation':'rep',
             'real_hard_pass':real_pass,'selection_score':.1}))
        (folder/'candidate.json').write_text(json.dumps({'candidate_id':name,'representation':'rep'}))
        rr=[]
        for i in range(8):
            row={'sample_id':f'ref{i}','role':'validation'}
            for s in SOURCES:
                row[s]=eur_prediction if s=='EUR' else (1-eur_prediction if s=='SSA' else 0)
                row['truth_'+s]=1 if s=='EUR' else 0
            rr.append(row)
        pd.DataFrame(rr).to_csv(folder/'predictions.tsv',sep='\t',index=False)
        ped=root/f'{name}.pedigree_predictions.tsv'
        p=[]
        for i in range(6):
            row={'sample_id':f'sim{i}','class':'f1','family_id':f'fam{i}',
                 'pedigree_design':'reduced','reference_disjoint':True,'founder_labels_validated':True,
                 'pred_EUR':.5,'pred_SSA':.5,'realized_EUR':.5,'realized_SSA':.5}
            for s in SOURCES:
                row.setdefault('pred_'+s,0.);row.setdefault('realized_'+s,0.)
            p.append(row)
        pd.DataFrame(p).to_csv(ped,sep='\t',index=False)
        return folder,ped
    def test_selection_uses_combined_tv_and_enforces_real_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            a,ap=self.make_candidate(root,'good',.99)
            b,bp=self.make_candidate(root,'bad',.75)
            cfg=root/'cfg.json';cfg.write_text(json.dumps({'retain_finalists':1,
                'min_reference_samples':8,'min_pedigree_samples':4,
                'required_scenarios':['reference_source','f1'],
                'max_gates':{'balanced_mean_tv':.3},'require_real_hard_gate':True}))
            ranks,chosen=choose([b,a],[bp,ap],cfg,root/'output')
            self.assertEqual(chosen.iloc[0]['candidate_id'],'good')
            self.assertEqual(chosen.iloc[0]['pedigree_design'],'reduced')
            self.assertEqual(len(ranks),2)
            self.assertTrue((root/'output'/'finalists'/'good.json').exists())
            b,bp=self.make_candidate(root,'regression',.999,False)
            ranks,chosen=choose([a,b],[ap,bp],cfg,root/'another')
            self.assertEqual(chosen.iloc[0]['candidate_id'],'good')
    def test_inconsistent_truth_errors(self):
        row=self.row('bad','f1','fam',{'EUR':.2,'SSA':.2},{'EUR':1})
        with self.assertRaises(ValueError): calculate(pd.DataFrame([row]),{})

    def test_historical_and_mixed_designs_cannot_be_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            a,ap=self.make_candidate(root,'reduced',.99)
            b,bp=self.make_candidate(root,'full',.99)
            full=pd.read_csv(bp,sep='\t')
            full['pedigree_design']='full'
            full.to_csv(bp,sep='\t',index=False)
            cfg=root/'cfg.json';cfg.write_text('{}')
            with self.assertRaisesRegex(ValueError,'Cannot compare'):
                choose([a,b],[ap,bp],cfg,root/'out')
            historical=pd.read_csv(ap,sep='\t').drop(columns=['reference_disjoint'])
            historical.to_csv(ap,sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'historical scores'):
                pedigree_frame(ap)

if __name__=='__main__': unittest.main()
