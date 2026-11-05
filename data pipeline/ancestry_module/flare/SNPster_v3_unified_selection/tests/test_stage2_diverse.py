import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/'pedigree_benchmark'/'scripts'/'select_stage2_diverse.py'
spec=importlib.util.spec_from_file_location('diverse_stage2',SCRIPT)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

class DiversityTests(unittest.TestCase):
    def test_representation_preservation(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            cfg=root/'space.json';cfg.write_text(json.dumps({'stage2':{'retain_top_candidates':3,'fixed_candidates':[]}}))
            for cid,rep,score in [('a1','A',.01),('a2','A',.02),('a3','A',.03),('b','B',.4),('c','C',.5)]:
                p=root/f'candidate_{cid}';p.mkdir()
                (p/'metrics.json').write_text(json.dumps({'candidate_id':cid,'representation':rep,'pass':True,
                     'real_hard_pass':True,'selection_score':score}))
                (p/'candidate.json').write_text(json.dumps({'candidate_id':cid,'representation':rep}))
            x=mod.select(root,cfg,root/'out',root/'selected.tsv')
            self.assertEqual(set(x.representation),{'A','B','C'})
            self.assertEqual(len(x),3)

if __name__=='__main__': unittest.main()
