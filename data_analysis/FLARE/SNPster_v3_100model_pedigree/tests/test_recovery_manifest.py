import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def load_script(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

class ManifestRegressionTests(unittest.TestCase):
    def make_case(self, root, count=4):
        manifest = root/'stage2_result_dirs.txt'
        folders=[]
        for i in range(count):
            folder=root/f'candidate_{i}'
            folder.mkdir()
            meta={'candidate_id':f'cid{i}', 'representation':f'rep{i%2}',
                  'pass':True, 'real_hard_pass':True, 'selection_score':i*.01}
            (folder/'metrics.json').write_text(json.dumps(meta))
            (folder/'candidate.json').write_text(json.dumps(meta))
            folders.append(folder)
        manifest.write_text(''.join(f'{x}\n' for x in folders))
        search=root/'space.json'
        search.write_text(json.dumps({'stage2':{'retain_top_candidates':4,'fixed_candidates':[]}}))
        return manifest,search,folders

    def test_diverse_manifest_happy_path(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            manifest,search,folders=self.make_case(root)
            mod=load_script(ROOT/'pedigree_benchmark/scripts/select_stage2_diverse.py','diverse_manifest')
            result=mod.select(None,search,root/'out',root/'selection.tsv',manifest)
            self.assertEqual(len(result),4)
            self.assertEqual(set(result.representation),{'rep0','rep1'})

    def test_manifest_missing_candidate_reports_path_and_line(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            manifest,search,folders=self.make_case(root)
            (folders[1]/'metrics.json').unlink()
            mod=load_script(ROOT/'pedigree_benchmark/scripts/select_stage2_diverse.py','diverse_bad')
            with self.assertRaisesRegex(ValueError,r'Manifest line 2.*metrics.json'):
                mod.select(None,search,root/'out',root/'selection.tsv',manifest)

    def test_legacy_selector_supports_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            manifest,search,folders=self.make_case(root)
            completed=subprocess.run([sys.executable,str(ROOT/'scripts/select_candidates.py'),
                '--manifest',str(manifest),'--search-space',str(search),
                '--outdir',str(root/'out'),'--table',str(root/'selected.tsv')],
                capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)
            self.assertTrue((root/'out'/'cid0.json').exists())

    def test_one_hundred_configs_exact_single_bayesian(self):
        search=json.loads((ROOT/'config/search_space_100.json').read_text())
        n=[]
        for strategy in search['flare_model_strategies']:
            count=1
            for values in strategy.get('grid',{}).values(): count*=len(values)
            n.append((strategy['family'],count))
        self.assertEqual(sum(c for _,c in n),100)
        self.assertEqual(sum(c for family,c in n if family=='bayesian_gaussian_mixture'),1)

    def test_recovery_skips_training_code_in_branch(self):
        nf=(ROOT/'main.nf').read_text()
        recovery=nf.split('if (params.resume_stage2_manifest) {',1)[1].split('} else {',1)[0]
        self.assertNotIn('BUILD_INITIAL_MODEL(',recovery)
        self.assertNotIn('FIT_STAGE2_CANDIDATE(',recovery)
        self.assertIn('SELECT_STAGE2_CANDIDATES(',recovery)

if __name__ == '__main__': unittest.main()
