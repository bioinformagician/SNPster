#!/usr/bin/env python3
"""Unit tests; no Ped-sim or Nextflow binary needed."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def module(name):
    p=ROOT/'scripts'/f'{name}.py'
    spec=importlib.util.spec_from_file_location(name,p)
    mod=importlib.util.module_from_spec(spec)
    import sys
    sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod

plan=module('make_plans')
truth=module('realized_ancestry')
agg=module('aggregate_truth')
verify=module('verify_founder_assignments')

class TestPedigrees(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.base=Path(self.tmp.name)
        self.meta=self.base/'source_validation_metadata.tsv'
        with self.meta.open('w') as f:
            f.write('sample\tfine\tsource_group\n')
            for src in plan.SOURCES:
                for j in range(100): f.write(f'{src}_donor_{j}\tX\t{src}\n')
        self.excl=self.base/'train_metadata.tsv'
        self.excl.write_text('sample\n'+''.join(f'{s}_donor_0\n' for s in plan.SOURCES))

    def tearDown(self): self.tmp.cleanup()

    def test_plan_mixes_sexes_donors_and_lineages(self):
        out=self.base/'plans'
        plan.generate(self.meta,[self.excl],out,18,4,1234,False)
        batch_dirs=sorted(out.glob('batch_*'))
        self.assertEqual(len(batch_dirs),5)
        roots=[];assigned=[]
        for directory in batch_dirs:
            info=json.loads((directory/'pedigrees.json').read_text())
            donors=[v['vcf_sample'] for v in info['founders'].values()]
            assigned+=donors
            self.assertEqual(len(set(donors)),len(donors))
            self.assertFalse(any(x.endswith('_0') for x in donors))
            people={n['sample_id']:n for n in info['nodes']}
            for n in info['nodes']:
                self.assertAlmostEqual(sum(n['expected'].values()),1.)
                if not n['is_founder']:
                    self.assertEqual(people[n['father_id']]['generation']+1,n['generation'])
                    self.assertEqual(people[n['mother_id']]['generation']+1,n['generation'])
                    self.assertNotEqual(n['father_id'],n['mother_id'])
            roots += [n for n in info['nodes'] if not any(n['sample_id'] in (q['father_id'],q['mother_id']) for q in info['nodes'])]
            d=(directory/'pedigrees.def').read_text()
            self.assertIn('def fam',d)
            self.assertIn('sM',d);self.assertIn('sF',d)
        self.assertEqual(len(roots),18)
        self.assertEqual(len(assigned),len(set(assigned)))
        f1=next(x for x in roots if x['class']=='f1')
        self.assertEqual(f1['expected'],{'EUR':0.5,'SEA':0.5})
        dilution=next(x for x in roots if x['class']=='dilution_3pct')
        self.assertAlmostEqual(min(dilution['expected'].values()),0.03125)
        self.assertEqual(dilution['generation'],6)
        repeated=next(x for x in roots if x['class']=='four_way_repeated_founders')
        self.assertEqual(repeated['expected'],{'EUR':0.25,'SSA':0.75})
        self.assertTrue(any(x['class']=='eight_way_admixed_x_admixed' for x in roots))
        rare=next(x for x in roots if x['class']=='rare_3pct_multisource')
        self.assertAlmostEqual(min(rare['expected'].values()),0.03125)
        serial=next(x for x in roots if x['class']=='serial_admixed_mating_3pct')
        self.assertEqual(serial['generation'],6)
        self.assertAlmostEqual(min(serial['expected'].values()),0.03125)

    def test_bp_truth_and_aggregation(self):
        out=self.base/'plans'
        plan.generate(self.meta,[self.excl],out,1,1,23,False)
        batch=out/'batch_00001'
        info=json.loads((batch/'pedigrees.json').read_text())
        founders=list(info['founders'])
        # Controlled mock BP: 2 founder samples, one admixed F1 (no recombination needed).
        f1,f2=founders
        kid=next(r['sample_id'] for r in info['nodes'] if not r['is_founder'])
        bp=self.base/'sim.bp'
        bp.write_text(
            f'{f1} s0 h0 1|100 0:1099\n'
            f'{f1} s0 h1 1|100 1:1099\n'
            f'{f2} s1 h0 1|100 2:1099\n'
            f'{f2} s1 h1 1|100 3:1099\n'
            f'{kid} s0 h0 1|100 0:1099\n'
            f'{kid} s0 h1 1|100 2:1099\n'
        )
        result=self.base/'truth1.tsv'
        truth.ancestry_for_batch(batch/'pedigrees.json',bp,'1',result)
        with result.open() as fh:
            rec={x['sample_id']:x for x in csv.DictReader(fh,delimiter='\t')}
        s1=info['founders'][f1]['source'];s2=info['founders'][f2]['source']
        self.assertAlmostEqual(float(rec[kid][f'realized_{s1}']),0.5)
        self.assertAlmostEqual(float(rec[kid][f'realized_{s2}']),0.5)
        # Chromosome 2 is a duplicate mock with same source ratios.
        bp2=self.base/'sim2.bp'
        bp2.write_text(bp.read_text().replace('1|100','2|100'))
        result2=self.base/'truth2.tsv'
        truth.ancestry_for_batch(batch/'pedigrees.json',bp2,'2',result2)
        final=self.base/'final.tsv'
        self.assertEqual(agg.run([result,result2],final,['1','2']),3)
        with final.open() as fh:
            row={x['sample_id']:x for x in csv.DictReader(fh,delimiter='\t')}[kid]
        self.assertAlmostEqual(float(row['tv_expected_vs_realized']),0.)
        self.assertEqual(int(row['bp_weight']),4000)
        with self.assertRaises(ValueError): agg.run([result],self.base/'missing.tsv',['1','2'])

    def test_founder_assignment_verification(self):
        expected=self.base/'a.ids'; actual=self.base/'b.ids'
        expected.write_text('fam1_g1-b1-i1\tEUR_donor_1\n')
        actual.write_text('fam1_g1-b1-i1\tEUR_donor_1\n')
        verify.run(expected,actual)
        actual.write_text('fam1_g1-b1-i1\tEUR_donor_2\n')
        with self.assertRaisesRegex(ValueError,'mismatch'):
            verify.run(expected,actual)

    def test_missing_source_is_fatal(self):
        meta=self.base/'partial.tsv'; meta.write_text('sample\tsource_group\nX\tEUR\n')
        with self.assertRaisesRegex(ValueError,'No eligible'):
            plan.generate(meta,[],self.base/'bad',1,1,1,True)

    def test_reduced_design_preserves_classes_and_founder_uniqueness(self):
        meta=self.base/'small.tsv'
        with meta.open('w') as fh:
            fh.write('sample\tfine\tsource_group\n')
            for source in plan.SOURCES:
                count=2 if source=='OCE' else 4 if source in ('SEA','NAM') else 8
                for i in range(count):
                    fh.write(f'{source}_{i}\t{source}_fine\t{source}\n')
        with self.assertRaisesRegex(ValueError,'Insufficient'):
            plan.generate(meta,[],self.base/'full',252,10,20261008,True,'full')
        for dirname in ('reduced1','reduced2'):
            out=self.base/dirname
            plan.generate(meta,[],out,252,10,20261008,True,'reduced')
            manifest=json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['design'],'reduced')
            self.assertEqual(set(manifest['classes']),set(plan.CLASSES))
            self.assertTrue(all(n==28 for n in manifest['classes'].values()))
            self.assertTrue(manifest['rejected_infeasible_proposals'])
            for path in out.glob('batch_*/pedigrees.json'):
                info=json.loads(path.read_text())
                families={}
                for node in info['nodes']:
                    if node['is_founder']:
                        families.setdefault(node['family_id'],[]).append(info['founders'][node['sample_id']])
                for founders in families.values():
                    donors=[f['vcf_sample'] for f in founders]
                    self.assertEqual(len(donors),len(set(donors)))
                    self.assertLessEqual(sum(f['source']=='OCE' for f in founders),2)
        self.assertEqual((self.base/'reduced1/expected_ancestry.tsv').read_text(),
                         (self.base/'reduced2/expected_ancestry.tsv').read_text())

if __name__=='__main__': unittest.main()
