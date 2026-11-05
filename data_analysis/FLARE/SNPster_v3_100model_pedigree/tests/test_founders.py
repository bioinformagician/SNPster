import sys, unittest, tempfile, csv, json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pedigree_benchmark'/'scripts'))
from prepare_founders import prepare, SOURCES, mosaic_donor_ids, check_existing_output

class FounderTests(unittest.TestCase):
    def test_additional_founders_are_validated_and_audited(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            val = d/'val.tsv'
            train = d/'train.tsv'
            extra = d/'extra.tsv'
            out = d/'eligible.tsv'
            labels = d/'labels.tsv'
            labels.write_text('sample\tfine\tsource_group\n'
                              'OCE_1\tPapuan\tOCE\nOCE_2\tPapuan\tOCE\n'
                              'OCE_3\tPapuan\tOCE\nTRAIN_1\tFrench\tEUR\n')
            val.write_text('sample\tfine\tsource_group\nOCE_1\tPapuan\tOCE\nOCE_3\tPapuan\tOCE\n')
            train.write_text('sample\tfine\tsource_group\nTRAIN_1\tFrench\tEUR\n')
            extra.write_text('sample\tfine\tsource_group\nOCE_2\tPapuan\tOCE\n')
            result = prepare(val, train, {'OCE_1', 'OCE_2', 'OCE_3'}, out, 2,
                             required_sources=('OCE',), reference_samples={'OCE_2'},
                             additional_metadata=extra, label_metadata=labels)
            self.assertEqual(result['selected_per_source']['OCE'], 2)
            self.assertEqual(result['additional_reference_overlap_ids'], ['OCE_2'])
            self.assertEqual(result['excluded_reference_overlap_ids'], ['OCE_2'])
            self.assertEqual(result['flare_reference_overlap_count'], 0)
            self.assertNotIn('OCE_2', out.read_text())
            with self.assertRaisesRegex(ValueError, 'absent'):
                prepare(val, train, {'OCE_1'}, out, additional_metadata=extra,
                        reference_samples=set(), label_metadata=labels)
            train.write_text('sample\tfine\tsource_group\nOCE_2\tPapuan\tOCE\n')
            with self.assertRaisesRegex(ValueError, 'training donor'):
                prepare(val, train, {'OCE_1', 'OCE_2'}, out, additional_metadata=extra,
                        reference_samples=set(), label_metadata=labels)

    def test_intersect_heldout_and_exclude_train(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            val=d/'val.tsv'; train=d/'train.tsv'; out=d/'eligible.tsv'
            with val.open('w') as f:
                f.write('sample\tfine\tsource_group\n')
                for source in SOURCES:
                    for i in range(3): f.write(f'{source}_{i}\tF{i}\t{source}\n')
                f.write('EUR_3\tF3\tEUR\n')
                f.write('AFFINITY_0\tF\tCAS\n')
            # Use source-specific fine labels so the ontology is unambiguous.
            text = val.read_text()
            for source in SOURCES:
                for i in range(4):
                    text = text.replace(f'{source}_{i}\tF{i}', f'{source}_{i}\t{source}_F{i}')
            val.write_text(text)
            labels=d/'labels.tsv'; labels.write_text(text)
            train.write_text('sample\tfine\tsource_group\nEUR_0\tEUR_F0\tEUR\n')
            vcf={f'{s}_{i}' for s in SOURCES for i in range(3)} | {'AFFINITY_0','EUR_3'}
            res=prepare(val,train,vcf,out,2,reference_samples={'EUR_1'},label_metadata=labels)
            self.assertEqual(res['flare_reference_overlap_count'],0)
            self.assertEqual(res['excluded_reference_overlap_ids'],['EUR_1'])
            self.assertEqual(res['selected_donors'],23)
            with out.open() as fh:
                ids={r['sample'] for r in csv.DictReader(fh,delimiter='\t')}
            self.assertNotIn('EUR_0',ids)
            self.assertNotIn('EUR_1',ids)
            self.assertNotIn('AFFINITY_0',ids)
            with self.assertRaisesRegex(ValueError,'Insufficient'):
                prepare(val,train,vcf,out,3,reference_samples={'EUR_1'},label_metadata=labels)

    def test_label_mismatch_and_required_audits_are_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            val=d/'val.tsv'; train=d/'train.tsv'; labels=d/'labels.tsv'
            val.write_text('sample\tfine\tsource_group\nDai\t1KG_CDX\tSEA\n')
            train.write_text('sample\tfine\tsource_group\nTRAIN\t1KG_CDX\tEAS\n')
            labels.write_text('sample\tfine\tsource_group\nDai\t1KG_CDX\tEAS\nTRAIN\t1KG_CDX\tEAS\n')
            kwargs=dict(required_sources=('EAS',),reference_samples=set(),label_metadata=labels)
            with self.assertRaisesRegex(ValueError,'label mismatch.*Dai'):
                prepare(val,train,{'Dai'},d/'out',1,**kwargs)
            val.write_text('sample\tfine\tsource_group\nDai\t1KG_CDX\tEAS\n')
            with self.assertRaisesRegex(ValueError,'cannot be disabled'):
                prepare(val,train,{'Dai'},d/'out',1,strict_reference_disjoint=False,**kwargs)
            with self.assertRaisesRegex(ValueError,'Reference sample IDs'):
                prepare(val,train,{'Dai'},d/'out',1,label_metadata=labels)
            with self.assertRaisesRegex(ValueError,'Authoritative'):
                prepare(val,train,{'Dai'},d/'out',1,reference_samples=set())
            with self.assertRaisesRegex(ValueError,'fine label mismatch'):
                prepare(val,train,{'Dai'},d/'out',1,required_sources=('EAS',),
                        reference_samples={'TRAIN'}, reference_labels={'TRAIN':'1KG_KHV'},
                        label_metadata=labels)

    def test_mosaic_donors_are_excluded(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            val=d/'val.tsv'; train=d/'train.tsv'; labels=d/'labels.tsv'; mosaic=d/'mosaic.tsv'
            val.write_text('sample\tfine\tsource_group\nA\tPapuan\tOCE\nB\tPapuan\tOCE\n')
            train.write_text('sample\tfine\tsource_group\nTRAIN\tFrench\tEUR\n')
            labels.write_text(val.read_text()+train.read_text().split('\n',1)[1])
            mosaic.write_text('donor_weights_json\n{"A": 1.0}\n')
            result=prepare(val,train,{'A','B'},d/'out',1,required_sources=('OCE',),
                           reference_samples=set(),label_metadata=labels,mosaic_metadata=mosaic)
            self.assertEqual(result['selected_donors'],1)
            self.assertEqual(result['excluded_mosaic_donor_count'],1)
            self.assertTrue(result['mosaic_training_exclusion_checked'])
            self.assertNotIn('A\t', (d/'out').read_text())
            mosaic.write_text('donor_weights_json\n{"A": -1.0}\n')
            with self.assertRaisesRegex(ValueError,'invalid mosaic donor weights'):
                mosaic_donor_ids(mosaic)

    def test_output_protection_and_concurrent_map_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            (out/'pedigree/maps').mkdir(parents=True)
            (out/'pedigree/maps/map').write_text('map process finished first')
            check_existing_output(out,'reduced')
            (out/'pedigree/evaluations').mkdir()
            (out/'pedigree/evaluations/old.tsv').write_text('old predictions')
            with self.assertRaisesRegex(ValueError,'no founder audit'):
                check_existing_output(out,'reduced')
            inputs=out/'pedigree/inputs';inputs.mkdir()
            audit_path=inputs/'eligible_heldout_founders.tsv.summary.json'
            audit_path.write_text(json.dumps({'strict_reference_disjoint':False,
                                             'flare_reference_overlap_count':5}))
            with self.assertRaisesRegex(ValueError,'Historical'):
                check_existing_output(out,'reduced')
            audit_path.write_text(json.dumps({'strict_reference_disjoint':True,
                'labels_validated':True,'mosaic_training_exclusion_checked':True,
                'flare_reference_overlap_count':0}))
            plans=out/'pedigree/plans/plans';plans.mkdir(parents=True)
            (plans/'manifest.json').write_text(json.dumps({'design':'full'}))
            with self.assertRaisesRegex(ValueError,'design differs'):
                check_existing_output(out,'reduced')
            check_existing_output(out,'full')

if __name__=='__main__': unittest.main()
