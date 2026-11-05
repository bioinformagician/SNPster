import sys, unittest, tempfile, csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pedigree_benchmark'/'scripts'))
from prepare_founders import prepare, SOURCES

class FounderTests(unittest.TestCase):
    def test_intersect_heldout_and_exclude_train(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            val=d/'val.tsv'; train=d/'train.tsv'; out=d/'eligible.tsv'
            with val.open('w') as f:
                f.write('sample\tfine\tsource_group\n')
                for source in SOURCES:
                    for i in range(3): f.write(f'{source}_{i}\tF{i}\t{source}\n')
                f.write('AFFINITY_0\tF\tCAS\n')
            train.write_text('sample\nEUR_0\n')
            vcf={f'{s}_{i}' for s in SOURCES for i in range(3)} | {'AFFINITY_0'}
            res=prepare(val,train,vcf,out,2,reference_samples={'EUR_1'})
            self.assertEqual(res['flare_reference_overlap_count'],1)
            self.assertEqual(res['selected_donors'],23)
            with out.open() as fh:
                ids={r['sample'] for r in csv.DictReader(fh,delimiter='\t')}
            self.assertNotIn('EUR_0',ids)
            self.assertNotIn('AFFINITY_0',ids)
            with self.assertRaisesRegex(ValueError,'Insufficient'):
                prepare(val,train,vcf,out,3)

if __name__=='__main__': unittest.main()
