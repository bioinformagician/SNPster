import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import joblib
from sklearn.dummy import DummyRegressor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pedigree_benchmark/scripts'))
from validate_frozen_candidates import SOURCES, sha256_file, validate, verify_reusable_inputs


class FrozenPilotTests(unittest.TestCase):
    def make_case(self, base):
        directory = base / 'candidate_c1'
        directory.mkdir()
        candidate = {'candidate_id': 'c1', 'representation': 'r1'}
        (directory / 'candidate.json').write_text(json.dumps(candidate))
        (directory / 'metrics.json').write_text(json.dumps({**candidate, 'feature_columns': ['anc_0']}))
        estimator = DummyRegressor(strategy='constant', constant=[1., 0., 0., 0., 0., 0., 0., 0.])
        estimator.fit([[1.]], [[1., 0., 0., 0., 0., 0., 0., 0.]])
        joblib.dump({'candidate': candidate, 'sources': SOURCES, 'feature_columns': ['anc_0'],
                     'decoder': estimator, 'calibrator': None}, directory / 'model.joblib')
        with (directory / 'predictions.tsv').open('w', newline='') as fh:
            fields = ['sample_id', 'role', *SOURCES, *[f'truth_{s}' for s in SOURCES]]
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerow({k: ('s1' if k == 'sample_id' else 'validation' if k == 'role' else
                                 1 if k in ('EUR', 'truth_EUR') else 0) for k in fields})
        model = base / 'r1.final.model'
        model.write_text('# frozen representation\nanc_0\n')
        config = base / 'r1.config.json'
        config.write_text(json.dumps({'rep_id': 'r1', 'min_maf': .005, 'min_mac': 50, 'flare_seed': 2}))
        manifest = base / 'input.tsv'
        manifest.write_text('candidate_id\trepresentation\tcandidate_dir\tmodel_path\tmodel_config\n'
                            f'c1\tr1\t{directory}\t{model}\t{config}\n')
        return manifest, model, directory

    def test_saved_pair_validation_and_identity_failure(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            manifest, model, directory = self.make_case(base)
            report = validate(manifest, base / 'validated.tsv', base / 'audit.json')
            self.assertEqual(report['candidate_count'], 1)
            self.assertFalse(report['models_retrained'])
            self.assertFalse(report['offspring_reused'])
            self.assertIn('model.joblib', report['frozen_artifact_sha256']['c1'])
            model.write_text('wrong_ancestry\n')
            with self.assertRaisesRegex(ValueError, 'ancestry columns'):
                validate(manifest, base / 'validated.tsv', base / 'audit.json')
            model.write_text('anc_0\n')
            (directory / 'candidate.json').write_text(json.dumps({'candidate_id': 'wrong', 'representation': 'r1'}))
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                validate(manifest, base / 'validated.tsv', base / 'audit.json')

    def test_reusable_inputs_reject_changed_sources_and_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            reuse = base / 'reuse'
            (reuse / 'maps').mkdir(parents=True)
            source_map = base / 'source.map'; source_map.write_text('map')
            source_ref = base / 'source.vcf'; source_ref.write_text('reference')
            (reuse / 'panel_sites.tsv').write_text('1\t100\t100\n')
            for chrom in range(1, 23):
                (reuse / f'maps/chr{chrom}.simmap').write_text(f'{chrom}\t100\t0\n{chrom}\t200\t1\n')
            files = [reuse / 'panel_sites.tsv', *sorted((reuse / 'maps').glob('*.simmap'))]
            manifest = {'genetic_map_sha256': sha256_file(source_map),
                        'reference_vcf_sha256': sha256_file(source_ref),
                        'files': {str(p.relative_to(reuse)): sha256_file(p) for p in files}}
            (reuse / 'manifest.json').write_text(json.dumps(manifest))
            self.assertEqual(verify_reusable_inputs(reuse, source_map, source_ref), manifest)
            source_ref.write_text('changed reference')
            with self.assertRaisesRegex(ValueError, 'source changed'):
                verify_reusable_inputs(reuse, source_map, source_ref)
            source_ref.write_text('reference')
            (reuse / 'maps/chr1.simmap').write_text('changed map')
            with self.assertRaisesRegex(ValueError, 'artifact changed'):
                verify_reusable_inputs(reuse, source_map, source_ref)

    def test_exact_ten_distinct_representation_manifest(self):
        with (ROOT / 'config/frozen_pilot_10.tsv').open() as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(len(rows), 10)
        self.assertEqual(len({r['representation'] for r in rows}), 10)
        self.assertIn('v1_frozen_k18', {r['candidate_id'] for r in rows})
        self.assertIn('gmm_diag_unrefined_K30_ef013d39e0__extratrees__ridge__cdeb886004',
                      {r['candidate_id'] for r in rows})


if __name__ == '__main__':
    unittest.main()
