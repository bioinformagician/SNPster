from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_dosage_models as dosage
import test_matched_global as matched_tests


class DosageModelTests(unittest.TestCase):
    def test_phased_unphased_and_swapped_calls_are_identical(self):
        for call, expected in (("0|0", 0), ("0|1", 1), ("1|0", 1), ("1|1", 2),
                               ("0/1", 1), ("1/0", 1)):
            self.assertEqual(dosage.genotype_dosage(call), expected)
            self.assertEqual(dosage.genotype_dosage(call[::-1]), expected)
        for invalid in (".|0", "./.", "1", "2|0", "0|1|0", "0:1"):
            with self.assertRaises(ValueError):
                dosage.genotype_dosage(invalid)

    def test_actual_bcftools_query_preserves_sample_allele_and_phase_axes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "fixture.vcf"
            source.write_text(
                '##fileformat=VCFv4.2\n##contig=<ID=1>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb\ta\n'
                '1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1|1\n'
                '1\t100010\t.\tC\tT\t.\tPASS\t.\tGT\t1|0\t0/0\n')
            compressed = root / "fixture.vcf.gz"
            subprocess.run(["bcftools", "view", "-Oz", "-o", compressed, source], check=True)
            subprocess.run(["bcftools", "index", "-t", compressed], check=True)
            regions = root / "positions.tsv"
            regions.write_text("1\t10\n1\t100010\n")
            keys = [("1", 10, "A", "G"), ("1", 100010, "C", "T")]
            x, altered = dosage.load_dosages(compressed, ["a", "b"], keys, regions, True)
            np.testing.assert_array_equal(x, [[2, 0], [1, 1]])
            np.testing.assert_array_equal(x, altered)
            with self.assertRaisesRegex(ValueError, "incompatible"):
                dosage.load_dosages(compressed, ["a", "b"],
                                     [("1", 10, "G", "A"), keys[1]], regions)

    def test_scaling_uses_supplied_reference_values_and_rejects_missing_calls(self):
        x = np.array([[0, 1], [2, 0]], dtype=np.int8)
        mean, scale = np.array([.4, .8]), np.array([.5, .6])
        np.testing.assert_allclose(dosage.standardize(x, mean, scale), (x - mean) / scale,
                                   atol=1e-7)
        for invalid in (np.array([[np.nan, 1]]), np.array([[3, 0]])):
            with self.assertRaises(ValueError):
                dosage.standardize(invalid, mean, scale)
        with self.assertRaises(ValueError):
            dosage.standardize(x, mean, np.array([0, 1]))

    def test_pca_model_bundle_roundtrip_and_unphased_genotype_equivalence(self):
        rng = np.random.default_rng(4)
        reference = rng.integers(0, 3, (120, 80), dtype=np.int8)
        mean = reference.mean(axis=0)
        scale = np.sqrt(2 * (mean / 2) * (1 - mean / 2))
        pca = PCA(n_components=8, whiten=True, random_state=3).fit(
            dosage.standardize(reference, mean, scale))
        scores = pca.transform(dosage.standardize(reference, mean, scale))
        y = np.eye(7)[np.arange(120) % 7]
        head = Ridge(alpha=.01).fit(scores, y)
        bundle = dict(kind="pca", mean=mean, scale=scale, pca=pca, head=head)
        before = dosage.predict(bundle, reference)
        after = dosage.predict(pickle.loads(pickle.dumps(bundle)), reference)
        np.testing.assert_array_equal(before, after)
        self.assertEqual(before.shape, (120, 7))
        np.testing.assert_allclose(before.sum(axis=1), 1, atol=1e-12)

    def test_all_supervised_heads_use_training_rows_only(self):
        _, _, pure, ped = matched_tests.MatchedGlobalTests().training_fixture()
        axes = pd.DataFrame(dict(
            group=["clean_training"] * 799 + ["imputed_training"] * 7 + ["real"],
            sample_id=list(pure.sample_id) + list(ped.sample_id)
            + list(ped.sample_id.iloc[:7]) + ["held_out"]))
        rng = np.random.default_rng(8)
        scores = rng.normal(size=(807, 8))
        from unittest.mock import patch
        with patch.object(dosage, "training_table", return_value=(pure, ped)):
            before = dosage.train_heads(scores, axes, {})
            changed = scores.copy()
            changed[-1] = 1e10
            after = dosage.train_heads(changed, axes, {})
        self.assertEqual(set(before), set(dosage.MODES[1:]))
        for name in before:
            np.testing.assert_array_equal(before[name].predict(scores[:20]),
                                          after[name].predict(scores[:20]))
        axes.loc[799, "sample_id"] = "held_out"
        with patch.object(dosage, "training_table", return_value=(pure, ped)):
            with self.assertRaisesRegex(ValueError, "nontraining"):
                dosage.train_heads(scores, axes, {})


if __name__ == "__main__":
    unittest.main()
