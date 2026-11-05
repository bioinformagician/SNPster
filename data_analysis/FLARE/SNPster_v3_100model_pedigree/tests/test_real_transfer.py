from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from audit_real_transfer import (
    crf_in_blocks, frequency_fits, genotype_posterior_qc, genotype_qc, randomize_phase, thin_sites,
)
import model_family_gnomix as adapter


class RealTransferTests(unittest.TestCase):
    def test_random_phase_preserves_diploid_genotypes_and_is_deterministic(self):
        x = np.random.default_rng(7).integers(0, 2, (10, 1000), dtype=np.int8)
        altered = randomize_phase(x, 19)
        np.testing.assert_array_equal(x.reshape(5, 2, -1).sum(axis=1),
                                      altered.reshape(5, 2, -1).sum(axis=1))
        np.testing.assert_array_equal(altered, randomize_phase(x, 19))
        self.assertFalse(np.array_equal(x, altered))
        with self.assertRaises(ValueError):
            randomize_phase(np.array([[0, 2], [1, 0]]), 19)

    def test_qc_distinguishes_missing_from_heterozygotes(self):
        result = genotype_qc(np.array([[0, 0, 1, 2], [0, 1, 0, 1]]))
        self.assertEqual(result["sites"][0], 4)
        self.assertEqual(result["missing_sites"][0], 1)
        self.assertEqual(result["heterozygous_sites"][0], 2)
        self.assertEqual(result["heterozygous_01"][0], 1)
        self.assertEqual(result["heterozygous_10"][0], 1)

    def test_frequency_fit_recovers_known_simplex_and_checks_support(self):
        rng = np.random.default_rng(5)
        frequencies = rng.uniform(.05, .95, (5000, 7))
        expected = np.array([[.5, .5, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1.]])
        fitted, support = frequency_fits(frequencies, 2 * expected @ frequencies.T)
        np.testing.assert_allclose(fitted, expected, atol=2e-5, rtol=0)
        self.assertEqual(support, 5000)
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            frequency_fits(frequencies[:50], np.zeros((1, 50)))
        with self.assertRaisesRegex(ValueError, "Invalid"):
            frequency_fits(frequencies, np.full((1, len(frequencies)), 3))

    def test_genotype_posterior_export_rounding_is_not_silently_repaired(self):
        gp = np.array([[.45, .45, .09], [.45, .45, .11], [1, 0, 0]])
        np.testing.assert_allclose(genotype_posterior_qc(gp), [.99, 1.01, 1])
        np.testing.assert_array_equal(gp[0], [.45, .45, .09])
        for values in ([[.4, .4, .4]], [[np.nan, 0, 1]], [[1.1, 0, -.1]]):
            with self.assertRaises(ValueError):
                genotype_posterior_qc(np.array(values))

    def test_observed_thinning_selects_first_eligible_site_in_each_bin(self):
        keys = [("1", pos, "A", "G") for pos in (10, 20, 100010, 100020)]
        result = thin_sites(keys, [False, True, False, True])
        np.testing.assert_array_equal(result, [False, True, False, True])
        with self.assertRaises(ValueError):
            thin_sites(keys, [True])

    def test_crf_blocks_keep_marker_order_and_restore_class_axis(self):
        class CRF:
            classes_ = [str(value) for value in (6, 0, 1, 2, 3, 4, 5)]

            def predict_proba(self, values):
                return values[:, :, [6, 0, 1, 2, 3, 4, 5]]
        class Smooth:
            model = CRF()
        class Model:
            mode = "fast"
            smooth = Smooth()
        x = np.random.default_rng(29).uniform(.1, 1, (4, 13, 7))
        x /= x.sum(axis=2, keepdims=True)
        for size in (1, 5, 13, 20):
            np.testing.assert_allclose(crf_in_blocks(Model(), x, size), x, atol=1e-15)

    def test_real_upstream_crf_blocks_match_full_inference(self):
        adapter.load_upstream()
        from src.Smooth.crf import CRF
        rng = np.random.default_rng(37)
        base = rng.uniform(.1, 1, (14, 13, 7))
        base /= base.sum(axis=2, keepdims=True)
        crf = CRF()
        crf.fit(base, np.repeat(np.repeat(np.arange(7), 2)[:, None], 13, axis=1))
        model = SimpleNamespace(smooth=SimpleNamespace(model=crf))
        order = np.argsort([int(value) for value in crf.classes_])
        expected = crf.predict_proba(base)[:, :, order]
        np.testing.assert_allclose(crf_in_blocks(model, base, 13), expected, atol=1e-12)
        blocked = crf_in_blocks(model, base, 5)
        self.assertEqual(blocked.shape, base.shape)
        np.testing.assert_allclose(blocked.sum(axis=2), 1, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
