from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_kernel_genotypes as kernels


class KernelGenotypeTests(unittest.TestCase):
    def test_weighted_kernel_intercept_matches_independent_linear_ridge(self):
        rng = np.random.default_rng(43)
        x, y, test = rng.normal(size=(23, 17)), rng.normal(size=(23, 7)), rng.normal(size=(4, 17))
        weights = rng.uniform(.1, 2, len(x))
        weights /= weights.sum()
        fitted = kernels.kernel_heads(x @ x.T, y, weights, [1e-5, 1e-3, .1])
        for alpha, model in fitted.items():
            expected = Ridge(alpha=alpha).fit(x, y, sample_weight=weights)
            actual = (test @ x.T) @ model["dual"] + model["intercept"]
            np.testing.assert_allclose(actual, expected.predict(test), atol=1e-9, rtol=0)

    def test_all_three_kernels_are_symmetric_positive_semidefinite(self):
        x = np.random.default_rng(5).normal(size=(30, 19))
        dot, norm = x @ x.T, (x * x).sum(axis=1)
        for kind in kernels.KINDS:
            matrix = kernels.kernel_values(dot, norm, norm, 19, kind)
            np.testing.assert_allclose(matrix, matrix.T, atol=1e-14)
            self.assertGreaterEqual(np.linalg.eigvalsh(matrix).min(), -1e-10)
            if kind.startswith("rbf"):
                np.testing.assert_allclose(matrix.diagonal(), 1, atol=1e-14)
        np.testing.assert_allclose(kernels.kernel_values(dot, norm, norm, 19, "poly2"),
                                   (1 + dot / 19) ** 2)

    def test_invalid_geometry_or_training_inputs_fail_explicitly(self):
        for kind in ("unsupported", "rbf2"):
            with self.assertRaises(ValueError):
                kernels.kernel_values(np.ones((3, 3)) * 4, np.ones(3), np.ones(3), 2, kind)
        y = np.eye(3)
        for weights, alphas in ((np.ones(3), [1]), (np.array([0, .5, .5]), [1]),
                                (np.repeat(1 / 3, 3), [0]), (np.repeat(1 / 3, 3), [])):
            with self.assertRaises(ValueError):
                kernels.kernel_heads(np.eye(3), y, weights, alphas)

    def test_self_contained_bundle_and_diploid_orientation_are_invariant(self):
        rng = np.random.default_rng(7)
        x = rng.integers(0, 3, size=(40, 20), dtype=np.int8)
        mean, scale = x.mean(axis=0), np.ones(20)
        z = kernels.standardize(x, mean, scale).astype(float)
        norm = (z * z).sum(axis=1)
        matrix = kernels.kernel_values(z @ z.T, norm, norm, 20, "rbf2")
        head = kernels.kernel_heads(matrix, np.eye(7)[np.arange(40) % 7],
                                    np.repeat(1 / 40, 40), [1e-4])[1e-4]
        bundle = dict(kernel_kind="rbf2", anchor_dosages=x, mean=mean, scale=scale, **head)
        before = kernels.predict(bundle, x[:4])
        restored = kernels.predict(pickle.loads(pickle.dumps(bundle)), x[:4].copy())
        np.testing.assert_array_equal(before, restored)
        np.testing.assert_allclose(before.sum(axis=1), 1, atol=1e-12)
        from screen_dosage_models import genotype_dosage
        calls = ["0|0", "0|1", "1|0", "1|1"] * 5
        original = np.array([[genotype_dosage(call) for call in calls]], dtype=np.int8)
        changed = np.array([[genotype_dosage(call[::-1].replace("|", "/")) for call in calls]],
                           dtype=np.int8)
        np.testing.assert_array_equal(kernels.predict(bundle, original),
                                      kernels.predict(bundle, changed))


if __name__ == "__main__":
    unittest.main()
