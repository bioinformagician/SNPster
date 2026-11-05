import json
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import sweep_phasefree_rich as rich
from sweep_phasefree_calibration import hgb
from verify_phasefree_rich import check_metrics


class PhasefreeRichTests(unittest.TestCase):
    def test_parent_output_provenance_uses_read_only_calibration_mount(self):
        self.assertEqual(rich.calibration_outputs({"/sweep/training.tsv": "digest"}),
                         {"/calibration/training.tsv": "digest"})
        with self.assertRaises(ValueError):
            rich.calibration_outputs({"/unexpected/training.tsv": "digest"})

    def fixture(self):
        protocol = json.loads(rich.PROTOCOL.read_text())
        protocol["boosting"]["iterations"] = 3
        protocol = json.loads(json.dumps(protocol, sort_keys=True))
        random = np.random.default_rng(12)
        x = np.column_stack([
            random.dirichlet(np.ones(7), 42), random.dirichlet(np.ones(7), 42),
            random.normal(size=(42, 64)), random.dirichlet(np.ones(100), 42)])
        return protocol, x, random.dirichlet(np.ones(7), 42), np.repeat(1 / 42, 42)

    def test_rich_features_keep_signed_pca_and_stabilize_one_ulp_changes(self):
        protocol, x, _, _ = self.fixture()
        expected = rich.canonical_features(x, protocol)
        self.assertTrue(np.any(expected[:, 14:78] < 0))
        for direction in (-np.inf, np.inf):
            np.testing.assert_array_equal(expected, rich.canonical_features(np.nextafter(x, direction), protocol))

    def test_missing_values_wrong_axes_and_broken_probability_blocks_fail(self):
        protocol, x, _, _ = self.fixture()
        bad = x.copy()
        bad[0, 20] = np.nan
        for values in (bad, x[:, :-1]):
            with self.assertRaises(ValueError):
                rich.canonical_features(values, protocol)
        bad = x.copy()
        bad[:, 78:] *= 2
        with self.assertRaisesRegex(ValueError, "fine-population simplex"):
            rich.canonical_features(bad, protocol)
        with self.assertRaises(ValueError):
            rich.genotype_features({}, np.array([[0, 1, np.nan]]), protocol)

    def test_ablation_columns_preserve_base_and_disjoint_new_blocks(self):
        protocol, _, _, _ = self.fixture()
        pca, fine = (set(rich.feature_columns(name, protocol)) for name in ("pca", "fine"))
        self.assertEqual(pca & fine, set(range(14)))
        self.assertEqual(pca | fine, set(range(178)))
        with self.assertRaises(ValueError):
            rich.feature_columns("unknown", protocol)

    def test_every_fixed_model_roundtrips_and_control_reproduces_old_head(self):
        protocol, x, y, weights = self.fixture()
        with threadpool_limits(limits=2):
            models = rich.fit_models(x, y, weights, protocol)
            original = hgb(protocol["head"], protocol).fit(
                rich.stable_features(x[:, :14]), y, sample_weight=weights)
            np.testing.assert_allclose(
                rich.predict(models[rich.CONTROL], x), rich.simplex_rows(original.predict(
                    rich.stable_features(x[:, :14]))), atol=1e-12, rtol=0)
            for model in models.values():
                before = rich.predict(model, x[:5])
                np.testing.assert_array_equal(before, rich.predict(pickle.loads(pickle.dumps(model)), x[:5]))
                for direction in (-np.inf, np.inf):
                    np.testing.assert_array_equal(before, rich.predict(model, np.nextafter(x[:5], direction)))
                self.assertTrue(np.all(before >= 0))
                np.testing.assert_allclose(before.sum(axis=1), 1, atol=1e-12, rtol=0)

    def test_linear_scaling_uses_only_supplied_training_rows_and_weights(self):
        protocol, x, y, weights = self.fixture()
        weights[:] = .5 / 41
        weights[0] = .5
        with threadpool_limits(limits=2):
            models = rich.fit_models(x, y, weights, protocol)
        values = rich.canonical_features(x, protocol)
        model = models["rich_combined_blend25"]
        np.testing.assert_allclose(model["mean"], weights @ values, atol=1e-12, rtol=0)
        self.assertTrue(np.all(model["scale"] >= .05))
        with self.assertRaises(ValueError):
            rich.fit_models(x, y, weights * 2, protocol)

    def test_independent_accuracy_checks_distinguish_tv_and_component_band(self):
        truth = np.array([[.5, .5, 0, 0, 0, 0, 0]])
        prediction = np.array([[.45, .45, .05, .05, 0, 0, 0]])
        report = dict(mean_tv=.1, fraction_all_components_within_10pp=1,
                      fraction_all_components_within_5pp=1)
        check_metrics(truth, prediction, report, 1e-12)
        report["mean_tv"] = .05
        with self.assertRaises(AssertionError):
            check_metrics(truth, prediction, report, 1e-12)


if __name__ == "__main__":
    unittest.main()
