from pathlib import Path
import pickle
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_supervised_genotypes as supervised
import test_matched_global as matched_tests


class SupervisedGenotypeTests(unittest.TestCase):
    def test_dual_solver_matches_sklearn_weighted_ridge_in_both_dimensions(self):
        rng = np.random.default_rng(17)
        for n, m in ((13, 31), (31, 13)):
            x, y = rng.normal(size=(n, m)), rng.normal(size=(n, 7))
            weights = rng.uniform(.1, 2, n)
            weights /= weights.sum()
            models = supervised.weighted_ridges(x, y, weights, [.1, 1, 10])
            for alpha, model in models.items():
                expected = Ridge(alpha=alpha, fit_intercept=True).fit(x, y, sample_weight=weights)
                np.testing.assert_allclose(model["coefficients"].T, expected.coef_, atol=1e-12)
                np.testing.assert_allclose(model["intercept"], expected.intercept_, atol=1e-12)
                self.assertLess(model["max_normal_equation_error"], 1e-10)

    def fixture(self, root):
        _, _, pure, ped = matched_tests.MatchedGlobalTests().training_fixture()
        ref = pd.DataFrame([
            dict(sample_id=f"reference{i}", source_group=supervised.SOURCES[i % 7],
                 fine=f"population{i}") for i in range(14)])
        ref.to_csv(root / "reference.tsv", sep="\t", index=False)
        axes = pd.DataFrame(dict(
            group=["reference"] * 14 + ["clean_training"] * 799
            + ["imputed_training"] * 7 + ["real"],
            sample_id=list(ref.sample_id) + list(pure.sample_id) + list(ped.sample_id)
            + list(ped.sample_id.iloc[:7]) + ["held_out"]))
        return axes, {"reference": {"labels_tsv": str(root / "reference.tsv")}}, pure, ped

    def test_policy_weights_and_training_labels_exclude_evaluation(self):
        with tempfile.TemporaryDirectory() as tmp:
            axes, cohort, pure, ped = self.fixture(Path(tmp))
            with patch.object(supervised, "training_table", return_value=(pure, ped)):
                rows, truth, weights, table = supervised.training(axes, cohort)
                self.assertNotIn(len(axes) - 1, rows)
                self.assertEqual(truth.shape, (820, 7))
                for values in weights.values():
                    self.assertAlmostEqual(values.sum(), 1)
                self.assertAlmostEqual(weights["reference"][:14].sum(), .5)
                self.assertAlmostEqual(weights["reference"][14:813].sum(), .5)
                self.assertAlmostEqual(weights["matched"][:14].sum(), .375)
                self.assertAlmostEqual(weights["matched"][-7:].sum(), .25)
                before = truth.copy()
                axes.loc[len(axes) - 1, "sample_id"] = "different_held_out"
                np.testing.assert_array_equal(supervised.training(axes, cohort)[1], before)
                axes.loc[813, "sample_id"] = "different_held_out"
                with self.assertRaisesRegex(ValueError, "nontraining"):
                    supervised.training(axes, cohort)

    def test_reference_identity_overlap_or_misordering_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            axes, cohort, pure, ped = self.fixture(root)
            with patch.object(supervised, "training_table", return_value=(pure, ped)):
                axes.loc[0, "sample_id"] = "wrong_reference"
                with self.assertRaisesRegex(ValueError, "Reference identities"):
                    supervised.training(axes, cohort)
                axes.loc[0, "sample_id"] = "reference0"
                axes.loc[len(axes) - 1, "sample_id"] = pure.sample_id.iloc[0]
                with self.assertRaisesRegex(ValueError, "Held-out"):
                    supervised.training(axes, cohort)

    def test_invalid_weights_labels_and_alpha_are_explicit_errors(self):
        x, y = np.eye(4), np.eye(4)
        for weights, alphas in ((np.ones(4), [1]), (np.array([0, .3, .3, .4]), [1]),
                                (np.repeat(.25, 4), [0]), (np.repeat(.25, 4), [])):
            with self.assertRaises(ValueError):
                supervised.weighted_ridges(x, y, weights, alphas)
        with self.assertRaises(ValueError):
            supervised.weighted_ridges(x * np.nan, y, np.repeat(.25, 4), [1])

    def test_saved_bundle_and_phase_independent_dosages_reproduce(self):
        rng = np.random.default_rng(31)
        x = rng.integers(0, 3, size=(40, 20), dtype=np.int8)
        mean, scale = x.mean(axis=0), np.ones(20)
        standardized = supervised.standardize(x, mean, scale).astype(float)
        y = np.eye(7)[np.arange(40) % 7]
        model = supervised.weighted_ridges(standardized, y, np.repeat(1 / 40, 40), [1])[1]
        bundle = dict(mean=mean, scale=scale, **model)
        restored = pickle.loads(pickle.dumps(bundle))
        np.testing.assert_array_equal(supervised.predict(bundle, x),
                                      supervised.predict(restored, x.copy()))
        self.assertEqual(supervised.predict(bundle, x).shape, (40, 7))
        np.testing.assert_allclose(supervised.predict(bundle, x).sum(axis=1), 1, atol=1e-12)
        with self.assertRaises(ValueError):
            supervised.predict(bundle, np.full((1, 20), .5))


if __name__ == "__main__":
    unittest.main()
