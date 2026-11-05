from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_reference_trees as trees


class ReferenceTreeTests(unittest.TestCase):
    def test_fixed_models_use_reproducible_training_not_validation_stopping(self):
        for parameters in trees.MODES.values():
            model = trees.model_for(parameters)
            if parameters["kind"] == "extratrees":
                self.assertEqual(model.n_estimators, 400)
                self.assertEqual(model.random_state, trees.adapter.SEED)
            else:
                self.assertFalse(model.estimator.early_stopping)
                self.assertEqual(model.estimator.random_state, trees.adapter.SEED)
                self.assertEqual(model.estimator.l2_regularization, 0)
        with self.assertRaises(ValueError):
            trees.model_for({"kind": "invalid"})

    def test_both_families_fit_weighted_global_truth_and_roundtrip(self):
        rng = np.random.default_rng(71)
        y = np.eye(7)[np.arange(140) % 7]
        x = y + rng.normal(0, .01, y.shape)
        weights = rng.uniform(.1, 2, len(x))
        weights /= weights.sum()
        with threadpool_limits(limits=2):
            for parameters in ({"kind": "extratrees", "leaf": 1},
                               {"kind": "boosting", "iterations": 150}):
                model = trees.model_for(parameters).fit(x, y, sample_weight=weights)
                predicted = model.predict(x[:5])
                np.testing.assert_allclose(predicted,
                    pickle.loads(pickle.dumps(model)).predict(x[:5]), atol=1e-12, rtol=0)
                self.assertEqual(predicted.shape, (5, 7))
                self.assertTrue(np.isfinite(predicted).all())


if __name__ == "__main__":
    unittest.main()
