import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from screen_v1_decoders import (
    family_donors, fractions, purged_folds, fit_stochastic, predict_decoder,
    training_weights, exact_join, reference_indices,
)


class DecoderScreenTests(unittest.TestCase):
    def test_shared_founders_are_purged_not_just_test_family(self):
        frame = pd.DataFrame({"family_id": ["a", "b", "c", "c"]})
        donors = {"a": {"one", "two"}, "b": {"two", "three"}, "c": {"four"}}
        folds = purged_folds(frame, donors)
        family, train, test, audit = folds[0]
        self.assertEqual(family, "a")
        self.assertEqual(train.tolist(), [2, 3])
        self.assertEqual(test.tolist(), [0])
        self.assertEqual(audit["donor_overlap"], 0)
        for family, train, test, audit in folds:
            self.assertFalse(set(audit["training_donors"]) & donors[family])

    def test_no_disjoint_training_families_fails(self):
        frame = pd.DataFrame({"family_id": ["a", "b"]})
        with self.assertRaisesRegex(ValueError, "no founder-disjoint"):
            purged_folds(frame, {"a": {"same"}, "b": {"same"}})

    def test_reference_donor_in_pedigree_training_is_excluded(self):
        reference = pd.DataFrame({"sample_id": ["heldout", "training_donor", "unrelated"]})
        self.assertEqual(reference_indices(reference, ["training_donor"], True).tolist(), [0, 2])
        self.assertEqual(reference_indices(reference, ["training_donor"], False).tolist(), [0, 1, 2])

    def test_missing_founders_and_within_family_reuse_fail(self):
        plan = pd.DataFrame({
            "family_id": ["a", "a"], "is_founder": [True, True],
            "founder_donor": ["same", "same"], "founder_source": ["EUR", "EUR"],
        })
        with self.assertRaisesRegex(ValueError, "reused"):
            family_donors(plan)
        with self.assertRaisesRegex(ValueError, "Missing"):
            family_donors(plan.assign(is_founder=False))

    def test_stochastic_decoder_is_simplex_and_mixture_consistent(self):
        x = np.eye(3)
        y = np.eye(3)
        model = fit_stochastic(x, y, np.ones(3), alpha=1e-8)
        pred = predict_decoder(model, x)
        np.testing.assert_allclose(pred, y, atol=1e-5)
        np.testing.assert_allclose(model["matrix"].sum(axis=1), 1, atol=1e-8)
        mixed = predict_decoder(model, np.array([[.25, .75, 0]]))
        np.testing.assert_allclose(mixed[0], .25 * pred[0] + .75 * pred[1], atol=1e-10)
        self.assertTrue((mixed >= 0).all())

    def test_weights_balance_sources_and_families(self):
        pure = pd.DataFrame({"source_group": ["EUR", "EUR", "SEA"],
                             "sample_id": ["a", "b", "c"]})
        ped = pd.DataFrame({"family_id": ["one", "two", "two"],
                            "sample_id": ["d", "e", "f"]})
        weights = training_weights(pure, ped)
        np.testing.assert_allclose(weights, [.125, .125, .25, .25, .125, .125])
        self.assertAlmostEqual(weights.sum(), 1)

    def test_invalid_fractions_and_mismatched_ids_fail(self):
        for bad in ([[np.nan, 0]], [[-.1, 1.1]], [[.1, .2]]):
            with self.assertRaises(ValueError):
                fractions(bad, "test")
        with self.assertRaisesRegex(ValueError, "IDs differ"):
            exact_join(pd.DataFrame({"sample_id": ["a"]}),
                       pd.DataFrame({"sample_id": ["b"]}))


if __name__ == "__main__":
    unittest.main()
