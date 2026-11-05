from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from audit_representations import (
    SOURCES, binary_metrics, collapsed_eas_sea_tv, expected_mixtures, fit_probe,
    reproduce,
)
from audit_feature_collisions import collisions, disagreements


class IdentityDecoder:
    def predict(self, x):
        return x


class RepresentationAuditTests(unittest.TestCase):
    def test_exact_collisions_and_lower_bound_do_not_round_features(self):
        frame = pd.DataFrame({
            "sample_id": ["a", "b", "c", "d"], "source_group": ["EAS", "EAS", "SEA", "SEA"],
            "anc_0": [1., 1., 1., .999999], "anc_1": [0., 0., 0., .000001],
        })
        groups, metrics = collisions(frame, ["anc_0", "anc_1"])
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["n_samples"], 3)
        self.assertEqual(metrics["exact_collision_minimum_mean_tv"], .25)
        self.assertEqual(metrics["n_samples_in_ambiguous_groups"], 3)

    def test_genotype_contrast_ignores_phase_orientation_not_dosage(self):
        self.assertEqual(disagreements(["0|1", "0|0"], ["1|0", "1|1"]), 1)
        with self.assertRaisesRegex(ValueError, "aligned"):
            disagreements([], [])

    def test_collapsing_sources_only_removes_within_pair_error(self):
        pred = np.zeros((2, 8))
        truth = np.zeros((2, 8))
        pred[:, SOURCES.index("EAS")] = 1
        truth[0, SOURCES.index("SEA")] = 1
        truth[1, SOURCES.index("EUR")] = 1
        np.testing.assert_allclose(collapsed_eas_sea_tv(pred, truth), [0, 1])

    def test_reproduction_aligns_ids_and_rejects_mismatches(self):
        pred = np.eye(2)
        old = pd.DataFrame({"sample_id": ["b", "a"], "one": [0, 1], "two": [1, 0]})
        reproduce(pred, ["a", "b"], old, ["one", "two"], "test")
        with self.assertRaisesRegex(ValueError, "cannot be reproduced"):
            reproduce(pred[::-1], ["a", "b"], old, ["one", "two"], "test")
        with self.assertRaisesRegex(ValueError, "IDs differ"):
            reproduce(pred, ["a", "c"], old, ["one", "two"], "test")

    def test_linear_decoder_has_zero_mixture_gap(self):
        cols = [f"anc_{i}" for i in range(8)]
        features = pd.DataFrame(np.eye(8)[:2], columns=cols)
        features["sample_id"] = ["donor1", "donor2"]
        plan = pd.DataFrame([
            {"sample_id": "f1", "generation": 1, "is_founder": True, "founder_donor": "donor1"},
            {"sample_id": "f2", "generation": 1, "is_founder": True, "founder_donor": "donor2"},
            {"sample_id": "child", "generation": 2, "is_founder": False,
             "father_id": "f1", "mother_id": "f2"},
        ])
        x, p = expected_mixtures(plan, features, cols,
                                 {"decoder": IdentityDecoder(), "calibrator": None}, ["child"])
        np.testing.assert_allclose(x, p, atol=1e-10)
        np.testing.assert_allclose(x[0, :2], [.5, .5])
        with self.assertRaisesRegex(ValueError, "Missing founder"):
            expected_mixtures(plan, features.iloc[:1], cols,
                              {"decoder": IdentityDecoder(), "calibrator": None}, ["child"])

    def test_fixed_probe_uses_both_sources_and_valid_metrics(self):
        x = np.array([[1., 0.], [.9, .1], [0., 1.], [.1, .9]])
        labels = np.array([0, 0, 1, 1])
        model = fit_probe(x, labels)
        result = binary_metrics(labels, model.predict_proba(x)[:, 1])
        self.assertEqual(result["balanced_accuracy"], 1.)
        self.assertEqual(result["roc_auc"], 1.)
        with self.assertRaisesRegex(ValueError, "both sources"):
            fit_probe(x, np.zeros(4, dtype=int))
        with self.assertRaisesRegex(ValueError, "both sources"):
            binary_metrics([1, 1], np.array([.9, .9]))


if __name__ == "__main__":
    unittest.main()
