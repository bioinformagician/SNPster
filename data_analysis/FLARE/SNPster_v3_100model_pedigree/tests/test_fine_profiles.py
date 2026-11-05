from pathlib import Path
import pickle
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_fine_profiles as fine


class FineProfileTests(unittest.TestCase):
    def fixture(self):
        profiles = np.tile(np.eye(14), (12, 1))
        x = np.repeat(2 * profiles.T, 2, axis=0).astype(np.int8)
        labels = pd.DataFrame([
            dict(sample_id=f"r{index}-{replicate}", source_group=fine.SOURCES[index // 2],
                 fine=f"population{index}")
            for index in range(14) for replicate in range(2)])
        axes = labels[["sample_id"]].assign(group="reference")
        axes = pd.concat([axes, pd.DataFrame([dict(sample_id="held_out", group="real")])],
                         ignore_index=True)
        return np.vstack([x, x[0]]), axes, labels

    def test_fine_profiles_and_count_based_shrinkage_use_reference_only(self):
        x, axes, labels = self.fixture()
        raw = fine.build_profiles(x, axes, labels, 0)
        shrunk = fine.build_profiles(x, axes, labels, 10)
        self.assertEqual(len(raw["profiles"]), 14)
        self.assertTrue(all(row["n_reference"] == 2 for row in raw["profiles"]))
        for index, row in enumerate(raw["profiles"]):
            group_rows = labels.fine.eq(row["fine"])
            source_rows = labels.source_group.eq(row["source_group"])
            pop = x[:-1][group_rows].mean(axis=0) / 2
            broad = x[:-1][source_rows].mean(axis=0) / 2
            np.testing.assert_array_equal(raw["frequencies"][:, index], pop)
            np.testing.assert_allclose(shrunk["frequencies"][:, index], pop / 6 + 5 * broad / 6)
        x[-1] = 2
        after = fine.build_profiles(x, axes, labels, 10)
        np.testing.assert_array_equal(shrunk["frequencies"], after["frequencies"])
        np.testing.assert_array_equal(raw["informative"], shrunk["informative"])

    def test_pure_and_equal_mixture_recovery_and_population_aggregation(self):
        x, axes, labels = self.fixture()
        bundle = fine.build_profiles(x, axes, labels, 0)
        mixture = ((x[0].astype(int) + x[8]) // 2).astype(np.int8)
        probabilities, weights, diagnostics = fine.predict(bundle, np.vstack([x[0], mixture]))
        expected = np.zeros((2, 7))
        expected[0, 0] = 1
        expected[1, [0, 2]] = .5
        np.testing.assert_allclose(probabilities, expected, atol=2e-6)
        np.testing.assert_allclose(weights.sum(axis=1), 1, atol=1e-12)
        self.assertLess(diagnostics.simplex_dual_gap.max(), 1e-6)

    def test_initial_mass_is_source_balanced_despite_unequal_profile_counts(self):
        x, axes, labels = self.fixture()
        labels.loc[labels.fine.eq("population1"), "fine"] = "population0"
        bundle = fine.build_profiles(x, axes, labels, 0)
        totals = np.bincount(bundle["source_indices"], weights=bundle["initial_weights"],
                             minlength=7)
        np.testing.assert_allclose(totals, np.repeat(1 / 7, 7), atol=1e-15)

    def test_invalid_metadata_genotypes_and_shrinkage_fail_explicitly(self):
        x, axes, labels = self.fixture()
        with self.assertRaisesRegex(ValueError, "identities"):
            fine.build_profiles(x, axes, labels.iloc[::-1], 0)
        bad = labels.copy()
        bad.loc[4, "fine"] = bad.loc[0, "fine"]
        with self.assertRaisesRegex(ValueError, "assignments"):
            fine.build_profiles(x, axes, bad, 0)
        for strength in (-1, np.nan):
            with self.assertRaises(ValueError):
                fine.build_profiles(x, axes, labels, strength)
        bundle = fine.build_profiles(x, axes, labels, 0)
        for invalid in (np.full((1, 168), np.nan), np.full((1, 168), .5),
                        np.full((1, 168), 3), x[:1, :-1]):
            with self.assertRaises(ValueError):
                fine.predict(bundle, invalid)

    def test_serialization_and_phase_relabeling_preserve_predictions(self):
        x, axes, labels = self.fixture()
        bundle = fine.build_profiles(x, axes, labels, 10)
        before, _, _ = fine.predict(bundle, x[:3])
        after, _, _ = fine.predict(pickle.loads(pickle.dumps(bundle)), x[:3].copy())
        np.testing.assert_array_equal(before, after)
        from screen_dosage_models import genotype_dosage
        calls = ["0|0", "0|1", "1|0", "1|1"] * 42
        original = np.array([[genotype_dosage(call) for call in calls]], dtype=np.int8)
        changed = np.array([[genotype_dosage(call[::-1].replace("|", "/"))
                             for call in calls]], dtype=np.int8)
        np.testing.assert_array_equal(fine.predict(bundle, original)[0],
                                      fine.predict(bundle, changed)[0])

    def test_optimizer_failure_is_not_a_success_shaped_fallback(self):
        x, axes, labels = self.fixture()
        bundle = fine.build_profiles(x, axes, labels, 0)
        from types import SimpleNamespace
        with patch.object(fine, "minimize",
                          return_value=SimpleNamespace(success=False, message="test failure")):
            with self.assertRaisesRegex(ValueError, "test failure"):
                fine.predict(bundle, x[:1])


if __name__ == "__main__":
    unittest.main()
