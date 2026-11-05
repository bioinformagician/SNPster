import hashlib
import json
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import phasefree_validation as validation
import sweep_phasefree_calibration as sweep
from verify_phasefree_sweep import independently_balanced_tv


class PhasefreeSweepTests(unittest.TestCase):
    def fixture(self):
        pure = [dict(sample_id=f"{source}_{index}", source_group=source, track="pure",
                     family_id=f"pure:{source}_{index}", scenario="reference_source")
                for source in validation.SOURCES for index in range(6)]
        ordered = {}
        for source in ("EUR", "SAS"):
            ordered[source] = sorted(
                [row["sample_id"] for row in pure if row["source_group"] == source],
                key=lambda identity: hashlib.sha256(f"123:{identity}".encode()).digest())
        donors, pedigree = {}, []
        for index in range(3):
            family = f"family{index}"
            donors[family] = {ordered["EUR"][index], ordered["SAS"][index]}
            pedigree += [dict(
                sample_id=f"{family}_child{child}", source_group="", track="pedigree",
                family_id=family, scenario="f1") for child in range(index + 1)]
        return pd.DataFrame(pure + pedigree), donors

    def test_donor_purging_and_exactly_once_out_of_fold_coverage(self):
        metadata, donors = self.fixture()
        folds = validation.donor_purged_splits(metadata, donors, 3, 123)
        self.assertEqual(len(folds), 6)
        tested = [person for row in folds for person in row["test_ids"]]
        self.assertEqual(len(tested), len(metadata))
        self.assertEqual(set(tested), set(metadata.sample_id))
        for row in folds:
            training = metadata[metadata.sample_id.isin(row["train_ids"])]
            used = set(training.loc[training.track.eq("pure"), "sample_id"])
            for family in training.loc[training.track.eq("pedigree"), "family_id"]:
                used.update(donors[family])
            self.assertFalse(used & set(row["heldout_donors"]))

    def test_weights_are_source_and_family_balanced_not_row_balanced(self):
        metadata, _ = self.fixture()
        metadata = metadata.sample(frac=1, random_state=7).reset_index(drop=True)
        weights = validation.fold_weights(metadata)
        self.assertAlmostEqual(weights[metadata.track.eq("pure")].sum(), .5)
        self.assertAlmostEqual(weights[metadata.track.eq("pedigree")].sum(), .5)
        for source in validation.SOURCES:
            self.assertAlmostEqual(weights[metadata.source_group.eq(source)].sum(), .5 / 7)
        for family in metadata.loc[metadata.track.eq("pedigree"), "family_id"].unique():
            self.assertAlmostEqual(weights[metadata.family_id.eq(family)].sum(), .5 / 3)

    def test_invalid_split_metadata_is_rejected(self):
        metadata, donors = self.fixture()
        with self.assertRaises(ValueError):
            validation.donor_purged_splits(pd.concat([metadata, metadata.iloc[:1]]), donors, 3, 123)
        donors["family0"].add("not_a_training_donor")
        with self.assertRaises(ValueError):
            validation.donor_purged_splits(metadata, donors, 3, 123)

    def test_feature_stabilization_is_sub_ulp_invariant_and_never_fills_missing_values(self):
        x = np.tile(np.r_[np.ones(7) / 7, np.ones(7) / 7], (4, 1))
        expected = validation.stable_features(x)
        np.testing.assert_array_equal(expected, validation.stable_features(np.nextafter(x, np.inf)))
        np.testing.assert_array_equal(expected, validation.stable_features(np.nextafter(x, -np.inf)))
        x[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite paired"):
            validation.stable_features(x)

    def test_all_fixed_model_kinds_fit_round_trip_and_preserve_probability_simplex(self):
        protocol = json.loads(sweep.PROTOCOL.read_text())
        protocol["boosting"]["iterations"] = 3
        protocol = json.loads(json.dumps(protocol, sort_keys=True))
        random = np.random.default_rng(10)
        x = np.column_stack([random.dirichlet(np.ones(7), 42), random.dirichlet(np.ones(7), 42)])
        y = random.dirichlet(np.ones(7), 42)
        models = sweep.fit_models(x, y, np.repeat(1 / len(x), len(x)), protocol)
        self.assertEqual(set(models), set(protocol["models"]))
        for name, model in models.items():
            before = sweep.predict(model, x[:5])
            after = sweep.predict(pickle.loads(pickle.dumps(model)), x[:5])
            np.testing.assert_array_equal(before, after)
            self.assertTrue(np.all(before >= 0), name)
            np.testing.assert_allclose(before.sum(axis=1), 1, atol=1e-12, rtol=0)
            if model["parameters"]["round_features"]:
                np.testing.assert_array_equal(before, sweep.predict(model, np.nextafter(x[:5], np.inf)))

    def test_balanced_score_does_not_give_large_families_extra_weight(self):
        metadata = pd.DataFrame(dict(
            sample_id=["a", "b", "c"], family_id=["family1", "family2", "family2"],
            scenario=["f1"] * 3))
        truth = np.tile([.5, .5, 0, 0, 0, 0, 0], (3, 1))
        prediction = np.array([[.4, .6, 0, 0, 0, 0, 0],
                               [.3, .7, 0, 0, 0, 0, 0], [.3, .7, 0, 0, 0, 0, 0]])
        self.assertAlmostEqual(independently_balanced_tv(metadata, truth, prediction, .01), .15)


if __name__ == "__main__":
    unittest.main()
