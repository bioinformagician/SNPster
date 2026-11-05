from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from screen_precision_retraining import paired_features, replace_features


class PrecisionRetrainingTests(unittest.TestCase):
    def test_paired_quantization_retains_precision_and_normalizes_rounded_rows(self):
        frame = pd.DataFrame({"anc_0": [.9997, .99746], "anc_1": [.0003, .00044],
                              "anc_2": [0., .00044], "anc_3": [0., .00166]})
        paired = paired_features(frame)
        np.testing.assert_allclose(paired["precision"], frame)
        np.testing.assert_allclose(paired["rounded"][0], [1., 0., 0., 0.])
        np.testing.assert_allclose(paired["rounded"][1], np.array([.997, 0., 0., .002]) / .999)
        for arm in paired.values():
            np.testing.assert_allclose(arm.sum(axis=1), 1.)
        for bad in (pd.DataFrame([[np.nan, 1.]]), pd.DataFrame([[-.1, 1.1]]),
                    pd.DataFrame([[0., 0.]])):
            with self.assertRaises(ValueError):
                paired_features(bad)

    def test_feature_replacement_aligns_ids_without_changing_truth_or_order(self):
        old = pd.DataFrame({"sample_id": ["a", "b"], "anc_0": [1., 0.],
                            "anc_1": [0., 1.], "truth": [.7, .3]})
        new = pd.DataFrame({"sample_id": ["b", "a"], "anc_0": [.1, .9], "anc_1": [.9, .1]})
        result = replace_features(old, new, ["anc_0", "anc_1"])
        self.assertEqual(result.sample_id.tolist(), ["a", "b"])
        np.testing.assert_array_equal(result.truth, old.truth)
        np.testing.assert_allclose(result.anc_0, [.9, .1])
        np.testing.assert_allclose(old.anc_0, [1., 0.])
        for bad in (new.iloc[:1], pd.concat([new, new.iloc[:1]])):
            with self.assertRaises(ValueError):
                replace_features(old, bad, ["anc_0", "anc_1"])


if __name__ == "__main__":
    unittest.main()
