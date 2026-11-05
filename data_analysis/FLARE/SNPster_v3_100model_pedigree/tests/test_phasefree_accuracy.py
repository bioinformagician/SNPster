import json
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import report_phasefree_accuracy as report


class PhasefreeAccuracyTests(unittest.TestCase):
    def setUp(self):
        self.targets = json.loads(
            (report.PROJECT / "config/phasefree_accuracy_targets.json").read_text())

    def test_ten_and_five_percentage_point_boundaries_are_inclusive(self):
        truth = np.array([[.5, .5, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0]])
        prediction = np.array([[.4, .6, 0, 0, 0, 0, 0], [.95, .05, 0, 0, 0, 0, 0]])
        summary, rows = report.accuracy_summary(truth, prediction, self.targets)
        self.assertAlmostEqual(summary["mean_tv"], .075)
        self.assertTrue(summary["mean_tv_acceptable"])
        self.assertFalse(summary["mean_tv_desirable"])
        self.assertEqual(summary["n_all_components_within_10pp"], 2)
        self.assertEqual(summary["n_all_components_within_5pp"], 1)
        np.testing.assert_array_equal(rows["all_components_within_5pp"], [False, True])

    def test_tv_and_each_component_bands_are_not_conflated(self):
        truth = np.array([[.25, .25, .25, .25, 0, 0, 0]])
        prediction = np.array([[.31, .31, .19, .19, 0, 0, 0]])
        summary, _ = report.accuracy_summary(truth, prediction, self.targets)
        self.assertAlmostEqual(summary["mean_tv"], .12)
        self.assertFalse(summary["mean_tv_acceptable"])
        self.assertEqual(summary["n_all_components_within_10pp"], 1)
        self.assertEqual(summary["n_all_components_within_5pp"], 0)

    def test_invalid_missing_or_empty_fractions_fail(self):
        truth = np.array([[1., 0, 0, 0, 0, 0, 0]])
        for prediction in (
            np.array([[float("nan"), 0, 0, 0, 0, 0, 0]]),
            np.zeros((1, 7)), np.array([[1.1, -.1, 0, 0, 0, 0, 0]]),
            np.zeros((0, 7)), np.zeros((1, 6)),
        ):
            with self.assertRaises(ValueError):
                report.accuracy_summary(truth, prediction, self.targets)


if __name__ == "__main__":
    unittest.main()
