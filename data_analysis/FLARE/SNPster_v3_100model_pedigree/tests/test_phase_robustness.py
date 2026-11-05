from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import model_family_gnomix as adapter
import screen_phase_robustness as robustness
from screen_phase_robustness import ARMS, phase_view, stable_seed


class PhaseRobustnessTests(unittest.TestCase):
    def test_exactly_two_fixed_arms(self):
        self.assertEqual(ARMS, {"phase_blocks_050": .5, "phase_blocks_200": 2.})

    def test_phase_labels_follow_genotypes_and_preserve_diploid_truth(self):
        rng = np.random.default_rng(17)
        x = rng.integers(0, 2, (8, 3000), dtype=np.int8)
        labels = np.repeat(np.arange(8)[:, None] % 7, 3000, axis=1)
        cm = np.linspace(0, 100, 3000)
        altered, target = phase_view(x, labels, cm, .5, 71)
        np.testing.assert_array_equal(altered.reshape(4, 2, -1).sum(axis=1),
                                      x.reshape(4, 2, -1).sum(axis=1))
        np.testing.assert_array_equal(np.sort(target.reshape(4, 2, -1), axis=1),
                                      np.sort(labels.reshape(4, 2, -1), axis=1))
        for i in range(4):
            swapped = target[2 * i] != labels[2 * i]
            np.testing.assert_array_equal(altered[2 * i, swapped], x[2 * i + 1, swapped])
            np.testing.assert_array_equal(altered[2 * i, ~swapped], x[2 * i, ~swapped])
        repeat, repeat_truth = phase_view(x, labels, cm, .5, 71)
        np.testing.assert_array_equal(altered, repeat)
        np.testing.assert_array_equal(target, repeat_truth)

    def test_no_switches_at_identical_genetic_positions_and_rates_are_measurable(self):
        x = np.tile([[0], [1]], (1, 50001))
        y = x.copy()
        for mean in ARMS.values():
            cm = np.linspace(0, 2000, x.shape[1])
            _, truth = phase_view(x, y, cm, mean, 19)
            switches = np.count_nonzero(np.diff(truth[0]))
            expected = (len(cm) - 1) * (-np.expm1(-2 * np.diff(cm)[0] / mean) / 2)
            self.assertLess(abs(switches - expected), 6 * np.sqrt(expected))
        _, truth = phase_view(x, y, np.zeros(x.shape[1]), .5, 19)
        self.assertEqual(np.count_nonzero(np.diff(truth[0])), 0)

    def test_invalid_genetic_map_and_missing_calls_fail(self):
        x = np.array([[0, 1], [1, 0]])
        for cm in ([2, 1], [np.nan, 1]):
            with self.assertRaises(ValueError):
                phase_view(x, x, np.asarray(cm), .5, 7)
        with self.assertRaises(ValueError):
            phase_view(x, x, np.array([0, 1]), 0, 7)
        with self.assertRaises(ValueError):
            phase_view(x * 2, x, np.array([0, 1]), 1, 7)
        self.assertEqual(stable_seed("a", 22), stable_seed("a", 22))
        self.assertNotEqual(stable_seed("a", 22), stable_seed("b", 22))

    def test_swapped_marker_truth_is_reduced_using_existing_window_mode(self):
        x = np.tile([[0], [1]], (1, 100))
        _, truth = phase_view(x, x, np.linspace(0, 10, 100), .5, 7)
        layout = dict(starts=np.array([0, 40]), ends=np.array([40, 100]))
        target = adapter.window_labels(truth, layout)
        expected = [[np.bincount(row[:40], minlength=7).argmax(),
                     np.bincount(row[40:], minlength=7).argmax()] for row in truth]
        np.testing.assert_array_equal(target, expected)

    def test_summary_exact_marker_weighting_and_separate_array_track(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            screen, pilot, out = [root / name for name in ("screen", "pilot", "out")]
            baseline = dict(pure={"mean_tv": .03}, pedigree={"mean_tv": .04})
            model_dir = screen / f"candidates/{robustness.BASE}"
            model_dir.mkdir(parents=True)
            adapter.write_json(model_dir / "metrics.json", baseline)
            pilot.mkdir()
            adapter.write_json(pilot / "summary.json",
                               {"stage_mean_tv": {"array_imputed": .25}})
            for arm in ARMS:
                for chrom in range(1, 23):
                    folder = out / f"candidates/{arm}/chr{chrom}"
                    folder.mkdir(parents=True)
                    n = 100000 if chrom == 1 else 4500 if chrom < 22 else 7260
                    eur = .2 if chrom == 1 else .6
                    predictions = []
                    for track, samples in (
                            ("pure", ["pure"]), ("pedigree", ["pedigree"]), ("real", ["real"]),
                            ("array_imputed", [f"control{i}" for i in range(5)])):
                        for sample in samples:
                            for estimand in ("marker", "physical"):
                                predictions.append(dict(
                                    dataset=track, sample_id=sample, estimand=estimand,
                                    n_markers=n, physical_span=1.,
                                    **dict(zip(adapter.SOURCES, [eur, 1 - eur, 0, 0, 0, 0, 0]))))
                    pd.DataFrame(predictions).to_csv(folder / "predictions.tsv", sep="\t", index=False)

            def fake_score(directory, arm):
                candidate = directory / f"candidates/{arm}"
                actual = pd.read_csv(candidate / "predictions.tsv", sep="\t")
                self.assertEqual(set(actual.dataset), {"pure", "pedigree", "real"})
                expected = (.2 * 100000 + .6 * 97260) / 197260
                np.testing.assert_allclose(actual.pred_EUR, expected)
                physical = pd.read_csv(candidate / "predictions_physical.tsv", sep="\t")
                np.testing.assert_allclose(physical.pred_EUR, (.2 + 21 * .6) / 22)
                adapter.write_json(candidate / "metrics.json", dict(
                    **baseline, all_gates_pass=False, gates={"real_world_regression": False},
                    real=[], real_constraints={}))

            with (patch.object(robustness, "load_experiment", return_value={
                    "input_sha256": {}, "augmentation_limitation": "synthetic phase only"}),
                  patch.object(robustness, "verified_checkpoint", return_value=True),
                  patch.object(robustness, "score", side_effect=fake_score)):
                robustness.summarize(screen, pilot, out)
            complete = json.loads((out / "completion.json").read_text())
            expected_tv = .5 - (.2 * 100000 + .6 * 97260) / 197260
            for result in complete["results"].values():
                self.assertAlmostEqual(result["array_imputed_mean_tv"], expected_tv)
                self.assertEqual(result["n_array_controls"], 5)
                self.assertEqual(result["failed_unchanged_gates"], ["real_world_regression"])


if __name__ == "__main__":
    unittest.main()
