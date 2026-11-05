import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import core
import evaluation
import final_comparison
import guided_followup
from make_smoke import fixture


class MetricGrouping(unittest.TestCase):
    def test_only_EAS_SEA_confusion_is_ignored_and_inputs_are_unchanged(self):
        sources = ["EUR", "EAS", "SSA", "SEA", "GME"]
        predicted = np.array([0, .3, 0, .7, 0])
        truth = np.array([0, 1, 0, 0, 0])
        before = predicted.copy(), truth.copy()
        self.assertAlmostEqual(core.ancestry_tv_pp(predicted, truth, sources), 70)
        self.assertEqual(core.ancestry_tv_pp(predicted, truth, sources, combine_eas_sea=True), 0)
        self.assertEqual(core.ancestry_tv_pp([1, 0, 0, 0, 0], [0, 0, 0, 0, 1], sources,
                                          combine_eas_sea=True), 100)
        np.testing.assert_array_equal(predicted, before[0])
        np.testing.assert_array_equal(truth, before[1])

    def test_grouping_option_requires_boolean_and_separate_axes(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = fixture(Path(temporary) / "inputs", {}, models=[{"id": "freq", "family": "nnls", "params": {}}])
            for value in (True, "true", 1):
                bad = copy.deepcopy(config)
                bad["combine_eas_sea_for_metrics"] = value
                with self.assertRaisesRegex(ValueError, "combine_eas_sea"):
                    core.validate_config(bad)
            config["combine_eas_sea_for_metrics"] = False
            core.validate_config(config)

    def test_grouped_metric_changes_ranking_not_eight_source_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = ["EUR", "EAS", "SEA"]
            samples = [{"sample_id": "val", "role": "validation", "family_id": "val_family", "track": "pure"}]
            for combined, expected in ((False, "B"), (True, "A")):
                config = {"version": 2, "mode": "development_screen", "seed": 42, "sources": sources,
                          "models": [{"id": name} for name in ("A", "B")],
                          "interpretations": [{"id": "native", "kind": "native"}],
                          "routes": [{"id": "mask", "stage": "harmonized"}],
                          "combine_eas_sea_for_metrics": combined, "independent_test": False}
                folders = []
                for name, values in (("A", [[0, .8, .2]]), ("B", [[.1, .9, 0]])):
                    folder = root / f"{combined}_{name}"
                    folder.mkdir()
                    np.savez_compressed(folder / "interpreted.npz", values=values, truth=[[0, 1, 0]])
                    core.write(folder / "interpretation.json", {
                        "config": config, "id": name, "route": config["routes"][0],
                        "interpreter": config["interpretations"][0], "samples": samples,
                        "sources": sources, "components": sources, "unscored_reason": None,
                        "license_track": "permissive_software", "sha256": core.digest(folder / "interpreted.npz")})
                    folders.append(folder)
                output = root / f"summary_{combined}"
                evaluation.select(folders, output)
                self.assertTrue(evaluation.verify(output)["verified"])
                self.assertEqual(core.read(output / "summary.json")["validation_selected_finalists"][0]["model"], expected)
                evidence = core.read(output / "prediction_evidence.json")
                self.assertEqual(evidence[0]["predictions"], [0, .8, .2])
                self.assertAlmostEqual(evidence[0]["eight_source_tv_percentage_points"], 20)
                self.assertAlmostEqual(evidence[0]["tv_percentage_points"], 0 if combined else 20)


class WinnerLock(unittest.TestCase):
    def test_test_is_planned_only_for_the_locked_validation_winner(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "genome/results"
            (results / "frozen/frozen").mkdir(parents=True)
            (results / "summary/summary").mkdir(parents=True)
            (root / "planning/design").mkdir(parents=True)
            (root / "testplanning").mkdir()
            original = root / "original.tsv"
            tests = [{"sample_id": f"test{i}", "role": "test", "track": "pure", "source": "EAS_SEA"}
                     for i in range(310)]
            guided_followup.write_rows(original, tests, list(tests[0]))
            development = root / "development.tsv"
            guided_followup.write_rows(development, [{"sample_id": "val", "role": "validation"}],
                                       ["sample_id", "role"])
            (root / "planning/design/segments.tsv").write_text(
                "sample_id\tchromosome\thaplotype\tstart\tend\tsource\n")
            core.write(results / "frozen/frozen/config.json", {
                "samples": str(development), "models": [{"id": name} for name in ("A", "B", "C")],
                "chromosomes": [22], "sources": ["EUR", "EAS", "SEA"]})
            leaderboard = [{"model": name, "validation_mean_person_tv_pp": value}
                           for name, value in (("A", 2), ("B", 3), ("C", 4))]
            summary = {"combine_eas_sea_for_metrics": True, "test_evaluated": False,
                       "validation_leaderboard": leaderboard}
            core.write(results / "summary/summary/summary.json", summary)
            core.write(results / "summary/summary/prediction_evidence.json", [{"role": "validation"}])
            spec = root / "spec.json"
            core.write(spec, {})
            authority = {row["sample_id"]: "EAS" for row in tests}
            def write_design(config, samples, segments, output, *args):
                output.mkdir()
                self.assertEqual(config["models"], [{"id": "A"}])
                self.assertEqual(sum(row["role"] == "test" for row in samples), 310)
                self.assertTrue(all(row["source"] == "EAS" for row in samples if row["role"] == "test"))
            with patch.object(final_comparison, "verify_manifest", return_value={"verified": True, "summary_sha256": "sha"}), \
                    patch.object(final_comparison, "source_inputs", return_value=(None, None, {"samples": str(original)},
                                                                                 None, {}, authority)), \
                    patch.object(final_comparison, "restore_segments", return_value=([], {})), \
                    patch.object(final_comparison, "write_design", side_effect=write_design):
                final_comparison.lock_and_test_plan(spec, root, root / "testplanning/design")
                self.assertEqual(core.read(root / "testplanning/design/winner-lock.json")["model"], {"id": "A"})
                summary["validation_leaderboard"][1]["validation_mean_person_tv_pp"] = 1
                core.write(results / "summary/summary/summary.json", summary)
                with self.assertRaisesRegex(ValueError, "Existing winner lock differs"):
                    final_comparison.lock_and_test_plan(spec, root, root / "changed")


class FinalReporting(unittest.TestCase):
    def test_report_requires_single_locked_model_and_cleanup_preserves_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = ["EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"]
            model = {"id": "A"}
            config = root / "test/results/frozen/frozen/config.json"
            config.parent.mkdir(parents=True)
            samples = [{"sample_id": f"test{i}", "source": sources[min(i, 7)]} for i in range(310)]
            guided_followup.write_rows(root / "samples.tsv", samples, list(samples[0]))
            core.write(config, {"models": [model], "sources": sources, "samples": str(root / "samples.tsv")})
            (root / "testplanning/design").mkdir(parents=True)
            core.write(root / "testplanning/design/winner-lock.json", {
                "model": model, "selection": {}, "development_summary_sha256": "sha"})
            summary = root / "test/results/summary/summary"
            summary.mkdir(parents=True)
            core.write(summary / "summary.json", {"combine_eas_sea_for_metrics": True,
                "validation_selected_finalists": [{"model": "A", "test_people": 310}]})
            rows = []
            for sample in samples:
                vector = np.eye(8)[sources.index(sample["source"])].tolist()
                for mask in ("one", "two"):
                    rows.append({"sample_id": sample["sample_id"], "role": "test", "track": "pure",
                                 "route": mask, "predictions": vector, "truth": vector,
                                 "eight_source_tv_percentage_points": 0})
            core.write(summary / "prediction_evidence.json", rows)
            with patch.object(final_comparison, "verify_manifest", return_value={"summary_sha256": "sha"}):
                (root / "reporting").mkdir()
                final_comparison.report(root, root / "reporting/report")
                report = core.read(root / "reporting/report/final.json")
                self.assertEqual(report["test_models_evaluated"], 1)
                self.assertEqual(report["test_eight_source_mean_TV_pp"], 0)
                self.assertEqual(report["EAS_SEA_discrimination"]["EAS"]["correct_per_mask_fraction"], 1)
                (root / "final-comparison").write_text("true\n")
                (root / "source.sha256").write_text("source evidence\n")
                for stage in ("planning", "genome", "testplanning", "test", "reporting"):
                    (root / stage / "work").mkdir(parents=True)
                final_comparison.cleanup(root, root / "cleanup.json")
                self.assertTrue(core.read(root / "cleanup.json")["durable_outputs_reverified"])
                self.assertTrue(config.is_file())
                self.assertTrue((root / "reporting/report/final.json").is_file())
                self.assertFalse((root / "genome/work").exists())
                core.write(config, {"models": [model, {"id": "B"}]})
                with self.assertRaisesRegex(ValueError, "locked single model"):
                    final_comparison.report(root, root / "wrong")


if __name__ == "__main__":
    unittest.main()
