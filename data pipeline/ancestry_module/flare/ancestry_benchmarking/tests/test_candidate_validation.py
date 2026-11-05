import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "fit_evaluate_candidate", ROOT / "scripts/fit_evaluate_candidate.py"
)
candidate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(candidate)


class CandidateValidationTests(unittest.TestCase):
    def test_weighted_fit_error_propagates_without_retry(self):
        class RejectWeights:
            calls = 0

            def fit(self, X, Y, sample_weight=None):
                self.calls += 1
                raise ValueError("invalid weights")

        model = RejectWeights()
        with self.assertRaisesRegex(ValueError, "invalid weights"):
            candidate.fit_with_optional_weight(model, [[1]], [[1]], [1])
        self.assertEqual(model.calls, 1)

    def test_unsupported_weights_are_rejected(self):
        class Unweighted:
            def fit(self, X, Y):
                raise AssertionError("must not fit without requested weights")

        with self.assertRaisesRegex(ValueError, "does not support required"):
            candidate.fit_with_optional_weight(Unweighted(), [[1]], [[1]], [1])

    def test_weighted_fit_preserves_weights(self):
        class Weighted:
            def fit(self, X, Y, sample_weight=None):
                self.weights = sample_weight

        model = Weighted()
        weights = np.array([2.0, 3.0])
        self.assertIs(candidate.fit_with_optional_weight(model, [[1], [2]], [[1], [2]], weights), model)
        np.testing.assert_array_equal(model.weights, weights)

    def test_unweighted_fit_remains_supported(self):
        class Unweighted:
            def fit(self, X, Y):
                self.fitted = True

        model = Unweighted()
        candidate.fit_with_optional_weight(model, [[1]], [[1]])
        self.assertTrue(model.fitted)

    def test_missing_enabled_features_fail_without_metrics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "candidate.json").write_text(json.dumps({"representation": "test_rep"}))
            pd.DataFrame([
                {"dataset_id": "disabled", "enabled": "false"},
                {"dataset_id": "validation", "enabled": "true"},
            ]).to_csv(root / "datasets.tsv", sep="\t", index=False)
            args = [
                "fit_evaluate_candidate.py",
                "--candidate", str(root / "candidate.json"),
                "--datasets", str(root / "datasets.tsv"),
                "--features-dir", str(root),
                "--real-benchmarks", str(root / "unused.tsv"),
                "--success-criteria", str(root / "unused.json"),
                "--outdir", str(root / "output"),
            ]
            with patch("sys.argv", args):
                with self.assertRaisesRegex(FileNotFoundError, "enabled dataset validation"):
                    candidate.main()
            self.assertFalse((root / "output/metrics.json").exists())

    def test_missing_selection_truth_fails_without_metrics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = {
                "representation": "test_rep", "candidate_id": "test",
                "decoder": {"type": "ridge"}, "calibrator": {"type": "none"},
            }
            (root / "candidate.json").write_text(json.dumps(config))
            pd.DataFrame({
                "sample_id": [f"s{i}" for i in range(8)],
                "source_group": candidate.SOURCES,
            }).to_csv(root / "truth.tsv", sep="\t", index=False)
            features = pd.DataFrame(np.eye(8), columns=[f"anc_{i}" for i in range(8)])
            features.insert(0, "sample_id", [f"s{i}" for i in range(8)])
            rows = []
            for did, role, truth in [
                ("train", "train_decoder", str(root / "truth.tsv")),
                ("validation", "validation", "."),
            ]:
                features.to_csv(root / f"test_rep__{did}.features.tsv", sep="\t", index=False)
                rows.append({
                    "dataset_id": did, "enabled": "true", "role": role,
                    "truth_path": truth, "truth_mode": "source_group",
                })
            pd.DataFrame(rows).to_csv(root / "datasets.tsv", sep="\t", index=False)
            (root / "benchmarks.tsv").write_text("benchmark_id\tdataset_id\thard_gate\n")
            (root / "criteria.json").write_text('{"selection_role": "validation"}')
            args = [
                "fit_evaluate_candidate.py",
                "--candidate", str(root / "candidate.json"),
                "--datasets", str(root / "datasets.tsv"),
                "--features-dir", str(root),
                "--real-benchmarks", str(root / "benchmarks.tsv"),
                "--success-criteria", str(root / "criteria.json"),
                "--outdir", str(root / "output"),
            ]
            with patch("sys.argv", args):
                with self.assertRaisesRegex(ValueError, "required selection role=validation"):
                    candidate.main()
            self.assertFalse((root / "output/metrics.json").exists())

    def test_fingerprint_changes_with_helper_contents(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "scripts").mkdir()
            (root / "config").mkdir()
            helper = root / "scripts/config_fingerprint.sh"
            shutil.copy2(ROOT / "scripts/config_fingerprint.sh", helper)
            (root / "config/search_space.json").write_text("{}")
            script = root / "scripts/fit_evaluate_candidate.py"
            script.write_text("version = 1\n")
            command = ["bash", str(helper)]
            before = subprocess.check_output(command, text=True)
            self.assertEqual(before, subprocess.check_output(command, text=True))
            script.write_text("version = 2\n")
            self.assertNotEqual(before, subprocess.check_output(command, text=True))


if __name__ == "__main__":
    unittest.main()
