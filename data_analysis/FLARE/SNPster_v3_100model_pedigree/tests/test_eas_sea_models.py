from pathlib import Path
import gzip
import json
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from screen_eas_sea_models import export_features, gates, population_split, source_diagnostics
from screen_v1_decoders import SOURCES


class EASSEAModelTests(unittest.TestCase):
    def test_export_executes_validation_and_preserves_sample_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            with gzip.open(directory / "run.global.anc.gz", "wt") as out:
                out.write("SAMPLE\tanc_0\tanc_1\nchild\t0.997\t0.002\n")
            self.assertEqual(export_features(directory, ["anc_0", "anc_1"]), 1)
            frame = pd.read_csv(directory / "features.tsv", sep="\t")
            self.assertEqual(frame.sample_id.tolist(), ["child"])
            self.assertAlmostEqual(frame.anc_0.iloc[0], .997 / .999)
            self.assertAlmostEqual(frame[["anc_0", "anc_1"]].sum(axis=1).iloc[0], 1.)
            audit = json.loads((directory / "output_audit.json").read_text())
            self.assertEqual(audit["n_samples"], 1)
            self.assertEqual(len(audit["output_sha256"]), 2)
            with self.assertRaises(ValueError):
                export_features(directory, ["anc_0"])

    def test_population_purge_removes_entire_shared_founder_families(self):
        pure = pd.DataFrame({"sample_id": list(SOURCES) + ["held"],
                             "source_group": list(SOURCES) + ["SEA"],
                             "fine": ["keep_" + s for s in SOURCES] + ["KHV"]})
        training = pd.DataFrame({"sample_id": ["c1", "c2", "c3"],
                                  "family_id": ["bad", "bad", "good"]})
        donors = {"bad": {"held", "EUR"}, "good": {"EAS", "SEA"}}
        p, t, audit = population_split(pure, training, donors, "KHV")
        self.assertNotIn("held", set(p.sample_id))
        self.assertEqual(t.sample_id.tolist(), ["c3"])
        self.assertEqual(audit["excluded_families"], ["bad"])
        self.assertNotIn("held", audit["training_founder_donors"])
        self.assertNotIn("KHV", audit["training_populations"])
        with self.assertRaises(ValueError):
            population_split(pure, training, {"bad": {"unknown"}, "good": {"EUR"}}, "KHV")
        with self.assertRaises(ValueError):
            population_split(pure, training, donors, "keep_EUR")

    def test_source_diagnostics_balance_sources_not_sample_counts(self):
        frame = pd.DataFrame({"source_group": ["EAS", "EAS", "SEA"],
                              "TV_distance": [0., 0., 1.],
                              "dominant_correct": [True, True, False]})
        result = source_diagnostics(frame)
        self.assertEqual(result["balanced_dominant_accuracy"], .5)
        self.assertEqual(result["sea_mean_tv"], 1.)
        with self.assertRaises(ValueError):
            source_diagnostics(frame.iloc[:2])

    def test_unchanged_gates_enforce_tail_real_and_minor_support(self):
        cfg = {"min_reference_samples": 100, "min_pedigree_samples": 20,
               "required_scenarios": ["reference_source", "f1"],
               "max_gates": {"balanced_mean_tv": .2, "balanced_p95_tv": .4},
               "min_gates": {}, "minor_gate_min_support": 20,
               "minor_recall_min": .4, "require_real_hard_gate": True}
        metrics = {"n_reference": 159, "n_pedigree": 123,
                   "scenarios": [{"scenario": "reference_source"}, {"scenario": "f1"}],
                   "balanced_mean_tv": .2, "balanced_p95_tv": .40001,
                   "minor_2_5_n": 19, "minor_2_5_recall": 0.}
        result, evaluable = gates(metrics, False, cfg)
        self.assertFalse(evaluable)
        self.assertNotIn("minor_2_5_recall", result)
        self.assertFalse(result["balanced_p95_tv"])
        self.assertFalse(result["real_world_regression"])
        metrics["minor_2_5_n"] = 20
        result, evaluable = gates(metrics, True, cfg)
        self.assertTrue(evaluable)
        self.assertFalse(result["minor_2_5_recall"])
        metrics["n_reference"] = 99
        with self.assertRaises(ValueError):
            gates(metrics, True, cfg)


if __name__ == "__main__":
    unittest.main()
