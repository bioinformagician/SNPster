import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import benchmark_observed_ancestry as benchmark


class ObservedAncestryTests(unittest.TestCase):
    def test_phase_changes_preserve_typed_diploid_calls(self):
        before = np.array([[0, 1, 0], [1, 0, 0]])
        after = np.array([[1, 1, 0], [0, 0, 0]])
        keys = [("22", i, "A", "G") for i in (1, 2, 3)]
        self.assertEqual(benchmark.validate_phase(before, after, ["a"], ["a"], keys, keys), 3)
        for bad in (np.array([[1, 1, 1], [0, 0, 0]]), np.array([[2, 1, 0], [0, 0, 0]])):
            with self.assertRaises(ValueError):
                benchmark.validate_phase(before, bad, ["a"], ["a"], keys, keys)
        with self.assertRaises(ValueError):
            benchmark.validate_phase(before, after, ["a"], ["b"], keys, keys)
        with self.assertRaises(ValueError):
            benchmark.validate_phase(before, after, ["a"], ["a"], keys, keys[::-1])

    def test_phasing_and_native_flare_contracts(self):
        command = benchmark.phasing_command("gt", "ref", "map", "out")
        self.assertIn("impute=false", command)
        cohort = dict(flare_jar="flare.jar", reference=dict(vcf="ref", broad_panel="panel"),
                      genetic_map="map")
        for array in (True, False):
            for em in (True, False):
                command = benchmark.flare_command(cohort, "gt", "out", dict(array=array, em=em))
                self.assertIn(f"array={str(array).lower()}", command)
                self.assertIn(f"em={str(em).lower()}", command)
                self.assertIn("update-p=false", command)
                self.assertFalse(any("gt-ancestries=" in part or "model=" in part for part in command))

    def test_exact_alleles_are_not_interchangeable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mask.tsv"
            path.write_text("22\t10\tG\tA\n22\t10\tG\tT\n")
            keys = benchmark.typed_keys(path)
            self.assertEqual(len(keys), 2)
            path.write_text("22\t10\tG\tA\n22\t10\tG\tA\n")
            with self.assertRaises(ValueError):
                benchmark.typed_keys(path)

    def test_aggregation_is_physical_not_marker_weighted_and_rejects_partial_results(self):
        rows = []
        for stage in ("inherited_phase", "array_phase_only"):
            for chrom, span, source in ((1, 100, 0), (2, 300, 1)):
                p = np.eye(7)[source]
                rows.append(dict(candidate="flare", mask="array", stage=stage, sample_id="a",
                                 chromosome=chrom, span_bp=span, track="pure", family_id="a",
                                 **{f"pred_{s}": p[i] for i, s in enumerate(benchmark.SOURCES)}))
        frame = pd.DataFrame(rows)
        result = benchmark.aggregate(frame, [1, 2], ["flare"], ["a"], ["array"])
        np.testing.assert_allclose(result.pred_EUR, .25)
        np.testing.assert_allclose(result.pred_EAS_SEA, .75)
        for bad in (frame.iloc[:-1], pd.concat([frame, frame.iloc[:1]])):
            with self.assertRaises(ValueError):
                benchmark.aggregate(bad, [1, 2], ["flare"], ["a"], ["array"])
        with self.assertRaises(ValueError):
            benchmark.aggregate(frame, [1, 2], ["flare", "gnomix"], ["a"], ["array"])


if __name__ == "__main__":
    unittest.main()
