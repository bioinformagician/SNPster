import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from diagnose_flare_native import anchored_model, command, select_controls
from model_family_gnomix import SOURCES


class FlareNativeDiagnosticTests(unittest.TestCase):
    def test_named_copying_axes_not_positional_identity(self):
        names = list(SOURCES)
        panels = list(reversed(SOURCES))
        template = dict(ancestries=names, panels=panels, p=np.full((7, 7), 1 / 7),
                        theta=np.full((7, 7), 1e-5), rho=np.arange(1, 8), generation=10)
        result = anchored_model(template, 1)
        for i, name in enumerate(names):
            self.assertEqual(result["p"][i, panels.index(name)], 1)
        np.testing.assert_array_equal(result["rho"], template["rho"])
        np.testing.assert_array_equal(result["theta"], template["theta"])
        np.testing.assert_array_equal(template["p"], np.full((7, 7), 1 / 7))
        result = anchored_model(template, .9, .001)
        np.testing.assert_allclose(result["p"].sum(axis=1), 1)
        np.testing.assert_array_equal(result["theta"], np.full((7, 7), .001))
        for value in (-1, 1.01):
            with self.assertRaises(ValueError):
                anchored_model(template, value)
        with self.assertRaises(ValueError):
            anchored_model(template, .9, 0)

    def test_hash_selection_is_source_complete_and_donor_disjoint(self):
        rows = [dict(sample_id=f"{source}{i}", source_group=source) for source in SOURCES for i in range(3)]
        a = select_controls(rows, {"evaluation"}, 5)
        b = select_controls(list(reversed(rows)), {"evaluation"}, 5)
        self.assertEqual(a, b)
        self.assertEqual({r["source_group"] for r in a}, set(SOURCES))
        with self.assertRaises(ValueError):
            select_controls(rows, {a[0]["sample_id"]}, 5)
        with self.assertRaises(ValueError):
            select_controls(rows + rows[:1], set(), 5)

    def test_native_model_is_explicit_and_no_truth_prior_leaks(self):
        cohort = dict(flare_jar="precision.jar", genetic_map="map",
                      reference=dict(vcf="reference", broad_panel="panel"))
        result = command(cohort, "target", "out", model="anchored.model", stock=True)
        self.assertEqual(result[result.index("-jar") + 1], "/opt/flare/flare.jar")
        self.assertIn("model=anchored.model", result)
        self.assertIn("em=false", result)
        self.assertIn("update-p=false", result)
        self.assertFalse(any(part.startswith("gt-ancestries=") for part in result))


if __name__ == "__main__":
    unittest.main()
