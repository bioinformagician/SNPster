from pathlib import Path
import io
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import audit_phasefree_errors as audit


class PhasefreeErrorAuditTests(unittest.TestCase):
    def config(self):
        return dict(minor_min=.02, minor_max=.05, minor_detection_threshold=.01,
                    absent_ancestry_max=.01, minor_gate_min_support=20)

    def frame(self):
        row = dict(sample_id="person", family_id="family")
        for source in audit.SOURCES:
            row.update({f"{prefix}_{source}": 0. for prefix in (
                "truth", "pred", "classifier", "frequency")})
        row.update(truth_EUR=.93, truth_GME=.02, truth_SAS=.05,
                   pred_EUR=.99, pred_GME=.01, classifier_EUR=1.,
                   frequency_EUR=.93, frequency_GME=.02, frequency_SAS=.05)
        return pd.DataFrame([row])

    def test_support_counts_contributions_not_independent_people(self):
        summary = audit.cell_summary(audit.ancestry_cells(self.frame(), self.config()),
                                     self.config())
        self.assertEqual(summary["minor_contributions"], 2)
        self.assertEqual(summary["minor_people"], 1)
        self.assertEqual(summary["minor_families"], 1)
        self.assertFalse(summary["minor_gate_evaluable"])
        self.assertEqual(summary["views"]["decoder"]["minor_detected"], 1)
        self.assertEqual(summary["views"]["decoder"]["minor_recall"], .5)

    def test_no_support_is_explicit_and_false_positives_use_frozen_threshold(self):
        frame = self.frame()
        frame["truth_GME"], frame["truth_SAS"], frame["truth_EUR"] = 0., 0., 1.
        summary = audit.cell_summary(audit.ancestry_cells(frame, self.config()), self.config())
        self.assertEqual(summary["minor_contributions"], 0)
        self.assertIsNone(summary["views"]["decoder"]["minor_recall"])
        self.assertEqual(summary["views"]["decoder"]["absent_flagged"], 1)

    def test_duplicate_samples_and_invalid_fractions_fail(self):
        frame = self.frame()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            audit.ancestry_cells(pd.concat([frame, frame]), self.config())
        frame["pred_EUR"] = float("nan")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            audit.ancestry_cells(frame, self.config())
        frame["pred_EUR"] = 1.1
        with self.assertRaisesRegex(ValueError, "Invalid"):
            audit.ancestry_cells(frame, self.config())

    def test_same_person_across_clean_and_imputed_conditions_is_not_a_duplicate(self):
        frame = pd.DataFrame([
            dict(group="pedigree", sample_id="person"),
            dict(group="array_imputed", sample_id="person"),
        ])
        audit.validate_condition_axis(frame)
        with self.assertRaisesRegex(ValueError, "condition/sample"):
            audit.validate_condition_axis(pd.concat([frame, frame.iloc[:1]]))

    def test_feature_tsv_values_round_trip_across_tree_split_boundaries(self):
        text = "value\n0.12345678901234568\n"
        result = audit.read_table(io.StringIO(text))
        self.assertEqual(result.iloc[0]["value"], float("0.12345678901234568"))


if __name__ == "__main__":
    unittest.main()
