import json
from pathlib import Path
import sys
import unittest

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from audit_production_ancestry import cohort_audit


class ProductionPreflightTests(unittest.TestCase):
    def test_original_eas_sea_labels_are_not_lost(self):
        records = [
            dict(sample_id="dai", source_group="EAS_SEA", original_source_group="EAS", fine="CDX"),
            dict(sample_id="viet", source_group="EAS_SEA", original_source_group="SEA", fine="KHV"),
        ]
        result = cohort_audit(records, {"dai"}, {"dai": ("EAS", "CDX"), "viet": ("SEA", "KHV")})
        self.assertEqual(result["reference_overlap_count"], 1)
        self.assertEqual(result["direct_id_disjoint_counts"], {"SEA": 1})
        self.assertFalse(result["alias_and_relatedness_audit_complete"])

    def test_duplicate_or_conflicting_labels_fail(self):
        row = dict(sample_id="person", source_group="NAM", fine="population")
        authority = {"person": ("NAM", "population")}
        for records, labels in (([row, row], authority), ([row], {}), ([row], {"person": ("EUR", "population")})):
            with self.assertRaises(ValueError):
                cohort_audit(records, set(), labels)

    def test_merged_source_without_original_label_is_not_eight_source_truth(self):
        row = dict(sample_id="person", source_group="EAS_SEA", fine="population")
        with self.assertRaisesRegex(ValueError, "eight-source"):
            cohort_audit([row], set(), {"person": ("EAS", "population")})

    def test_launch_gate_cannot_silently_succeed(self):
        import subprocess
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "preflight.json"
            report.write_text(json.dumps({
                "status": "blocked_missing_independent_source_coverage",
                "missing_direct_id_disjoint_development_sources": ["NAM"],
            }))
            result = subprocess.run(
                [sys.executable, str(PROJECT / "pedigree_benchmark/scripts/audit_production_ancestry.py"),
                 "gate", "--output", str(report)], text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Scientific launch blocked", result.stderr)
            self.assertIn("NAM", result.stderr)


if __name__ == "__main__":
    unittest.main()
