from pathlib import Path
import os
import sys
import unittest
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import audit_phasefree_readiness as readiness


class PhasefreeReadinessTests(unittest.TestCase):
    def test_unused_affinity_samples_are_not_one_hot_truth(self):
        pools = [
            ("reference_fit", [{"sample_id": "ref", "fine": "1KG_CEU", "source_group": "EUR"}],
             "sample_id"),
            ("previous_founders", [{"sample": "previous"}], "sample"),
        ]
        roles = [
            dict(sample="ref", fine="1KG_CEU", role="source"),
            dict(sample="previous", fine="1KG_CEU", role="source"),
            dict(sample="new", fine="1KG_CEU", role="source"),
            dict(sample="admixed", fine="1KG_MXL", role="affinity_only"),
        ]
        report, rows = readiness.unused_donors(pools, roles, {"ref", "new", "admixed"})
        self.assertEqual(report["unused_current_known_source_dense_identities"], 1)
        self.assertEqual(report["unused_annotated_dense_identities"], 2)
        self.assertEqual([row["sample_id"] for row in rows
                          if row["eligible_for_known_source_truth"]], ["new"])

    def test_inconsistent_source_labels_fail_but_reporting_merge_is_allowed(self):
        records = [
            dict(sample_id="a", fine="CDX", source_group="EAS"),
            dict(sample_id="b", fine="CDX", source_group="SEA"),
        ]
        readiness.unused_donors([("pool", records, "sample_id")], [], set())
        records[1]["source_group"] = "EUR"
        with self.assertRaisesRegex(ValueError, "Inconsistent current population"):
            readiness.unused_donors([("pool", records, "sample_id")], [], set())
        with self.assertRaisesRegex(ValueError, "identity overlap"):
            readiness.unused_donors([
                ("reference", records[:1], "sample_id"),
                ("validation", records[:1], "sample_id"),
            ], [], set())

    def test_query_requires_exact_alleles_and_sample_and_complete_axis(self):
        keys = [("22", 10, "A", "G")]
        valid = "22\t10\tA\tG\tperson\t0|1\t1\t0,1,0\n"
        with patch.object(readiness.subprocess, "check_output", return_value=valid):
            frame = readiness.qc_frame(Path("input.vcf.gz"), "person", keys)
        self.assertEqual(frame.GT.tolist(), ["0|1"])
        for text, message in (
            (valid.replace("\tG\t", "\tC\t"), "Missing"),
            (valid.replace("person", "other"), "incorrect sample"),
            (valid + valid, "Duplicate"),
            (valid.replace("0|1", "./."), "Complete biallelic"),
            (valid.replace("0,1,0", "0,nan,0"), "Invalid biallelic"),
        ):
            with (patch.object(readiness.subprocess, "check_output", return_value=text),
                  self.assertRaisesRegex(ValueError, message)):
                readiness.qc_frame(Path("input.vcf.gz"), "person", keys)

    def test_actual_production_qc_boundaries_and_zero_kept(self):
        if "SNPSTER_QC_MODULE" in os.environ:
            directory = Path(os.environ["SNPSTER_QC_MODULE"])
        else:
            directory = Path(__file__).resolve().parents[4] / "data pipeline/qc_module"
        with patch.dict(os.environ, {
            "GP_MIN": ".90", "DS_TOL": ".05", "SNPS_ONLY": "True", "BIALLELIC_ONLY": "True",
        }):
            method, thresholds = readiness.load_production_qc(directory)
        frame = pd.DataFrame([
            ["A", "G", "0|0", .05, .95, .05, 0],
            ["A", "G", "0|0", .051, .95, .05, 0],
            ["A", "G", "0|1", .95, .05, .90, .05],
            ["A", "G", "0|1", 1, .05, .89, .06],
            ["A", "G", "1|1", 1.95, 0, .05, .95],
            ["A", "G", "1|1", 1.949, 0, .05, .95],
            ["A", "GG", "0|0", 0, 1, 0, 0],
        ], columns=["REF", "ALT", "GT", "DS", "GP_00", "GP_01", "GP_11"])
        self.assertEqual(readiness.qc_kept_indices(frame, method, thresholds), {0, 2, 4})
        self.assertEqual(readiness.qc_kept_indices(frame.iloc[[1]], method, thresholds), set())
        self.assertEqual(len(frame), 7)


if __name__ == "__main__":
    unittest.main()
