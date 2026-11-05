from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import verify_phasefree_arrays as verification


class PhasefreeArrayVerificationTests(unittest.TestCase):
    def test_dosage_is_phase_invariant_and_rejects_missing_or_haploid(self):
        self.assertEqual([verification.dosage(value) for value in ("0|1", "1|0", "1/0")],
                         [1, 1, 1])
        for value in ("./.", "1", "0|2", "0|1|1"):
            with self.assertRaisesRegex(ValueError, "verification genotype"):
                verification.dosage(value)

    def test_vcf_format_order_and_sample_identity_are_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.vcf"
            header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson\n"
            body = "22\t10\t.\tA\tG\t.\tPASS\t.\tDS:GT\t1:1/0\n"
            path.write_text(header + body)
            self.assertEqual(list(verification.one_sample_vcf(path, "person")),
                             [(("22", 10, "A", "G"), 1)])
            with self.assertRaisesRegex(ValueError, "sample-axis"):
                list(verification.one_sample_vcf(path, "other"))
            path.write_text(header + body + body)
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                list(verification.one_sample_vcf(path, "person"))
            path.write_text(body)
            with self.assertRaisesRegex(ValueError, "Invalid single-sample"):
                list(verification.one_sample_vcf(path, "person"))


if __name__ == "__main__":
    unittest.main()
