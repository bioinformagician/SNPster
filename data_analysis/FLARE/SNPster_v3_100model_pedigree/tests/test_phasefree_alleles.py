from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from recover_phasefree_alleles import load_exact


class ExactAlleleRecoveryTests(unittest.TestCase):
    def vcf(self, directory, duplicate=False):
        path = Path(directory) / "one.vcf"
        rows = ["4\t164203107\t.\tG\tA\t.\tPASS\t.\tGT\t0|1",
                "4\t164203107\t.\tG\tT\t.\tPASS\t.\tGT\t0/0"]
        if duplicate:
            rows.append(rows[0])
        path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson\n"
                        + "\n".join(rows) + "\n")
        return path

    def test_same_position_distinct_alleles_keep_order_and_never_substitute(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.vcf(tmp)
            a, t = ("4", 164203107, "G", "A"), ("4", 164203107, "G", "T")
            np.testing.assert_array_equal(load_exact(path, "person", [t, a]), [0, 1])
            np.testing.assert_array_equal(load_exact(path, "person", [a]), [1])
            np.testing.assert_array_equal(load_exact(path, "person", [t]), [0])

    def test_missing_exact_alt_fails_even_when_position_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "Missing 1 exact-allele"):
                load_exact(self.vcf(tmp), "person", [("4", 164203107, "G", "C")])

    def test_duplicate_exact_records_and_wrong_samples_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            key = ("4", 164203107, "G", "A")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                load_exact(self.vcf(tmp, duplicate=True), "person", [key])
            with self.assertRaisesRegex(ValueError, "sample-axis"):
                load_exact(self.vcf(tmp), "wrong", [key])


if __name__ == "__main__":
    unittest.main()
