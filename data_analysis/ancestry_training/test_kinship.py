import unittest
from pathlib import Path
import tempfile

import numpy as np

import core

from prepare_kinship import ancestry_totals, choose_roles, expected_vector, pedigree_closure


class KinshipDatasetTests(unittest.TestCase):
    def test_raw_same_position_alternates_require_exact_allele_intersection(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            header = ('##fileformat=VCFv4.2\n##contig=<ID=22>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdonor\n')
            (root / "raw.vcf").write_text(
                header + '22\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n'
                '22\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\n'
                '22\t20\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\n')
            (root / "sites.vcf").write_text(
                header + '22\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n'
                '22\t20\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\n')
            for name in ("raw", "sites"):
                core.run(["bcftools", "view", "-Oz", "-o", root / f"{name}.vcf.gz",
                          root / f"{name}.vcf"], root / name)
                core.run(["bcftools", "index", root / f"{name}.vcf.gz"], root / f"{name}.index")
            core.run(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz", "-o",
                      root / "exact.vcf.gz", root / "raw.vcf.gz", root / "sites.vcf.gz"], root / "exact")
            keys, calls = core.vcf(root / "exact.vcf.gz", 22, ["donor"], True)
            self.assertEqual(keys, [("22", 10, "A", "C"), ("22", 20, "G", "A")])
            np.testing.assert_array_equal(calls, [[[0, 1], [1, 0]]])

    def test_pinned_king_output_has_explicit_pair_and_marker_columns(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            calls = np.random.default_rng(7).integers(0, 2, (3, 2, 50), dtype=np.int8)
            core.write_vcf(root / "test.vcf.gz",
                           [("22", index + 1, "A", "C") for index in range(50)],
                           ["a", "b", "c"], calls)
            core.run(["plink2", "--threads", "1", "--memory", "640", "--vcf", root / "test.vcf.gz",
                      "--double-id", "--make-king-table", "counts", "--out", root / "king"], root / "run")
            rows = core.rows(root / "king.kin0")
            self.assertEqual(len(rows), 3)
            self.assertTrue({"IID1", "IID2", "NSNP", "KINSHIP"} <= rows[0].keys())
            self.assertTrue(all(int(row["NSNP"]) == 50 for row in rows))

    def test_parent_closure_preserves_intermediate_nodes(self):
        nodes = {"child": {"father_id": "father", "mother_id": "mother"},
                 "father": {"father_id": "", "mother_id": ""},
                 "mother": {"father_id": "", "mother_id": ""}}
        self.assertEqual(pedigree_closure(["child"], nodes), set(nodes))

    def test_unknown_and_cyclic_parents_fail(self):
        with self.assertRaises(ValueError):
            pedigree_closure(["child"], {"child": {"father_id": "absent", "mother_id": "absent"}})
        with self.assertRaises(ValueError):
            pedigree_closure(["child"], {"child": {"father_id": "child", "mother_id": "child"}})
        with self.assertRaises(ValueError):
            pedigree_closure(["child"], {"child": {"father_id": "", "mother_id": "absent"}})

    def test_expected_eas_sea_collapse_is_explicit(self):
        names = ["EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"]
        row = {f"expected_{source}": .5 if source in ("EAS", "SEA") else 0. for source in names}
        result = expected_vector(row, ["EUR", "EAS_SEA", "SAS", "GME", "NAM", "OCE", "SSA"])
        self.assertEqual(result["EAS_SEA"], 1.)
        self.assertEqual(sum(result.values()), 1.)

    def test_physical_truth_uses_two_homologs_not_parent_expectation(self):
        spans = {"1": [10, 20], "2": [0, 30]}
        segments = []
        for chrom, (start, end) in spans.items():
            for hap, source in [("0", "A"), ("1", "B")]:
                segments.append({"sample_id": "child", "chromosome": chrom, "haplotype": hap,
                                 "start": start, "end": end, "source": source})
        self.assertEqual(ancestry_totals(segments, ["A", "B"], spans)["child"], {"A": .5, "B": .5})
        segments[0]["start"] = 11
        with self.assertRaises(ValueError):
            ancestry_totals(segments, ["A", "B"], spans)

    def test_split_is_deterministic_and_never_splits_family(self):
        candidates = {}
        groups = []
        for index in range(4):
            members = [f"mixed{index}", f"donor{index}"]
            groups.append(members)
            for person in members:
                candidates[person] = {"track": "mixed" if person.startswith("mixed") else "pure",
                                      "source": "" if person.startswith("mixed") else "A", "role": "support"}
        for index in range(25):
            person = f"pure{index}"
            groups.append([person])
            candidates[person] = {"track": "pure", "source": "A", "role": "reference"}
        roles = choose_roles(groups, candidates, 1234, ["A"])
        self.assertEqual(roles, choose_roles(groups, candidates, 1234, ["A"]))
        self.assertEqual(set(roles.values()), {"reference", "train", "calibration", "validation", "test"})
        for group in groups:
            self.assertEqual(len({roles[person] for person in group}), 1)
        with self.assertRaises(ValueError):
            choose_roles(groups[:3], candidates, 1234, ["A"])


if __name__ == "__main__":
    unittest.main()
