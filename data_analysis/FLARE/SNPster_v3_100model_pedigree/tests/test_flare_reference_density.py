import csv
import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from screen_flare_reference_density import (
    SOURCES, candidates, chromosome_truth, command, reference_sets, thin_indices,
)
from verify_flare_reference_density import digest, parent_calls


class FlareReferenceDensityTests(unittest.TestCase):
    def test_reference_control_matches_source_counts_and_selection_is_order_invariant(self):
        rows = [dict(sample_id=f"{s}{p}_{i}", source_group=s, fine=f"{s}{p}")
                for s in SOURCES for p in range(2) for i in range(5 + p)]
        selected = reference_sets(rows, 3, 7)
        self.assertEqual(selected, reference_sets(list(reversed(rows)), 3, 7))
        lookup = {row["sample_id"]: row for row in rows}
        for source in SOURCES:
            counts = [sum(lookup[x]["source_group"] == source for x in selected[name])
                      for name in ("finecap16", "size_control")]
            self.assertEqual(counts, [6, 6])
        self.assertEqual(len(selected["full"]), len(rows))
        self.assertTrue(set(selected["finecap16"]) <= set(selected["full"]))
        self.assertTrue(set(selected["size_control"]) <= set(selected["full"]))
        self.assertNotEqual(selected["finecap16"], selected["size_control"])
        for invalid in (rows + rows[:1], [r for r in rows if r["source_group"] != "OCE"]):
            with self.assertRaises(ValueError):
                reference_sets(invalid, 3, 7)
        with self.assertRaises(ValueError):
            reference_sets(rows, 0, 7)

    def test_thinning_is_nested_exact_count_and_spans_endpoints(self):
        positions = [1, 2, 3, 4, 100, 101, 500, 501, 502, 999]
        for count in range(2, len(positions) + 1):
            indices = thin_indices(positions, count)
            self.assertEqual(len(set(indices)), count)
            self.assertEqual(sorted(indices), indices)
            self.assertEqual(indices[0], 0)
            self.assertEqual(indices[-1], len(positions) - 1)
        for positions, count in (([1, 1, 3], 2), ([1, 2], 1), ([1, 2], 3)):
            with self.assertRaises(ValueError):
                thin_indices(positions, count)

    def segments(self, directory, rows):
        path = Path(directory) / "segments.tsv"
        with path.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=(
                "sample_id", "chrom", "haplotype", "start_bp", "end_bp", "source_group"), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_truth_uses_full_physical_span_both_copies_and_specific_chromosome(self):
        rows = [dict(sample_id="child", chrom=chrom, haplotype=hap, start_bp=start,
                     end_bp=end, source_group=source)
                for chrom, hap, start, end, source in (
                    ("22", "1", 10, 30, "EUR"), ("22", "1", 30, 110, "EAS"),
                    ("22", "2", 10, 110, "SEA"), ("1", "1", 10, 110, "SSA"))]
        with tempfile.TemporaryDirectory() as directory:
            path = self.segments(directory, rows)
            actual = chromosome_truth(path, {"child"}, 22, [10, 110])["child"]
            np.testing.assert_allclose(actual, [.1, .9, 0, 0, 0, 0, 0])
            for invalid in (rows[:2], [dict(r, start_bp=11) if i == 0 else r for i, r in enumerate(rows)],
                            [dict(r, end_bp=111) if i == 2 else r for i, r in enumerate(rows)]):
                path = self.segments(directory, invalid)
                with self.assertRaises(ValueError):
                    chromosome_truth(path, {"child"}, 22, [10, 110])

    def test_native_commands_never_pass_target_truth_and_use_fixed_models(self):
        cohort = dict(flare_jar="jar", genetic_map="map", reference=dict(vcf="ref", broad_panel="panel"))
        for em in (False, True):
            result = command(cohort, "target", "prefix", em, "template.model")
            self.assertIn(f"em={str(em).lower()}", result)
            self.assertIn("update-p=false", result)
            self.assertIn("model=template.model", result)
            self.assertFalse(any(arg.startswith("gt-ancestries=") for arg in result))
        protocol = dict(reference_representations=["full", "finecap16", "size_control"], em_modes=[False, True])
        self.assertEqual(len(candidates(protocol)), 6)
        self.assertEqual(len({name for name, _, _ in candidates(protocol)}), 6)

    def test_independent_verifier_checks_named_parent_gt_and_stream_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.vcf.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("##fileformat=VCFv4.2\n")
                stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA\n")
                stream.write("22\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\n")
                stream.write("22\t20\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|0\n")
            keys = [("22", 10, "A", "G"), ("22", 20, "C", "T")]
            actual = parent_calls(path, ["A", "B"], keys)
            self.assertEqual(list(actual["A"]), [2, 0])
            self.assertEqual(list(actual["B"]), [1, 3])
            self.assertEqual(digest(path), hashlib.sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(ValueError):
                parent_calls(path, ["A"], [("22", 10, "A", "T")])


if __name__ == "__main__":
    unittest.main()
