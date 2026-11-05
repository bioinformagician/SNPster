import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))

from freeze_dense_split import complete_families, local_truth, related_to
from dense_split_preprocessing import aligned_training_labels, synthesise
import dense_split_preprocessing as preprocessing
from verify_production_matched import aggregate_truth, error_pct
from finalize_production_matched import verify_references, verify
from audit_production_ancestry import digest, write_json


class DenseSplitTests(unittest.TestCase):
    def test_purge_preserves_roots_and_removes_only_detected_partners(self):
        pairs = [dict(first="a", second="b"), dict(first="b", second="c")]
        self.assertEqual(related_to(["a"], pairs), {"a", "b"})
        self.assertEqual(complete_families({"x": ["a"], "y": ["a", "b"]}, ["a"]), {"x"})

    def test_native_truth_preserves_both_distinct_founder_homologs(self):
        plan = {"founders": {"f": {"vcf_sample": "eur"}, "g": {"vcf_sample": "sea"}}}
        bp = {"f": {"h0": [(1, 20, 1)], "h1": [(1, 20, 2)]},
              "g": {"h0": [(1, 20, 3)], "h1": [(1, 20, 4)]},
              "child": {"h0": [(1, 10, 1), (11, 20, 4)], "h1": [(1, 20, 2)]}}
        segments, truth = local_truth(plan, bp, "child", {"eur": "EUR", "sea": "SEA"}, [1, 21])
        self.assertEqual(truth, [.75, 0, 0, .25, 0, 0, 0, 0])
        self.assertEqual([row["donor_haplotype"] for row in segments], [0, 1, 1])
        self.assertAlmostEqual(sum(truth), 1)
        with self.assertRaises(ValueError):
            local_truth(plan, bp, "child", {"eur": "EUR"}, [1, 21])
        bp["child"]["h0"][1] = (12, 20, 4)
        with self.assertRaises(ValueError):
            local_truth(plan, bp, "child", {"eur": "EUR", "sea": "SEA"}, [1, 21])

    def test_synthetic_pedigree_uses_both_correct_homologs(self):
        cohort = dict(local_ancestry={"22": {"person": [
            dict(haplotype=0, start=1, end=3, donor="a", donor_haplotype=1, source="EUR"),
            dict(haplotype=0, start=3, end=5, donor="b", donor_haplotype=0, source="SEA"),
            dict(haplotype=1, start=1, end=5, donor="a", donor_haplotype=0, source="EUR")]}})
        keys = [("22", i, "A", "C") for i in range(1, 5)]
        founder_calls = np.asarray([[[0, 0, 0, 0], [1, 1, 1, 1]],
                                    [[0, 1, 0, 1], [1, 0, 1, 0]]], dtype=np.int8)
        calls, labels = synthesise(cohort, 22, keys, ["a", "b"], founder_calls, [dict(sample_id="person")])
        np.testing.assert_array_equal(calls, [[[1, 1, 0, 1], [0, 0, 0, 0]]])
        np.testing.assert_array_equal(labels, [[[0, 0, 3, 3], [0, 0, 0, 0]]])
        cohort["local_ancestry"]["22"]["person"][0]["end"] = 2
        with self.assertRaisesRegex(ValueError, "exactly once"):
            synthesise(cohort, 22, keys, ["a", "b"], founder_calls, [dict(sample_id="person")])

    def test_training_only_label_alignment_tracks_flip_and_switch(self):
        truth = np.asarray([[[0, 0, 0, 1], [1, 0, 1, 0]]], dtype=np.int8)
        inferred = truth[:, ::-1].copy()
        labels = np.asarray([[[0, 0, 0, 0], [3, 3, 3, 3]]], dtype=np.int8)
        actual, evidence = aligned_training_labels(truth, inferred, labels, [1, 2, 3, 4])
        np.testing.assert_array_equal(actual, labels[:, ::-1])
        self.assertEqual(evidence, [3])
        inferred[:, :, 2:] = truth[:, :, 2:]
        actual, _ = aligned_training_labels(truth, inferred, labels, [1, 2, 3, 4])
        np.testing.assert_array_equal(actual, [[[3, 3, 0, 0], [0, 0, 3, 3]]])
        with self.assertRaisesRegex(ValueError, "informative"):
            aligned_training_labels(truth * 0, truth * 0, labels, [1, 2, 3, 4])

    def test_mixed_truth_uses_realized_physical_lengths_not_founder_average(self):
        cohort = dict(truth_design="synthetic_reference_haplotype_development",
                      truth_by_chromosome={"1": {"person": [.75, 0, 0, .25, 0, 0, 0, 0]},
                                           "22": {"person": [.25, 0, 0, .75, 0, 0, 0, 0]}})
        truth = aggregate_truth(cohort, dict(sample_id="person"), [1, 22], [3, 1])
        np.testing.assert_array_equal(truth, [.625, 0, 0, .375, 0, 0, 0, 0])
        self.assertEqual(error_pct(truth, truth), 0)
        self.assertEqual(error_pct([.5, 0, 0, .5, 0, 0, 0, 0], truth), 12.5)

    def test_raw_mask_loci_are_not_thinned_to_common_validation_truth(self):
        with tempfile.TemporaryDirectory() as directory:
            cohort = Path(directory) / "cohort.json"
            write_json(cohort, dict(input_sha256={}, chromosome_spans={"22": [1, 1000]},
                                   protocol=dict(masks=["lonni", "frederik"])))
            canonical = [("22", 100, "A", "C")]
            mask = {("22", 200, "G", "T")}
            with patch.object(preprocessing, "query_keys", return_value=canonical), \
                    patch.object(preprocessing, "mask_keys", return_value=mask), \
                    patch.object(preprocessing.base, "bind_rsids") as bind:
                preprocessing.rsids(cohort, [22], Path(directory) / "output")
            self.assertEqual(set(bind.call_args.args[2]), set(canonical) | mask)

    def test_reference_publication_is_reverified_before_and_after_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "references_smoke"
            run.mkdir()
            references = root / "references"
            folder = references / "beagle_references"
            folder.mkdir(parents=True)
            filename = preprocessing.base.reference(22).name
            target, index = folder / filename, folder / (filename + ".csi")
            target.write_bytes(b"retained reference")
            index.write_bytes(b"index")
            cohort_path = root / "cohort.json"
            write_json(cohort_path, dict(imputation_reference_ids=["fit"], imputation_excluded=["test"],
                       template_sha256={str(preprocessing.base.reference(22)): "original"}))
            (root / "source_sha256.txt").write_text(f"{digest(cohort_path)}  {cohort_path}\n")
            metadata = dict(status="all_retained_reference_alleles_and_phased_genotypes_preserved",
                            cohort_sha256=digest(cohort_path), retained_ids=["fit"], excluded_ids=["test"],
                            source_sha256="original", copy_sha256=digest(target), index_sha256=digest(index))
            path = folder / "chr22.reference.json"
            write_json(path, metadata)
            verify_references(run, references, cohort_path, [22])
            verify(run, post_cleanup=True)
            metadata["retained_ids"] = ["fit", "test"]
            write_json(path, metadata)
            with self.assertRaisesRegex(ValueError, "frozen purge"):
                verify_references(run, references, cohort_path, [22])


if __name__ == "__main__":
    unittest.main()
