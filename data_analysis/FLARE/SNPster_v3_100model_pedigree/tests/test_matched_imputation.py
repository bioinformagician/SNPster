from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from pilot_matched_imputation import ancestor_donors, genotype_strata, pedigree_haplotypes, phase_metrics


class MatchedImputationTests(unittest.TestCase):
    def test_ancestors_follow_the_selected_node_not_unrelated_family_members(self):
        plan = {
            "nodes": [dict(sample_id="child", father_id="a", mother_id="b")],
            "founders": {person: dict(vcf_sample=donor)
                         for person, donor in (("a", "A"), ("b", "B"), ("other", "UNUSED"))},
        }
        self.assertEqual(ancestor_donors(plan, "child"), {"A", "B"})

    def test_breakpoint_reconstruction_preserves_haplotypes_and_inclusive_ends(self):
        plan = {"founders": {"a": {"vcf_sample": "A"}, "b": {"vcf_sample": "B"}}}
        bp = {
            "a": {"h0": [(10, 40, 1)], "h1": [(10, 40, 2)]},
            "b": {"h0": [(10, 40, 3)], "h1": [(10, 40, 4)]},
            "child": {"h0": [(10, 20, 1), (21, 40, 2)], "h1": [(10, 40, 4)]},
        }
        donors = np.array([[0, 0, 0, 0], [1, 1, 1, 1],
                           [0, 1, 0, 1], [1, 0, 1, 0]])
        got = pedigree_haplotypes(plan, bp, "child", ["A", "B"], donors,
                                  np.array([10, 20, 21, 40]))
        np.testing.assert_array_equal(got, [[0, 0, 1, 1], [1, 0, 1, 0]])
        with self.assertRaisesRegex(ValueError, "cover"):
            pedigree_haplotypes(plan, bp, "child", ["A", "B"], donors,
                                np.array([9, 20, 21, 40]))

    def test_switch_rate_allows_whole_chromosome_swap(self):
        truth = np.array([[0, 0, 1, 0, 1], [1, 1, 0, 1, 0]])
        result = phase_metrics(truth, truth[::-1])
        self.assertEqual(result["switches"], 0)
        self.assertEqual(result["phase_hamming_errors"], 0)
        self.assertEqual(result["switch_pairs"], 4)
        self.assertEqual(result["n_genotype_errors"], 0)
        altered = truth.copy()
        altered[:, 2:] = altered[::-1, 2:]
        result = phase_metrics(truth, altered)
        self.assertEqual(result["switches"], 1)
        self.assertEqual(result["phase_hamming_errors"], 2)

    def test_phase_denominator_excludes_genotype_errors_and_single_het_is_unevaluable(self):
        truth = np.array([[0, 0, 1], [1, 1, 0]])
        inferred = np.array([[0, 0, 0], [1, 0, 0]])
        result = phase_metrics(truth, inferred)
        self.assertEqual(result["n_genotype_errors"], 2)
        self.assertEqual(result["n_correct_heterozygotes"], 1)
        self.assertEqual(result["switch_pairs"], 0)
        with self.assertRaises(ValueError):
            phase_metrics(truth, np.full_like(truth, 2))

    def test_observed_and_imputed_metrics_are_separate(self):
        truth = np.array([[0, 0, 1], [1, 1, 0]])
        inferred = truth.copy()
        inferred[:, 1] = 0
        result = genotype_strata(truth, inferred, np.array([True, False, True]))
        self.assertEqual(result["observed"]["n_genotype_errors"], 0)
        self.assertEqual(result["imputed"]["n_genotype_errors"], 1)
        self.assertEqual(result["all"]["n_sites"], 3)


if __name__ == "__main__":
    unittest.main()
