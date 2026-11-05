import json
from pathlib import Path
import sys
import tempfile
import unittest
import gzip
from unittest.mock import patch

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from freeze_production_cohort import psam_ids
from production_matched_preprocessing import transmission, validate_concordance, write_vcf, read_call_matrix
from production_matched_models import constrained_frequencies
from verify_production_matched import error_pct, gate, validate_qp
import production_matched_preprocessing as preprocessing
from finalize_production_matched import verify as verify_published
from audit_production_ancestry import digest
import prepare_production_truth_support as truth_support


class ProductionMatchedTests(unittest.TestCase):
    def test_real_plink_psam_iid_and_fid_variants(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.psam"
            for text in ("#IID\tSEX\na\tNA\nb\tNA\n",
                         "#FID\tIID\tSEX\n0\ta\tNA\n0\tb\tNA\n"):
                path.write_text(text)
                self.assertEqual(psam_ids(path), ["a", "b"])
            path.write_text("#IID SEX\na NA\na NA\n")
            with self.assertRaisesRegex(ValueError, "duplicated"):
                psam_ids(path)
            path.write_text("#FID SEX\n0 NA\n")
            with self.assertRaisesRegex(ValueError, "IID"):
                psam_ids(path)

    def test_model_plan_has_exactly_ten_and_both_resolutions(self):
        protocol = json.loads((PROJECT / "config/production_matched_protocol.json").read_text())
        contract = json.loads((PROJECT / "config/ancestry_benchmark_contract_v2.json").read_text())
        self.assertEqual(protocol["sources"], contract["sources"]["prediction"])
        self.assertEqual(len(protocol["candidates"]), 10)
        self.assertEqual(sum(row["family"] == "flare" for row in protocol["candidates"].values()), 6)
        self.assertEqual(sum(row["family"] == "gnomix" for row in protocol["candidates"].values()), 2)
        self.assertEqual(sum(row["family"] == "phasefree" for row in protocol["candidates"].values()), 2)
        self.assertEqual(len(protocol["masks"]), 2)
        self.assertTrue(protocol["production"]["unchanged"])
        self.assertFalse(protocol["production_promotion"])
        self.assertFalse(protocol["independent_final_test"])

    def test_same_tv_for_pure_mixed_and_eas_sea_merge(self):
        pure = np.eye(8)[0]
        self.assertAlmostEqual(error_pct(.9 * pure + .1 * np.eye(8)[1], pure), 10)
        truth = (np.eye(8)[0] + np.eye(8)[1] + np.eye(8)[3]) / 3
        prediction = .4 * np.eye(8)[1] + .4 * np.eye(8)[3] + .2 * np.eye(8)[7]
        self.assertAlmostEqual(error_pct(prediction, truth), 100 / 3)
        self.assertEqual(error_pct(np.eye(8)[1], np.eye(8)[3], collapse=True), 0)
        self.assertEqual(error_pct(np.eye(8)[1], np.eye(8)[3]), 100)
        with self.assertRaisesRegex(ValueError, "Invalid"):
            error_pct(np.ones(8), pure)

    def test_per_person_untyped_gate_cannot_hide_bad_individual(self):
        policy = dict(minimum_untyped_validation_markers=100, minimum_observed_genotype_concordance=.995,
                      minimum_untyped_genotype_concordance=.9)
        truth = np.zeros((32, 2000), dtype=int)
        actual = truth.copy()
        observed = np.arange(2000) < 1000
        self.assertEqual(validate_concordance(truth, actual, observed, policy)["untyped"]["concordance"], 1)
        actual[0, 1000:1200] = 1
        with self.assertRaisesRegex(ValueError, "minimum-person"):
            validate_concordance(truth, actual, observed, policy)
        actual[0, 0] = -1
        with self.assertRaisesRegex(ValueError, "missing"):
            validate_concordance(truth, actual, observed, policy)

    def test_real_vcf_round_trip_does_not_confuse_unphased_gt_with_phase(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            keys = [("22", 100, "A", "C"), ("22", 200, "G", "T")]
            paired = np.asarray([[[0, 1], [1, 0]]], dtype=np.int8)
            path = write_vcf(root / "observed.vcf", keys, ["person"], paired)
            np.testing.assert_array_equal(read_call_matrix(path, keys, ["person"], root), paired)
            with self.assertRaisesRegex(ValueError, "unphased"):
                read_call_matrix(path, keys, ["person"], root, phased=True)
            path = write_vcf(root / "phased.vcf", keys, ["person"], paired, phased=True)
            np.testing.assert_array_equal(read_call_matrix(path, keys, ["person"], root, phased=True), paired)
            self.assertIn(transmission(20261022, 22, "f1", "parent"), (0, 1))
            self.assertEqual(transmission(20261022, 22, "f1", "parent"),
                             transmission(20261022, 22, "f1", "parent"))

    def test_phasefree_simplex_fit_and_independent_convex_optimality(self):
        rng = np.random.default_rng(20261022)
        frequencies = rng.uniform(.05, .95, size=(2000, 8))
        truth = np.eye(8)[2]
        dosage = rng.binomial(2, frequencies @ truth)
        result = constrained_frequencies(frequencies, dosage)
        self.assertGreater(result[2], .85)
        self.assertAlmostEqual(result.sum(), 1)
        validate_qp(frequencies, dosage, result)
        with self.assertRaisesRegex(ValueError, "optimality"):
            validate_qp(frequencies, dosage, np.eye(8)[0])
        with self.assertRaisesRegex(ValueError, "at least100"):
            constrained_frequencies(frequencies[:20], dosage[:20])

    def test_preprocessing_gate_rejects_missing_and_duplicate_batches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cohort = root / "cohort.json"
            cohort.write_text(json.dumps({"people": [{}] * 32}))
            report = root / "report.json"
            report.write_text(json.dumps(dict(chromosome=22, mask="lonni", people=32, production_changed=False,
                                               status="production_preprocessing_genotypes_and_coverage_verified",
                                               input_sha256={})))
            gate(cohort, [22], ["lonni"], [report], root / "passed.json")
            with self.assertRaisesRegex(ValueError, "All declared"):
                gate(cohort, [22], ["lonni", "frederik"], [report], root / "passed.json")
            with self.assertRaisesRegex(ValueError, "duplicated"):
                gate(cohort, [22], ["lonni"], [report, report], root / "passed.json")

    def test_missing_input_alleles_are_explicitly_intersected_never_filled(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cohort = root / "cohort.json"
            candidates = [dict(sample_id=f"evaluation{i}", public=i == 5) for i in range(8)]
            training = [dict(sample_id=f"training{i}") for i in range(16)]
            public = root / "public.vcf.gz"
            cohort.write_text(json.dumps(dict(evaluation_founders=candidates, training_founders=training,
                                               normalized_vcfs={"22": str(public)})))
            keys = [("22", 100 + i, "A", "C") for i in range(120)]
            markers = root / "markers.tsv"
            markers.write_text("".join(f"{c}\t{p}\t{r}\t{a}\trs{p}\n" for c, p, r, a in keys))
            def available(path, chromosome, region):
                return keys[:-1] if path == public else keys
            def calls(path, selected, ids, output):
                return np.zeros((len(ids), 2, len(selected)), dtype=np.int8)
            with patch.object(preprocessing, "PROJECT", root), \
                    patch.object(preprocessing, "query_keys", side_effect=available), \
                    patch.object(preprocessing, "read_call_matrix", side_effect=calls):
                preprocessing.founders(cohort, 22, markers, root / "founders")
            record = json.loads((root / "founders/founders.json").read_text())
            self.assertEqual(record["source_loci"], 120)
            self.assertEqual(record["retained_loci"], 119)
            self.assertEqual(record["unavailable_exact_loci_by_input"][str(public)], 1)
            self.assertEqual(np.load(root / "founders/observed.npz")["calls"].shape, (24, 2, 119))

    def test_finalizer_binds_independently_checked_published_models_before_and_after_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "smoke"
            for name in ("preparation", "founders", "founder_phase", "truth", "raw", "standardized",
                         "combined", "harmonized", "merged", "imputed", "split", "qc", "indexed_qc",
                         "preprocessing", "preprocessing_gate", "panels", "models", "summary"):
                (run / name).mkdir(parents=True)
            model = run / "models/known_model.tsv"
            model.write_text("independently audited prediction\n")
            independent = run / "summary/summary/independent_verification.json"
            independent.parent.mkdir()
            independent.write_text(json.dumps(dict(
                status="smoke_passed_not_a_whole_genome_leaderboard",
                independently_verified_model_sha256={str(model): digest(model)})))
            source = root / "source.txt"
            source.write_text("frozen source\n")
            (root / "source_sha256.txt").write_text(f"{digest(source)}  {source}\n")
            model.write_text("changed before durable verification\n")
            with self.assertRaises(ValueError):
                verify_published(run)
            model.write_text("independently audited prediction\n")
            verify_published(run)
            verify_published(run, post_cleanup=True)
            model.write_text("changed after cleanup\n")
            with self.assertRaises(ValueError):
                verify_published(run, post_cleanup=True)

    def test_original_assay_orientation_excludes_evaluation_donor_and_rejects_gt_changes(self):
        original = np.array([0, 0, 2, 2, 1, 0], dtype=np.int8)
        reverse = 2 - original
        self.assertEqual(truth_support.reconstruct_assay(original, reverse, 5), 2)
        changed = reverse.copy()
        changed[5] = 1
        with self.assertRaisesRegex(ValueError, "differs"):
            truth_support.reconstruct_assay(original, changed, 5)
        changed[0] = 1
        with self.assertRaisesRegex(ValueError, "disagree"):
            truth_support.reconstruct_assay(original, changed, 5)
        ambiguous = np.array([1, 1, 1, 1, 1, 0], dtype=np.int8)
        self.assertIsNone(truth_support.reconstruct_assay(ambiguous, ambiguous, 5))
        ambiguous[5] = 1
        self.assertEqual(truth_support.reconstruct_assay(ambiguous, ambiguous, 5), 1)
        changed_heterozygote = ambiguous.copy()
        changed_heterozygote[5] = 0
        with self.assertRaisesRegex(ValueError, "differs"):
            truth_support.reconstruct_assay(ambiguous, changed_heterozygote, 5)

    def test_completeness_uses_post_assay_qc_without_filling_missing_calls(self):
        calls = np.zeros((24, 2, 1000), dtype=np.int8)
        calls[3, :, :130] = -1
        calls[0, :, 130:150] = -1
        eligible = np.ones(1000, dtype=bool)
        eligible[:130] = False
        before = calls.copy()
        retained, report = truth_support.observed_completeness(calls, eligible)
        self.assertEqual(report["source_compatible_loci"], 1000)
        self.assertEqual(report["assay_eligible_loci"], 870)
        self.assertEqual(report["complete_loci"], 850)
        self.assertAlmostEqual(report["complete_fraction_of_eligible"], 850 / 870)
        self.assertAlmostEqual(report["complete_fraction_of_source_compatible"], .85)
        self.assertFalse(retained[:150].any())
        np.testing.assert_array_equal(calls, before)
        calls[0, :, 150:300] = -1
        _, bad = truth_support.observed_completeness(calls, eligible)
        self.assertLess(bad["complete_fraction_of_eligible"], .9)
        with self.assertRaisesRegex(ValueError, "No loci"):
            truth_support.observed_completeness(calls, np.zeros(1000, dtype=bool))
        with self.assertRaisesRegex(ValueError, "outside"):
            truth_support.observed_completeness(
                np.zeros((24, 2, 1000), dtype=np.int8), eligible)

    def test_original_assay_quality_filter_never_uses_filled_genotypes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "original.txt.gz"
            names = ["Bajo25KN"] + [f"anchor{i}" for i in range(195)]
            header = ["Name"] + [field for name in names for field in (name + ".GType", name + ".Score")]
            values = ["rs100"] + [value for i in range(196) for value in ("AA", ".9" if i else ".7")]
            with gzip.open(raw, "wt") as handle:
                handle.write("\t".join(header) + "\n" + "\t".join(values) + "\n")
            with patch.object(truth_support, "RAW", raw):
                truth_support.assay(root / "assay")
            axis, calls = truth_support.raw_dosage(root / "assay/assay.tsv.gz", {"rs100"})
            self.assertEqual(axis, names)
            self.assertEqual(int(calls["rs100"][0]), -1)
            self.assertTrue((calls["rs100"][1:] == 0).all())

if __name__ == "__main__":
    unittest.main()
