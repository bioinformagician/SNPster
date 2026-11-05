import hashlib
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from select_unified_ancestry import SOURCES, calculate


class BenchmarkContractV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = json.loads(
            (PROJECT / "config/ancestry_benchmark_contract_v2.json").read_text()
        )

    def collapse(self, values):
        return np.array([
            sum(values.get(source, 0.0) for source in members)
            for members in self.contract["sources"]["seven_category_members"].values()
        ])

    def test_examples_match_existing_eight_source_tv_and_named_collapse(self):
        self.assertEqual(tuple(self.contract["sources"]["prediction"]), SOURCES)
        for example in self.contract["metric"]["acceptance_examples"]:
            with self.subTest(example=example["id"]):
                prediction, truth = example["prediction"], example["truth"]
                for values in (prediction, truth):
                    self.assertTrue(set(values).issubset(SOURCES))
                    self.assertAlmostEqual(sum(values.values()), 1.0)
                row = dict(sample_id="person", scenario="test", family_id="family")
                for source in SOURCES:
                    row[f"pred_{source}"] = prediction.get(source, 0.0)
                    row[f"truth_{source}"] = truth.get(source, 0.0)
                score = calculate(pd.DataFrame([row]), {})
                self.assertAlmostEqual(
                    100 * score["unstratified_mean_tv"],
                    example["expected_eight_error_pct"],
                )
                merged = 50 * np.abs(self.collapse(prediction) - self.collapse(truth)).sum()
                self.assertAlmostEqual(merged, example["expected_seven_error_pct"])

    def test_merge_is_a_partition_and_cannot_increase_tv(self):
        sources = self.contract["sources"]
        groups = sources["seven_category_members"]
        self.assertEqual(list(groups), sources["seven_category"])
        self.assertEqual(groups["EAS_SEA"], ["EAS", "SEA"])
        members = [source for group in groups.values() for source in group]
        self.assertCountEqual(members, SOURCES)
        self.assertEqual(len(members), len(set(members)))
        rng = np.random.default_rng(20261021)
        for prediction, truth in zip(rng.dirichlet(np.ones(8), 50),
                                     rng.dirichlet(np.ones(8), 50)):
            merged = np.abs(
                self.collapse(dict(zip(SOURCES, prediction)))
                - self.collapse(dict(zip(SOURCES, truth)))
            ).sum()
            self.assertLessEqual(merged, np.abs(prediction - truth).sum() + 1e-12)

    def test_same_loss_for_pure_and_mixed_without_scenario_reweighting(self):
        rows = []
        for sample, scenario, truth, prediction in [
            ("pure", "reference_source", {"EUR": 1}, {"EUR": 0.9, "EAS": 0.1}),
            ("mixed", "mixed", {"EUR": 0.5, "SSA": 0.5}, {"EUR": 0.6, "SSA": 0.4}),
        ]:
            row = dict(sample_id=sample, scenario=scenario, family_id=sample)
            for source in SOURCES:
                row[f"truth_{source}"] = truth.get(source, 0.0)
                row[f"pred_{source}"] = prediction.get(source, 0.0)
            rows.append(row)
        self.assertAlmostEqual(calculate(pd.DataFrame(rows), {})["unstratified_mean_tv"], 0.1)
        extra = dict(rows[0], sample_id="another_pure", pred_EUR=0.7, pred_EAS=0.3)
        scores = calculate(pd.DataFrame(rows + [extra]), {})
        self.assertAlmostEqual(scores["unstratified_mean_tv"], 1 / 6)
        self.assertNotAlmostEqual(scores["balanced_mean_tv"], 1 / 6)
        selection = self.contract["selection"]
        self.assertEqual(selection["person_weighting"], "equal")
        self.assertEqual(selection["headline"], "mean_person_error_pct")
        self.assertEqual(selection["order"], [
            "aggregate_chromosome_fractions_by_physical_length",
            "calculate_person_mask_tv", "average_masks_within_person",
            "average_unique_people",
        ])

    def test_metric_units_targets_and_strict_validation(self):
        metric = self.contract["metric"]
        self.assertEqual(metric["formula"],
                         "100 * 0.5 * sum(abs(prediction[source] - truth[source]))")
        self.assertEqual(metric["same_formula_for"], ["pure", "mixed"])
        self.assertEqual(metric["range"], [0, 100])
        validation = metric["vector_validation"]
        self.assertEqual(validation["sum_absolute_tolerance"], 1e-10)
        self.assertEqual(validation["sum_relative_tolerance"], 0)
        self.assertTrue(validation["finite"])
        self.assertTrue(validation["nonnegative"])
        self.assertFalse(validation["silently_normalize"])
        self.assertFalse(validation["silently_clip"])
        self.assertEqual(validation["missing_source_columns"], "error")
        self.assertEqual(self.contract["selection"]["mean_error_targets_pct"],
                         dict(acceptable_max=10, desirable_max=5))

    def test_both_resolutions_share_a_single_primary_objective(self):
        selection = self.contract["selection"]
        self.assertEqual(selection["primary_resolution"], "eight_category")
        self.assertEqual(selection["reported_resolutions"], ["eight_category", "seven_category"])
        self.assertEqual(selection["missing_or_duplicate_predictions"], "error_not_a_partial_score")
        self.assertEqual(selection["tied_primary_scores"],
                         "report_tie_no_posthoc_secondary_objective")
        self.assertIn("Seven-category baselines only",
                      self.contract["sources"]["legacy_seven_category_models"])

    def test_production_source_snapshot_has_not_changed(self):
        preprocessing = self.contract["production_preprocessing"]
        root = PROJECT.parents[2] / preprocessing["source_root_relative_to_repository"]
        self.assertEqual(preprocessing["policy"], "unchanged")
        for relative, expected in preprocessing["source_sha256"].items():
            with self.subTest(source=relative):
                self.assertEqual(hashlib.sha256((root / relative).read_bytes()).hexdigest(), expected)
        self.assertEqual(set(preprocessing["allowed_ancestry_inputs"]),
                         {"raw_harmonized", "production_phased_imputed"})
        self.assertEqual(preprocessing["not_selection_tracks"],
                         ["historical_array_phase_only", "inherited_phase"])
        self.assertIn("reference-backed imputation", preprocessing["beagle"])
        self.assertIn("do not modify that production reference",
                      self.contract["cohort_and_truth"]["reference_overlap"])

    def test_all_autosomes_and_no_automatic_promotion(self):
        self.assertEqual(self.contract["completion"]["chromosomes"], list(range(1, 23)))
        for key in ("production_pipeline_changed", "pgs_changed", "frozen_prior_protocols_changed"):
            self.assertFalse(self.contract[key])
        self.assertEqual(self.contract["comparison"]["maximum_initial_configurations"], 10)
        self.assertEqual(self.contract["comparison"]["required_comparators"],
                         ["genuine_upstream_gnomix", "best_available_phase_independent_control"])

    def test_containerized_controller_progress_and_cleanup_are_required(self):
        orchestration = self.contract["orchestration"]
        self.assertEqual(orchestration["engine"], "Nextflow DSL2")
        controller = orchestration["controller"]
        self.assertEqual(controller["runtime"], "docker")
        self.assertTrue(controller["named_container_required"])
        self.assertEqual(controller["monitor_command"], "docker logs -f <container_name>")
        self.assertIn("-ansi-log false", controller["progress"])
        self.assertIn("real task-container smoke test", controller["task_engine"])
        self.assertIn("no global pruning", orchestration["cleanup"])
        instructions = (PROJECT.parents[2] /
                        ".github/instructions/analysis-orchestration.instructions.md").read_text()
        self.assertTrue(instructions.startswith('---\ndescription: "'))
        self.assertIn('applyTo: "data_analysis/**"\n---', instructions)
        self.assertIn("Nextflow orchestrator itself in a named Docker container", instructions)
        self.assertIn("docker logs -f <container_name>", instructions)

    def test_contract_does_not_claim_an_unperformed_preflight(self):
        self.assertEqual(self.contract["version"], 2)
        self.assertEqual(self.contract["status"], "contract_only_not_execution_ready")
        self.assertEqual(set(self.contract["required_before_scientific_launch"]), {
            "resolved_production_runtime_and_reference_identity_manifest",
            "evaluation_founder_overlap_and_authoritative_eight_source_label_audit",
            "frozen_cohort_case_mix_masks_truth_and_split_manifest",
            "training_only_candidate_plan_with_at_most_ten_configurations",
            "versioned_workflow_and_independent_scorer_enforcing_this_contract",
            "named_containerized_nextflow_controller_and_real_task_container_smoke_test",
        })


if __name__ == "__main__":
    unittest.main()
