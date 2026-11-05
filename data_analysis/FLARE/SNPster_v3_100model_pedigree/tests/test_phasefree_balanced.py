import json
from pathlib import Path
import pickle
import sys
import unittest

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import sweep_phasefree_balanced as balanced
import sweep_phasefree_rich as rich


class PhasefreeBalancedTests(unittest.TestCase):
    def fixture(self):
        labels = pd.DataFrame([
            dict(sample_id=f"{source}{i}", source_group=source, fine=f"pop{source}{i}")
            for source in balanced.SOURCES for i in range(2)])
        settings = json.loads(balanced.PROTOCOL.read_text())
        spans = {"1": [10, 110], "2": [0, 200]}
        keys = [("1", pos, "A", "G") for pos in (10, 11, 12, 13, 100, 109)]
        keys += [("2", pos, "A", "G") for pos in (1, 100, 199)]
        founders = np.random.default_rng(1).integers(0, 2, (len(labels), 2, len(keys)), dtype=np.uint8)
        return labels, settings, spans, keys, founders

    def test_plan_bank_balances_roles_fractions_and_repairs_donor_safe_coverage(self):
        labels, settings, _, _, _ = self.fixture()
        specs = balanced.templates(settings)
        self.assertEqual(len(specs), 371)
        self.assertEqual(sum(row["replicates"] for row in specs), 1428)
        splits = [dict(heldout_donors=labels[labels.sample_id.str.endswith("0")].sample_id.tolist())]
        plans, initial = balanced.make_plans(labels, splits, settings)
        self.assertEqual(initial, 1428)
        self.assertLessEqual(len(plans), settings["maximum_mosaics"])
        kept = balanced.safe_indices(plans, splits[0]["heldout_donors"])
        for template in specs:
            self.assertGreaterEqual(sum(plans[i]["stratum"] == template["stratum"] for i in kept), 2)
        for plan in plans:
            self.assertEqual(len(plan["donors"]), len(set(plan["donors"])))
            self.assertTrue(set(plan["donors"]) <= set(labels.sample_id))
        settings["maximum_mosaics"] = 10
        with self.assertRaisesRegex(ValueError, "fixed mosaic budget"):
            balanced.make_plans(labels, splits, settings)

    def test_physical_truth_does_not_follow_unequal_snp_density(self):
        labels, _, spans, keys, founders = self.fixture()
        eur = labels.index[labels.sample_id.eq("EUR0")][0]
        sas = labels.index[labels.sample_id.eq("SAS0")][0]
        founders[eur] = 0
        founders[sas] = 1
        segments = np.array([
            [1, 0, 10, 110, eur, 0, 0], [1, 1, 10, 110, eur, 1, 0],
            [2, 0, 0, 200, eur, 0, 0], [2, 1, 0, 200, sas, 1, 2]])
        values, truth = balanced.render(segments, founders, keys, spans)
        np.testing.assert_allclose(truth, [2 / 3, 0, 1 / 3, 0, 0, 0, 0])
        self.assertAlmostEqual(values.sum() / (2 * len(keys)), 1 / 6)
        self.assertNotAlmostEqual(truth[2], values.sum() / (2 * len(keys)))

    def test_actual_founder_haplotypes_are_inherited_and_target_bp_is_exact(self):
        labels, settings, spans, keys, founders = self.fixture()
        rng = np.random.default_rng(3)
        for template in balanced.templates(settings)[::37]:
            plan = balanced.new_plan(template, labels, rng, 1)
            segments = balanced.segments_for(plan, labels, spans, 30)
            values, truth = balanced.render(segments, founders, keys, spans)
            self.assertTrue(np.isin(values, [0, 1, 2]).all())
            target = np.zeros(7)
            for source, weight in zip(plan["sources"], plan["target_weights"]):
                target[balanced.SOURCES.index(source)] = weight
            np.testing.assert_allclose(truth, target, atol=2 / 600, rtol=0)
            if plan["scenario"] == "f1":
                self.assertEqual(len(segments), 4)
                np.testing.assert_allclose(sorted(truth[truth > 0]), [.5, .5])

    def test_bad_donor_labels_gaps_and_unknown_chromosomes_fail(self):
        labels, settings, spans, keys, founders = self.fixture()
        template = balanced.templates(settings)[0]
        plan = balanced.new_plan(template, labels, np.random.default_rng(1), 0)
        plan["donors"][0] = "not_training"
        with self.assertRaisesRegex(ValueError, "identity/source"):
            balanced.segments_for(plan, labels, spans, 30)
        plan = balanced.new_plan(template, labels, np.random.default_rng(1), 0)
        segments = balanced.segments_for(plan, labels, spans, 30)
        bad = segments.copy()
        bad[0, 2] += 1
        with self.assertRaisesRegex(ValueError, "gap"):
            balanced.render(bad, founders, keys, spans)
        bad = segments.copy()
        bad[0, 0] = 3
        with self.assertRaisesRegex(ValueError, "chromosome"):
            balanced.render(bad, founders, keys, spans)

    def test_same_founder_on_both_copies_preserves_original_diploid_genotype(self):
        labels, settings, spans, keys, founders = self.fixture()
        founders[:, 0] = 0
        founders[:, 1] = 1
        template = next(row for row in balanced.templates(settings) if row["scenario"] == "minor")
        plan = balanced.new_plan(template, labels, np.random.default_rng(4), 0)
        segments = balanced.segments_for(plan, labels, spans, 30)
        values, _ = balanced.render(segments, founders, keys, spans)
        selected = np.full((2, len(keys)), -1, dtype=int)
        for chrom, hap, start, end, donor, _, _ in segments:
            for i, key in enumerate(keys):
                if key[0] == str(chrom) and start <= key[1] < end:
                    selected[hap, i] = donor
        same = selected[0] == selected[1]
        self.assertTrue(same.any())
        np.testing.assert_array_equal(values[same], np.ones(same.sum(), dtype=np.uint8))

    def test_augmentation_weights_balance_scenarios_and_strata_not_rows(self):
        _, settings, _, _, _ = self.fixture()
        plans = []
        for scenario in settings["profiles"]["balanced"]:
            for stratum in range(2):
                plans += [dict(sample_id=f"{scenario}{stratum}_{i}", scenario=scenario,
                               stratum=f"{scenario}:{stratum}") for i in range(stratum + 1)]
        for profile in settings["profiles"].values():
            weights = balanced.augmentation_weights(plans, profile)
            frame = pd.DataFrame(plans)
            for scenario, expected in profile.items():
                self.assertAlmostEqual(weights[frame.scenario.eq(scenario)].sum(), expected)
                for stratum in range(2):
                    self.assertAlmostEqual(weights[frame.stratum.eq(f"{scenario}:{stratum}")].sum(), expected / 2)

    def test_all_weight_policies_roundtrip_and_only_control_omits_augmentation(self):
        _, settings, _, _, _ = self.fixture()
        protocol = {**json.loads(rich.PROTOCOL.read_text()), "augmentation": settings,
                    "models": settings["models"], "seed": settings["decoder_seed"]}
        protocol["boosting"]["iterations"] = 3
        protocol = json.loads(json.dumps(protocol, sort_keys=True))
        random = np.random.default_rng(2)
        x = np.column_stack([
            random.dirichlet(np.ones(7), 72), random.dirichlet(np.ones(7), 72),
            random.normal(size=(72, 64)), random.dirichlet(np.ones(100), 72)])
        truth = random.dirichlet(np.ones(7), 72)
        scenarios = list(settings["profiles"]["balanced"])
        plans = [dict(sample_id=f"aug{i}", scenario=scenarios[i % 5], stratum=scenarios[i % 5])
                 for i in range(30)]
        with threadpool_limits(limits=2):
            models = balanced.fit_models(x[:42], truth[:42], np.repeat(1 / 42, 42),
                                         x[42:], truth[42:], plans, protocol)
            for name, model in models.items():
                self.assertEqual(model["n_augmented"], 0 if name == balanced.CONTROL else 30)
                values = balanced.predict(model, x[:5])
                np.testing.assert_array_equal(values, balanced.predict(pickle.loads(pickle.dumps(model)), x[:5]))
                np.testing.assert_array_equal(values, balanced.predict(model, np.nextafter(x[:5], np.inf)))
                np.testing.assert_allclose(values.sum(axis=1), 1, atol=1e-12, rtol=0)
                self.assertTrue(np.all(values >= 0))


if __name__ == "__main__":
    unittest.main()
