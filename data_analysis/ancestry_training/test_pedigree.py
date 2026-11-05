import copy
import unittest
from collections import Counter

from prepare_pedigree import allocate, build_trees, pair_generation
from evaluation import person_metrics


class PedigreeDatasetTests(unittest.TestCase):
    def fixture(self):
        sources = ["A", "B"]
        nodes = {}
        for source in sources:
            for index in range(40):
                identity = f"{source}{index}"
                nodes[identity] = {
                    "sample_id": identity, "kind": "real_person", "source": source,
                    "genotype_available": True, "kinship_component": identity,
                    "expected_ancestry": {s: float(s == source) for s in sources},
                }
        spec = {"version": 1, "seed": 1234, "generations": 3, "minimum_reference_per_source": 5,
                "founders": {role: {"A": 8, "B": 8}
                             for role in ("train", "calibration", "validation", "test")}}
        return nodes, spec, sources

    def test_whole_components_and_reference_exclusion(self):
        nodes, spec, sources = self.fixture()
        nodes["relative"] = {**nodes["A0"], "sample_id": "relative"}
        selected, reference = allocate(nodes, spec, sources)
        self.assertEqual((selected, reference), allocate(nodes, spec, sources))
        active = reference + [identity for values in selected.values() for identity in values]
        self.assertEqual(len(active), len({nodes[identity]["kinship_component"] for identity in active}))
        self.assertTrue(all(len(values) == 16 for values in selected.values()))

    def test_no_reuse_three_generations_and_expected_truth(self):
        nodes, spec, sources = self.fixture()
        selected, _ = allocate(nodes, spec, sources)
        result, trees = build_trees(selected, nodes, sources, spec["seed"])
        for role, roots in trees.items():
            self.assertEqual(len(roots), 2)
            generations = Counter(node["generation"] for node in result.values() if node.get("role") == role)
            self.assertEqual(generations, {0: 16, 1: 8, 2: 4, 3: 2})
            used = set()
            for identity in roots:
                donors = set(result[identity]["founder_donors"])
                self.assertEqual(len(donors), 8)
                self.assertFalse(used & donors)
                used |= donors
            self.assertEqual(used, set(selected[role]))
            for node in result.values():
                if node.get("role") == role:
                    self.assertEqual(sum(node["expected_ancestry"].values()), 1)
        self.assertEqual(nodes["A0"].get("generation"), None)

    def test_insufficient_founders_and_reference_fail(self):
        nodes, spec, sources = self.fixture()
        for node in nodes.values():
            node["kinship_component"] = node["source"]
        with self.assertRaises(ValueError):
            allocate(nodes, spec, sources)
        nodes, spec, sources = self.fixture()
        spec["minimum_reference_per_source"] = 9
        with self.assertRaises(ValueError):
            allocate(nodes, spec, sources)

    def test_nonhalving_and_invalid_spec_fail(self):
        nodes, spec, sources = self.fixture()
        bad = copy.deepcopy(spec)
        bad["founders"]["test"]["A"] += 1
        with self.assertRaises(ValueError):
            allocate(nodes, bad, sources)
        bad = copy.deepcopy(spec)
        bad["generations"] = 4
        with self.assertRaises(ValueError):
            allocate(nodes, bad, sources)

    def test_source_pure_first_generation_fails(self):
        nodes, _, _ = self.fixture()
        with self.assertRaises(ValueError):
            pair_generation(["A0", "A1"], nodes, 123)

    def test_generation_precision_keeps_people_and_family_units_separate(self):
        people = {
            "pure": [{"family_id": "tree1", "track": "pure", "sources": ["A", "B"],
                      "predictions": [.9, .1], "truth": [1., 0.], "tv_percentage_points": 10.}],
            "child": [{"family_id": "tree1", "track": "mixed_g1", "sources": ["A", "B"],
                       "predictions": [.5, .5], "truth": [.5, .5], "tv_percentage_points": 0.}],
            "other": [{"family_id": "tree2", "track": "mixed_g3", "sources": ["A", "B"],
                       "predictions": [.75, .25], "truth": [.5, .5], "tv_percentage_points": 25.}],
        }
        metrics = person_metrics(people, ["A", "B"], 1234)
        self.assertEqual(metrics, person_metrics(people, ["A", "B"], 1234))
        self.assertEqual(metrics["people"], 3)
        self.assertEqual(metrics["family_components"], 2)
        self.assertAlmostEqual(metrics["mean_person_tv_pp"], 35 / 3)
        self.assertAlmostEqual(metrics["per_source"]["A"]["mae_pp"], 35 / 3)
        self.assertEqual(metrics["by_track"]["mixed_g3"]["people"], 1)
        self.assertIsNone(person_metrics({"pure": people["pure"]}, ["A", "B"], 1234)[
            "family_bootstrap_mean_tv_95pct_ci_pp"])


if __name__ == "__main__":
    unittest.main()
