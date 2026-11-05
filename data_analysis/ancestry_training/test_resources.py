import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import core
from master_resources import assign_roles, base_panel, canonical, project_chromosome, read_gt, write_gt
from prepare_kinship import SAMPLE_FIELDS, release
from prepare_observed import write_rows
from prepare_pedigree import pair_generation


class ResourceTests(unittest.TestCase):
    def nodes(self):
        result = {}
        for source in ("A", "B"):
            for index in range(12):
                identity = f"{source}{index}"
                result[identity] = {"sample_id": identity, "resource": "hgdp1kg", "source": source, "label_role": "source"}
        result["relative"] = {"sample_id": "relative", "resource": "gse73996", "source": "A", "label_role": "source"}
        result["unknown"] = {"sample_id": "unknown", "resource": "kazakh", "source": "", "label_role": "affinity_only",
                             "expected_ancestry": None}
        return result

    def design(self):
        return {"founders": {role: {"A": 1, "B": 1} for role in ("train", "calibration", "validation", "test")}}

    def test_canonical_ids_do_not_accept_ambiguous_vcf_delimiters(self):
        self.assertEqual(canonical("lebbo 3"), "lebbo_3")
        for identity in ("a,b", "a/b", "a:b", ""):
            with self.assertRaises(ValueError):
                canonical(identity)

    def test_true_relatives_share_role_and_are_retained(self):
        nodes = self.nodes()
        groups = [["A0", "relative"]] + [[identity] for identity in nodes if identity not in {"A0", "relative"}]
        assign_roles(nodes, groups, self.design(), 123)
        self.assertEqual(nodes["A0"]["component_role"], nodes["relative"]["component_role"])
        self.assertIn(nodes["relative"]["role"], core.ROLES)
        self.assertEqual(nodes["A0"]["kinship_component"], nodes["relative"]["kinship_component"])

    def test_affinity_is_not_coerced_to_source_truth(self):
        nodes = self.nodes()
        groups = [[identity] for identity in nodes]
        assign_roles(nodes, groups, self.design(), 123)
        self.assertEqual(nodes["unknown"]["role"], "affinity_only")
        self.assertIsNone(nodes["unknown"]["expected_ancestry"])

    def test_global_split_is_deterministic(self):
        first = self.nodes()
        second = copy.deepcopy(first)
        groups = [[identity] for identity in first]
        a = assign_roles(first, groups, self.design(), 123)
        b = assign_roles(second, groups, self.design(), 123)
        self.assertEqual(a, b)
        self.assertEqual(first, second)
        self.assertEqual(len({person for people in a[0].values() for person in people}), 8)

    def test_founder_pairing_keeps_a_feasible_cross_source_remainder(self):
        nodes = {identity: {"source": source, "expected_ancestry": {
            label: float(label == source) for label in ("A", "B", "C")}}
                 for source, people in (("A", ("a1", "a2", "a3")), ("B", ("b1", "b2")), ("C", ("c1",)))
                 for identity in people}
        for seed in range(100):
            pairs = pair_generation(list(nodes), nodes, seed)
            self.assertEqual(len(pairs), 3)
            self.assertEqual({person for pair in pairs for person in pair}, set(nodes))
            self.assertTrue(all(nodes[a]["source"] != nodes[b]["source"] for a, b in pairs))

    def test_missing_source_calls_roundtrip_without_filling(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.vcf.gz"
            keys = [("22", 100, "A", "G"), ("22", 200, "C", "T")]
            calls = np.array([[[0, -1], [1, -1]], [[-1, 1], [0, 1]]], dtype=np.int8)
            write_gt(path, keys, ["a", "b"], calls)
            actual_keys, actual_calls = read_gt(path, 22, ["a", "b"])
            self.assertEqual(actual_keys, keys)
            np.testing.assert_array_equal(actual_calls, calls)

    def test_every_source_gt_and_separator_roundtrips(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            keys = [("22", pos, "A", "G") for pos in (100, 200, 300)]
            pairs = np.asarray([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], dtype=np.int8)
            calls = pairs.reshape(3, 3, 2).transpose(0, 2, 1)
            for phased in (False, True):
                path = root / f"source{phased}.vcf.gz"
                write_gt(path, keys, ["a", "b", "c"], calls, phased=phased)
                actual_keys, actual_calls = read_gt(path, 22, ["a", "b", "c"])
                self.assertEqual(actual_keys, keys)
                np.testing.assert_array_equal(actual_calls, calls)

    def test_model_writer_rejects_missing_or_noninteger_gt(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.vcf.gz"
            for calls in (np.array([[[-1], [0]]], dtype=np.int8), np.array([[[0.], [1.]]])):
                with self.assertRaisesRegex(ValueError, "complete biallelic"):
                    core.write_vcf(path, [("22", 100, "A", "G")], ["a"], calls)

    def test_gt_decoder_preserves_phase_and_missingness_contracts(self):
        np.testing.assert_array_equal(core.decode_gt("0|1\t1|0", 2, phased=True), [[0, 1], [1, 0]])
        np.testing.assert_array_equal(core.decode_gt(".\t0/.", 2, allow_missing=True), [[-1, -1], [0, -1]])
        for value in ("0/.", ".", "0/2", "0:1", "0/1\t0/1"):
            with self.assertRaises(ValueError):
                core.decode_gt(value, 1)
        with self.assertRaisesRegex(ValueError, "Unphased"):
            core.decode_gt("0/1", 1, phased=True)
        self.assertFalse(core.decode_gt("0/1", 1).flags.writeable)

    def test_projection_output_and_temporary_directories_are_distinct(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = {"specification": {"observed_dataset": str(root)}, "resources": {}}
            with patch("master_resources.core.read", return_value=data), \
                    patch("master_resources.sites", side_effect=RuntimeError("reached source reader")):
                with self.assertRaisesRegex(RuntimeError, "reached source reader"):
                    project_chromosome(root, 22, root / "projection")
            self.assertTrue((root / "projection").is_dir())
            self.assertTrue((root / "projection_work").is_dir())

    def test_custom_assay_mask_reaches_query_export(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            master, cohort = root / "master", root / "cohort"
            master.mkdir()
            cohort.mkdir()
            people = ["ref", "train", "cal", "val", "test"]
            keys = [("22", 100, "A", "G"), ("22", 200, "C", "T")]
            core.write_vcf(master / "all.phased.vcf.gz", keys, people, np.zeros((5, 2, 2), dtype=np.int8))
            write_rows(master / "assay.sites.tsv", [{"chromosome": "22", "position": 100, "ref": "A", "alt": "G"}],
                       ["chromosome", "position", "ref", "alt"])
            write_rows(cohort / "samples.tsv", [
                dict(zip(SAMPLE_FIELDS, (person, role, person, person, "A", "pure", "population")))
                for person, role in zip(people, core.ROLES)], SAMPLE_FIELDS)
            release(master, cohort, 22, root / "bundle", masks=["assay"])
            actual, calls = core.vcf(root / "bundle/query_assay.vcf.gz", 22, people[2:], False)
            self.assertEqual(actual, keys[:1])
            self.assertEqual(calls.shape, (3, 2, 1))

    def test_external_reference_role_does_not_overwrite_phase_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            projection, cohort, inventory = (root / name for name in ("projection", "cohort", "inventory"))
            for folder in (projection, cohort, inventory):
                folder.mkdir()
            keys = [("22", 100, "A", "G"), ("22", 200, "C", "T")]
            external = [f"external_{role}" for role in core.ROLES]
            core.write_vcf(projection / "hgdp1kg.vcf.gz", keys, ["hgdp_reference"], np.zeros((1, 2, 2), dtype=np.int8))
            core.write_vcf(projection / "gse73996.vcf.gz", keys, external, np.ones((5, 2, 2), dtype=np.int8), phased=False)
            people = ["hgdp_reference"] + external
            write_rows(cohort / "samples.tsv", [
                {"sample_id": person, "role": role, "family_id": person, "donors": person,
                 "source": "A", "track": "pure", "fine": "population"}
                for person, role in zip(people, ("reference", *core.ROLES))], SAMPLE_FIELDS)
            core.write(cohort / "plan.json", {
                "panel": "gse73996", "masks": ["gse73996"], "spans": {"22": [1, 300]},
                "nodes": {person: {"kind": "real_person", "resource": "hgdp1kg" if person == "hgdp_reference" else "gse73996"}
                          for person in people},
                "template": {"tools": {"beagle_jar": "fixture.jar"}, "maps": {"22": "fixture.map"}},
            })
            core.write(inventory / "inventory.json", {"specification": {"minimum_pipeline_markers": 1, "seed": 123}})
            original_run = core.run
            invocations = []

            def phase_fixture(command, prefix):
                if command[0] != "java":
                    return original_run(command, prefix)
                options = dict(argument.split("=", 1) for argument in command if "=" in str(argument))
                self.assertNotEqual(options["gt"], options["ref"])
                self.assertEqual(core.rows(cohort / "samples.tsv")[0]["sample_id"], "hgdp_reference")
                ref_keys, ref_calls = core.vcf(options["ref"], 22, ["hgdp_reference"], True)
                np.testing.assert_array_equal(ref_calls, np.zeros((1, 2, 2), dtype=np.int8))
                ids = [Path(options["gt"]).stem.split(".")[0]]
                role = ids[0]
                actual_keys, calls = core.vcf(options["gt"], 22, [f"external_{role}"], False)
                self.assertEqual(actual_keys, ref_keys)
                core.write_vcf(Path(options["out"] + ".vcf.gz"), actual_keys, [f"external_{role}"], calls)
                invocations.append(role)

            with patch("master_resources.core.run", side_effect=phase_fixture):
                base_panel(projection, cohort, inventory, 22, root / "base")
            self.assertEqual(invocations, list(core.ROLES))


if __name__ == "__main__":
    unittest.main()
