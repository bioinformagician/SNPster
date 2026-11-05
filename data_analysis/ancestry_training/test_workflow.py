import copy
import gzip
import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np

import core
import evaluation
import models
from prepare_observed import components
from make_smoke import fixture
from prepare_comparison import derive


class Contracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch = tempfile.TemporaryDirectory()
        cls.root = Path(cls.scratch.name)
        cls.config = fixture(cls.root / "inputs", {}, models=[{"id": "freq", "family": "nnls", "params": {}}])

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def test_valid_config_and_roles(self):
        core.validate_config(self.config)
        core.validate_roles(core.rows(self.config["samples"]), self.config["sources"])

    def test_unknown_hyperparameters_and_cluster_misuse(self):
        for params in ({"k": 3}, {"foo": 10}):
            config = copy.deepcopy(self.config)
            config["models"][0]["params"] = params
            with self.assertRaisesRegex(ValueError, "hyperparameters"):
                core.validate_config(config)

    def test_invalid_hyperparameters(self):
        for params in ({"n_estimators": 0}, {"crf_c2": -1}, {"em": "false"}, {"context_ratio": 2},
                       {"alpha": float("nan")}, {"gen": True}, {"k": 1.5}):
            with self.subTest(params=params), self.assertRaises(ValueError):
                core.validate_parameters(params)

    def test_full_scope_required(self):
        config = copy.deepcopy(self.config)
        config["mode"] = "full"
        with self.assertRaisesRegex(ValueError, "Partial-autosome"):
            core.validate_config(config)

    def test_smoke_not_certification(self):
        config = copy.deepcopy(self.config)
        config["independent_test"] = True
        with self.assertRaisesRegex(ValueError, "smoke"):
            core.validate_config(config)

    def test_comparison_smoke_preserves_models_and_split(self):
        config = copy.deepcopy(self.config)
        config.update(mode="full", chromosomes=list(range(1, 23)))
        for key in ("reference", "training", "maps", "spans"):
            config[key] = {str(chrom): config[key]["22"] for chrom in config["chromosomes"]}
        for route in config["routes"]:
            route["vcfs"] = {str(chrom): route["vcfs"]["22"] for chrom in config["chromosomes"]}
        original = copy.deepcopy(config)
        smoke = derive(config, "smoke")
        self.assertEqual(config, original)
        self.assertEqual(smoke["chromosomes"], [22])
        for key in ("reference", "training", "maps", "spans"):
            self.assertEqual(set(smoke[key]), {"22"})
        for route in smoke["routes"]:
            self.assertEqual(set(route["vcfs"]), {"22"})
        for key in ("models", "interpretations", "samples", "segments", "seed"):
            self.assertEqual(smoke[key], config[key])
        self.assertFalse(smoke["independent_test"])

    def test_focused_comparison_grid_and_base_reuse(self):
        config = copy.deepcopy(self.config)
        config.update(mode="full", chromosomes=list(range(1, 23)), allow_research_only=True)
        for key in ("reference", "training", "maps", "spans"):
            config[key] = {str(chrom): config[key]["22"] for chrom in config["chromosomes"]}
        for route in config["routes"]:
            route["vcfs"] = {str(chrom): route["vcfs"]["22"] for chrom in config["chromosomes"]}
        config["tools"] = {"flare_jar": "/absolute/flare.jar", "gnomix_root": "/absolute/gnomix"}
        config["models"] = [
            {"id": f"flare_{reference}_{mode}_g010", "family": "flare",
             "params": {"gen": 10, "em": mode == "em", **({"reference_cap": 16} if reference == "cap16" else {})}}
            for reference in ("full", "cap16") for mode in ("fixed", "em")]
        config["models"] += [
            {"id": f"gnomix_fast_w050_c2_{label}", "family": "gnomix",
             "params": {"window_size_cM": .5, "crf_c2": value}}
            for label, value in (("0p1", .1), ("1", 1))]
        config["models"] += [
            {"id": "frequency_control", "family": "nnls", "params": {}},
            {"id": "pca_ridge_control", "family": "pca_ridge", "params": {}}]
        config["interpretations"] = [{"id": kind, "kind": kind, "params": {}}
                                    for kind in ("native", "simplex", "ridge", "rf", "xgboost")]
        focused = derive(config, "focused16")
        self.assertEqual(len(focused["models"]), 16)
        self.assertEqual(len(models.plan(focused)["bases"]), 6)
        fixed = [model for model in focused["models"] if model["family"] == "flare" and not model["params"]["em"]]
        self.assertEqual(len(fixed), 8)
        self.assertEqual({model["params"]["gen"] for model in fixed}, {2, 5, 10, 20})
        self.assertEqual({model["params"].get("reference_cap") for model in fixed}, {None, 16})
        self.assertEqual(focused["samples"], config["samples"])
        self.assertEqual(focused["interpretations"], config["interpretations"])
        with self.assertRaisesRegex(ValueError, "Missing reviewed"):
            derive({**config, "models": config["models"][:-1]}, "focused16")

    def test_roles_and_donors_purged(self):
        samples = core.rows(self.config["samples"])
        for field in ("family_id", "donors"):
            altered = copy.deepcopy(samples)
            altered[-1][field] = altered[0][field]
            with self.assertRaisesRegex(ValueError, "leakage"):
                core.validate_roles(altered, self.config["sources"])
        altered = copy.deepcopy(samples)
        altered[-1]["donors"] = altered[0]["sample_id"]
        with self.assertRaisesRegex(ValueError, "leakage"):
            core.validate_roles(altered, self.config["sources"])

    def test_founder_components_are_transitive_and_deterministic(self):
        groups = [["a", "b"], ["c", "d"], ["b", "c"], ["e", "f"]]
        expected = [["a", "b", "c", "d"], ["e", "f"]]
        self.assertEqual(components(groups), expected)
        self.assertEqual(components(groups[::-1]), expected)

    def test_missing_and_duplicate_people_rejected(self):
        samples = core.rows(self.config["samples"])
        for altered in (samples + samples[:1], [row for row in samples if row["role"] != "calibration"]):
            with self.assertRaises(ValueError):
                core.validate_roles(altered, self.config["sources"])

    def test_physical_endpoints(self):
        np.testing.assert_array_equal(core.physical_weights([10, 20, 90], [0, 100]), [15, 40, 45])
        for pos in ([10, 10], [20, 10], [10, 100]):
            with self.assertRaises(ValueError):
                core.physical_weights(pos, [0, 100])

    def test_invalid_native_probabilities_not_repaired(self):
        for value in ([.4, .4], [float("nan"), 1], [-.1, 1.1], [.5, .5, 0]):
            with self.assertRaises(ValueError):
                core.probability(value, 2)
        np.testing.assert_allclose(core.project_simplex([[1.3, -.2]]), [[1, 0]])

    def test_reference_base_reused_across_crf_l2(self):
        config = copy.deepcopy(self.config)
        config["models"] = [
            {"id": "a", "family": "gnomix", "params": {"window_size_cM": .5, "crf_c2": .1}},
            {"id": "b", "family": "gnomix", "params": {"window_size_cM": .5, "crf_c2": 1}},
        ]
        self.assertEqual(len(models.plan(config)["bases"]), 1)

    def test_flare_array_alias_not_double_counted(self):
        config = copy.deepcopy(self.config)
        config["tools"] = {"flare_jar": "/absolute/flare.jar"}
        config["models"] = [
            {"id": "a", "family": "flare", "params": {"array": True, "em": False}},
            {"id": "b", "family": "flare", "params": {"array": False, "em": False}},
        ]
        with self.assertRaisesRegex(ValueError, "alias"):
            core.validate_config(config)

    def test_vcf_phase_and_identity_order(self):
        path = self.root / "inputs/query.vcf.gz"
        ids = ["test_2", "calibration_0"]
        keys, calls = core.vcf(path, 22, ids, True)
        _, reversed_calls = core.vcf(path, 22, ids[::-1], True)
        np.testing.assert_array_equal(calls[::-1], reversed_calls)
        self.assertEqual(len(keys), 256)
        self.assertEqual(calls.dtype, np.dtype(np.int8))
        with self.assertRaisesRegex(ValueError, "Unphased"):
            core.vcf(self.root / "inputs/harmonized.vcf.gz", 22, ids, True)
        with self.assertRaisesRegex(ValueError, "identity"):
            core.vcf(path, 22, ["unknown"], True)

    def test_truth_realized_and_scope(self):
        samples = core.rows(self.config["samples"])
        truth, local = core.truth(self.config, samples, 22, np.array([501, 1501]))
        self.assertEqual(local.shape, (20, 2, 2))
        self.assertTrue(any(np.array_equal(row, [.5, .5]) for row in truth))
        config = copy.deepcopy(self.config)
        config["spans"]["22"] = [0, 256001]
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            core.truth(config, samples, 22)

    def test_missing_genotype_not_filled(self):
        path = self.root / "missing.vcf"
        path.write_text(
            '##fileformat=VCFv4.2\n##contig=<ID=22>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n'
            '22\t100\t.\tA\tG\t.\t.\t.\tGT\t./.\n')
        with self.assertRaisesRegex(ValueError, "filling/DS substitution"):
            core.vcf(path, 22, ["s"], False)

    def test_phase_adapter_rejects_changed_gt_or_key(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panel = root / "panel"
            panel.mkdir()
            keys = [("22", 10, "A", "G"), ("22", 90, "A", "G")]
            query = np.array([[[0, 0], [0, 1]]], dtype=np.int8)
            np.savez(panel / "data.npz", query=query)
            config = copy.deepcopy(self.config)
            config["models"] = [{"family": "flare"}]
            config["tools"] = {"beagle_jar": "fake"}
            core.write(panel / "panel.json", {
                "queries": [{"sample_id": "s"}], "keys": keys, "config": config,
                "chromosome": 22, "route": {"stage": "harmonized"}, "tool_sha256": {},
            })
            for index, changed in enumerate(((keys, np.ones_like(query)),
                                             ([("22", 10, "A", "T"), keys[1]], query))):
                with patch("core.run"), patch("core.vcf", return_value=changed):
                    with self.assertRaisesRegex(ValueError, "changed observed keys or diploid GT"):
                        core.phase(panel, root / str(index))

    def test_exact_allele_panel_intersection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = copy.deepcopy(self.config)
            ids = [row["sample_id"] for row in core.rows(config["samples"])
                   if row["role"] in ("calibration", "validation", "test")]
            keys, calls = core.vcf(config["routes"][1]["vcfs"]["22"], 22, ids, True)
            keys[0] = (keys[0][0], keys[0][1], "A", "T")
            target = root / "allele_mismatch.vcf.gz"
            core.write_vcf(target, keys, ids, calls)
            config["routes"] = [{"id": "mismatch", "stage": "imputed", "vcfs": {"22": str(target)}}]
            config["reference_filters"]["minimum_markers"] = 256
            core.write(root / "config.json", config)
            core.freeze(root / "config.json", root / "frozen")
            with self.assertRaisesRegex(ValueError, "Too few exact common.*255"):
                core.prepare(root / "frozen", "mismatch", 22, root / "panel")

    def test_flare_rounding_and_exact_keys(self):
        path = self.root / "posterior.vcf.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("22\t10\t.\tA\tG\t.\t.\t.\tANP1:ANP2\t0.333:0.333\n")
        with self.assertRaises(ValueError):
            models.flare_posterior(path, [("22", 10, "A", "G"), ("22", 20, "A", "G")], ["A", "B"], [0, 100])
        with gzip.open(path, "wt") as handle:
            for pos in (10, 90):
                handle.write(f"22\t{pos}\t.\tA\tG\t.\t.\t.\tANP1:ANP2\t0.333,0.666:0.333,0.666\n")
        value, correction = models.flare_posterior(path, [("22", 10, "A", "G"), ("22", 90, "A", "G")],
                                                   ["A", "B"], [0, 100])
        np.testing.assert_allclose(value, [1 / 3, 2 / 3])
        self.assertGreater(correction, 0)

    def test_native_vs_training_only_component_interpretation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            aggregate = root / "aggregate"
            aggregate.mkdir()
            config = copy.deepcopy(self.config)
            samples = [row for row in core.rows(config["samples"]) if row["role"] in ("calibration", "validation", "test")]
            true = np.array([[1, 0], [0, 1], [.5, .5], [.5, .5]] * 3)
            # Raw component axis is the reverse of named sources.
            np.savez(aggregate / "aggregate.npz", values=true[:, ::-1], truth=true)
            core.write(aggregate / "aggregate.json", {
                "id": "freq", "sources": ["B", "A"], "components": ["anc0", "anc1"],
                "config": config, "samples": samples,
                "sha256": core.digest(aggregate / "aggregate.npz"),
            })
            evaluation.interpret(aggregate, "native", root / "native")
            self.assertIsNotNone(core.read(root / "native/interpretation.json")["unscored_reason"])
            evaluation.interpret(aggregate, "simplex", root / "mapped")
            trained = np.load(root / "mapped/interpreted.npz")["values"]
            np.testing.assert_allclose(trained, true, atol=.002)
            # Alter every evaluation label, preserving calibration labels: learned mapping must not change.
            true[4:] = true[4:, ::-1]
            np.savez(aggregate / "aggregate.npz", values=np.load(root / "mapped/interpreted.npz")["truth"][:, ::-1], truth=true)
            meta = core.read(aggregate / "aggregate.json")
            meta["sha256"] = core.digest(aggregate / "aggregate.npz")
            core.write(aggregate / "aggregate.json", meta)
            evaluation.interpret(aggregate, "simplex", root / "mapped_again")
            np.testing.assert_allclose(np.load(root / "mapped_again/interpreted.npz")["values"], trained, atol=1e-12)

    def test_seven_source_flare_export_rounding_is_component_bounded(self):
        raw = [0, .98, 0, 0, 0, 0, .01]
        np.testing.assert_array_equal(core.rounded_flare_probability(raw, 7), raw)
        for invalid in ([.6, .42, 0, 0, 0, 0, 0], [1.01, 0, 0, 0, 0, 0, 0],
                        [0] * 7, [-.001, 1, 0, 0, 0, 0, 0]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                core.rounded_flare_probability(invalid, 7)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "native.vcf.gz"
            with gzip.open(path, "wt") as handle:
                for pos in (10, 90):
                    handle.write(f"22\t{pos}\t.\tA\tG\t.\t.\t.\tANP1:ANP2\t0,0.98,0,0,0,0,0.01:0,0.98,0,0,0,0,0.01\n")
            value, correction = models.flare_posterior(
                path, [("22", 10, "A", "G"), ("22", 90, "A", "G")],
                [f"source{i}" for i in range(7)], [0, 100])
            np.testing.assert_allclose(value, np.asarray(raw) / .99)
            self.assertAlmostEqual(correction, .01)

    def test_aggregation_rejects_missing_and_duplicate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prediction = root / "prediction"
            prediction.mkdir()
            core.write(prediction / "prediction.json", {"config": self.config, "chromosome": 1})
            for index, folders in enumerate(([prediction], [prediction, prediction])):
                with self.assertRaisesRegex(ValueError, "Incomplete or duplicate"):
                    evaluation.aggregate(folders, root / str(index))

    def test_physical_chromosome_aggregation_not_mean_tv(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = copy.deepcopy(self.config)
            config["chromosomes"] = [1, 22]
            config["spans"] = {"1": [0, 10], "22": [0, 90]}
            folders = []
            for chrom, values, true in ((1, [[1, 0]], [[0, 1]]), (22, [[0, 1]], [[1, 0]])):
                folder = root / str(chrom)
                folder.mkdir()
                np.savez(folder / "predictions.npz", values=values, truth=true)
                core.write(folder / "prediction.json", {
                    "id": "freq", "family": "nnls", "license_track": "permissive_software",
                    "route": config["routes"][0], "components": config["sources"], "sources": config["sources"],
                    "config": config, "chromosome": chrom, "span": config["spans"][str(chrom)],
                    "samples": [{"sample_id": "s"}], "sha256": core.digest(folder / "predictions.npz"),
                })
                folders.append(folder)
            evaluation.aggregate(folders, root / "aggregate")
            data = np.load(root / "aggregate/aggregate.npz")
            np.testing.assert_allclose(data["values"], [[.1, .9]])
            self.assertAlmostEqual(float(.5 * np.abs(data["values"] - data["truth"]).sum()), .8)

    def test_finalists_selected_without_test_labels(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = copy.deepcopy(self.config)
            config["models"] = [{"id": name, "family": "nnls", "params": {}} for name in ("good", "bad")]
            config["interpretations"] = [{"id": "native", "kind": "native", "params": {}}]
            samples = [row for row in core.rows(config["samples"])
                       if row["role"] in ("calibration", "validation", "test")]
            true = np.array([[1, 0], [0, 1], [1, 0], [0, 1]] * 3, dtype=float)
            folders = []
            for route in config["routes"]:
                for model in config["models"]:
                    folder = root / f"{route['id']}_{model['id']}"
                    folder.mkdir()
                    values = true.copy()
                    values[4:8] = true[4:8] if model["id"] == "good" else true[4:8, ::-1]
                    values[8:] = true[8:, ::-1] if model["id"] == "good" else true[8:]
                    np.savez(folder / "interpreted.npz", values=values, truth=true)
                    core.write(folder / "interpretation.json", {
                        "id": model["id"], "route": route, "interpreter": config["interpretations"][0],
                        "sources": config["sources"], "components": config["sources"], "samples": samples,
                        "license_track": "permissive_software", "unscored_reason": None, "config": config,
                        "sha256": core.digest(folder / "interpreted.npz"),
                    })
                    folders.append(folder)
            evaluation.select(folders, root / "selected")
            summary = core.read(root / "selected/summary.json")
            self.assertEqual({row["model"] for row in summary["validation_selected_finalists"]}, {"good"})
            self.assertTrue(all(row["selected_test_mean_person_tv_pp"] == 100
                                for row in summary["validation_selected_finalists"]))
            self.assertTrue(evaluation.verify(root / "selected")["verified"])

    def test_published_hash_failure(self):
        from audit import verify_manifest
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "artifact").write_text("original")
            core.write(root / "independent_verification.json", {"files": {"artifact": core.digest(root / "artifact")}})
            (root / "artifact").write_text("changed")
            with self.assertRaisesRegex(ValueError, "Durable artifact changed"):
                verify_manifest(root)


if __name__ == "__main__":
    unittest.main()
