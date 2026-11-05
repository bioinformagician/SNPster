import copy
import csv
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

import core
import evaluation
from make_smoke import fixture
import models
import projection
import reuse
import staged_comparison as staged


def cohort():
    sources = [f"S{i}" for i in range(7)]
    samples = []
    for role, families in (("reference", 7), ("train", 7), ("calibration", 7), ("validation", 14), ("test", 14)):
        for family in range(families):
            tracks = ("pure",) if role in ("reference", "train") else (
                "pure", "pure", "pure", "pure", "mixed_g1", "mixed_g1", "mixed_g2", "mixed_g3")
            for index, track in enumerate(tracks):
                identity = f"{role}_{family}_{index}"
                samples.append({
                    "sample_id": identity, "role": role, "family_id": f"{role}_family_{family}",
                    "donors": f"{role}_donor_{family}", "source": sources[family % 7] if track == "pure" else "",
                    "track": track, "fine": f"fine_{family}_{index % 2}" if track == "pure" else "",
                })
    return samples, sources


class DevelopmentScope(unittest.TestCase):
    def test_deterministic_source_family_generation_scope_without_test(self):
        samples, sources = cohort()
        original = copy.deepcopy(samples)
        selected, report = staged.development_samples(samples, sources, 19)
        reordered, _ = staged.development_samples(list(reversed(samples)), sources, 19)
        self.assertEqual(samples, original)
        self.assertEqual({row["sample_id"] for row in selected}, {row["sample_id"] for row in reordered})
        self.assertEqual(report["calibration"]["people"], 32)
        self.assertEqual(report["validation"]["people"], 64)
        self.assertEqual(report["calibration"]["families"], 7)
        self.assertEqual(report["validation"]["families"], 14)
        for role in ("reference", "train"):
            self.assertEqual([row for row in selected if row["role"] == role],
                             [row for row in samples if row["role"] == role])
        self.assertFalse(any(row["role"] == "test" for row in selected))
        for role in ("calibration", "validation"):
            self.assertEqual(set(report[role]["pure_sources"]), set(sources))
            self.assertEqual(set(report[role]["by_track"]), {"pure", "mixed_g1", "mixed_g2", "mixed_g3"})
        with self.assertRaises(ValueError):
            core.validate_roles(selected, sources)
        core.validate_roles(selected, sources, development_only=True)
        with self.assertRaises(ValueError):
            core.validate_roles(samples, sources, development_only=True)

    def test_missing_generation_is_an_explicit_failure(self):
        samples, sources = cohort()
        samples = [row for row in samples if row["track"] != "mixed_g3"]
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            staged.development_samples(samples, sources, 19)

    def test_measured_call_reduction(self):
        selected, _ = staged.development_samples(*cohort(), 19)
        config = {"routes": [{}, {}], "chromosomes": list(range(1, 23)),
                  "models": [{"family": "flare"}] * 5}
        self.assertEqual(staged.call_count(config, selected), 21120)
        config["models"] *= 2
        old = [{"role": "calibration"}] * 857
        self.assertEqual(staged.call_count(config, old), 377080)

    def test_exactly_five_distinct_validation_selected_models_and_reserves(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = fixture(Path(temporary) / "inputs", {"flare_jar": "/absolute/pinned.jar"})
            config.update(version=2, mode="development_full", chromosomes=list(range(1, 23)))
            for key in ("reference", "training", "maps", "spans"):
                config[key] = {str(chrom): config[key]["22"] for chrom in config["chromosomes"]}
            for route in config["routes"]:
                route["vcfs"] = {str(chrom): route["vcfs"]["22"] for chrom in config["chromosomes"]}
            config["models"] = [{"id": f"f{i}", "family": "flare", "params": {"gen": i + 1}}
                                for i in range(10)] + [{"id": "control", "family": "nnls", "params": {}}]
            summary = {"scope": "development_screen", "test_evaluated": False,
                       "validation_leaderboard": [
                           {"model": f"f{i}", "interpreter": head, "stage": "harmonized",
                            "validation_mean_person_tv_pp": 10 - i + penalty}
                           for i in range(10) for head, penalty in (("native", 1), ("ridge", 0))]}
            config["interpretations"] = [{"id": name, "kind": name, "params": {}} for name in ("native", "ridge")]
            result, ranked = staged.finalists(config, summary)
            self.assertEqual({model["id"] for model in result["models"]},
                             {"f9", "f8", "f7", "f6", "f5", "control"})
            self.assertEqual(len(ranked[5:]), 5)
            self.assertEqual(config["models"][0]["id"], "f0")
            with self.assertRaisesRegex(ValueError, "exactly five"):
                staged.finalists(config, summary, 6)
            with self.assertRaisesRegex(ValueError, "Only an audited"):
                staged.finalists(config, {**summary, "test_evaluated": True})
            with self.assertRaisesRegex(ValueError, "interpretation coverage"):
                staged.finalists(config, {**summary, "validation_leaderboard": summary["validation_leaderboard"][1:]})
            config["chromosomes"] = [22]
            with self.assertRaisesRegex(ValueError, "Partial-autosome"):
                core.validate_config(config)
            config["mode"] = "development_screen"
            config["independent_test"] = True
            with self.assertRaisesRegex(ValueError, "Development"):
                core.validate_config(config)


class DevelopmentSelection(unittest.TestCase):
    def test_no_test_access_or_test_results_and_independent_metric_reconstruction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = fixture(root / "inputs", {}, models=[{"id": "freq", "family": "nnls", "params": {}}])
            config.update(version=2, mode="development_screen")
            config["interpretations"] = [{"id": "native", "kind": "native", "params": {}}]
            samples = [row for row in core.rows(config["samples"]) if row["role"] in ("calibration", "validation")]
            true = np.tile(np.eye(2), (4, 1))
            folders = []
            for route in config["routes"]:
                folder = root / route["id"]
                folder.mkdir()
                np.savez(folder / "interpreted.npz", values=true, truth=true)
                core.write(folder / "interpretation.json", {
                    "id": "freq", "route": route, "interpreter": config["interpretations"][0],
                    "sources": config["sources"], "components": config["sources"], "samples": samples,
                    "license_track": "permissive_software", "unscored_reason": None, "config": config,
                    "sha256": core.digest(folder / "interpreted.npz"),
                })
                folders.append(folder)
            evaluation.select(folders, root / "selected")
            summary = core.read(root / "selected/summary.json")
            self.assertFalse(summary["test_evaluated"])
            self.assertTrue(all("selected_test_metrics" not in row for row in summary["validation_selected_finalists"]))
            self.assertTrue(evaluation.verify(root / "selected")["verified"])
            changed = core.read(folders[0] / "interpretation.json")
            changed["samples"][0]["role"] = "test"
            core.write(folders[0] / "interpretation.json", changed)
            with self.assertRaisesRegex(ValueError, "only calibration and validation"):
                evaluation.select(folders, root / "bad")


class ProjectedPanels(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        tools = {"flare_jar": os.environ["SNPSTER_TEST_FLARE_JAR"]} if "SNPSTER_TEST_FLARE_JAR" in os.environ else {}
        self.config = fixture(self.root / "inputs", tools, models=[{"id": "freq", "family": "nnls", "params": {}}])
        self.config["routes"] = [self.config["routes"][1]]
        core.write(self.root / "source.json", self.config)
        core.freeze(self.root / "source.json", self.root / "source-frozen")
        self.parent = self.root / "parent"
        core.prepare(self.root / "source-frozen", "imputed", 22, self.parent)
        core.phase(self.parent, self.root / "parent-phase")
        samples = [row for row in core.rows(self.config["samples"])
                   if row["role"] in ("reference", "train")
                   or (row["role"] in ("calibration", "validation") and int(row["sample_id"].split("_")[-1]) < 2)]
        with (self.root / "subset.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(samples[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(samples)
        self.target = {**self.config, "version": 2, "mode": "development_screen",
                       "samples": str(self.root / "subset.tsv")}
        core.write(self.root / "target.json", self.target)
        core.freeze(self.root / "target.json", self.root / "target-frozen")
        self.record = {
            "version": 2, "config": self.target, "projections": {"imputed/22": {
                "panel": {"directory": str(self.parent), "files": reuse.inventory(self.parent)},
                "phase": {"directory": str(self.root / "parent-phase"), "files": reuse.inventory(self.root / "parent-phase")},
                "input_sha256": core.read(self.root / "source-frozen/manifest.json")["input_sha256"],
            }},
        }

    def project(self):
        projection.prepare(self.record, self.root / "target-frozen", "imputed", 22, self.root / "projected")
        projection.phase(self.record, self.root / "projected", self.root / "projected-phase")

    def test_subset_equals_fresh_preprocessing_and_preserves_all_training_axes(self):
        self.project()
        core.prepare(self.root / "target-frozen", "imputed", 22, self.root / "fresh")
        core.phase(self.root / "fresh", self.root / "fresh-phase")
        with np.load(self.root / "projected/data.npz") as actual, np.load(self.root / "fresh/data.npz") as expected:
            self.assertEqual(actual.files, expected.files)
            for name in actual.files:
                np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)
        np.testing.assert_array_equal(np.load(self.root / "projected-phase/query.npy"),
                                      np.load(self.root / "fresh-phase/query.npy"))
        self.assertEqual(core.read(self.root / "projected/panel.json")["config"], self.target)

    def test_source_hash_and_input_changes_are_rejected(self):
        self.record["projections"]["imputed/22"]["input_sha256"][self.config["segments"]] = "changed"
        with self.assertRaisesRegex(ValueError, "different scientific input"):
            self.project()
        self.record["projections"]["imputed/22"]["input_sha256"][self.config["segments"]] = core.digest(self.config["segments"])
        (self.parent / "reference.panel").write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "source artifact changed"):
            self.project()

    def test_phase_policy_and_reference_parameter_changes_are_rejected(self):
        self.project()
        params = {"gen": 10, "em": False}
        before = projection.signature(self.parent, params)
        self.assertEqual(before, projection.signature(self.root / "projected", params))
        self.assertNotEqual(before, projection.signature(self.root / "projected", {**params, "gen": 20}))
        self.assertNotEqual(before, projection.signature(self.root / "projected", {**params, "reference_cap": 2}))
        meta = core.read(self.root / "projected/panel.json")
        meta["route"]["stage"] = "harmonized"
        meta["config"]["models"] = [{"family": "flare"}]
        core.write(self.root / "projected/panel.json", meta)
        with self.assertRaisesRegex(ValueError, "phase-only semantics"):
            projection.phase(self.record, self.root / "projected", self.root / "bad-phase")

    @unittest.skipUnless(os.environ.get("SNPSTER_TEST_FLARE_JAR"), "requires the pinned FLARE jar")
    def test_real_fixed_and_em_native_subset_reuse_equals_fresh(self):
        self.project()
        for mode in (False, True):
            identity = "em" if mode else "fixed"
            params = {"em": mode, "gen": 10, "states": 20}
            for panel in (self.parent, self.root / "projected"):
                meta = core.read(panel / "panel.json")
                meta["config"]["models"] = [{"id": identity, "family": "flare", "params": params}]
                meta["base_plan"] = models.plan(meta["config"])
                core.write(panel / "panel.json", meta)
                base_id = next(iter(meta["base_plan"]["bases"]))
                suffix = "parent" if panel == self.parent else "subset"
                models.fit_base(panel, base_id, self.root / f"{identity}-{suffix}-base")
                models.fit_model(panel, self.root / f"{identity}-{suffix}-base", identity,
                                 self.root / f"{identity}-{suffix}-fit")
            source = self.root / f"{identity}-source"
            models.infer(self.parent, self.root / "parent-phase", self.root / f"{identity}-parent-base",
                         self.root / f"{identity}-parent-fit", source)
            meta = core.read(self.root / "projected/panel.json")
            key = f"imputed/22/{identity}"
            self.record.update(config=meta["config"], native={key: {
                "signature": projection.signature(self.parent, params),
                "samples": {row["sample_id"]: projection.capture(source, meta, params, row["sample_id"],
                    str(self.root / f"{identity}-parent-base"), str(source)) for row in meta["queries"]},
            }})
            reuse_path = self.root / f"{identity}-reuse.json"
            core.write(reuse_path, self.record)
            for label, manifest in (("fresh", None), ("reused", reuse_path)):
                models.infer(self.root / "projected", self.root / "projected-phase",
                             self.root / f"{identity}-subset-base", self.root / f"{identity}-subset-fit",
                             self.root / f"{identity}-{label}", reuse_manifest=manifest)
            np.testing.assert_array_equal(np.load(self.root / f"{identity}-fresh/predictions.npz")["values"],
                                          np.load(self.root / f"{identity}-reused/predictions.npz")["values"])
            self.assertEqual(len(core.read(self.root / f"{identity}-reused/prediction.json")["reused_singletons"]), 4)


class SuccessfulCacheCleanup(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.parent = Path(self.scratch.name)
        self.root = self.parent / "dedicated"
        self.root.mkdir()
        (self.root / "staged-development").write_text("true\n")
        (self.root / "source.sha256").write_text("retained source hashes\n")
        (self.root / "experiment.json").write_text("{}\n")
        for stage in ("planning", "screen", "advancement", "genome"):
            work = self.root / stage / "work"
            work.mkdir(parents=True)
            (work / "disposable").write_text("cache\n")
        for stage in ("screen", "genome"):
            results = self.root / stage / "results"
            results.mkdir()
            (results / "science").write_text("durable\n")
            core.write(results / "independent_verification.json",
                       {"verified": True, "files": {"science": core.digest(results / "science")}})
        self.outside = self.parent / "superseded"
        self.outside.mkdir()
        (self.outside / "failure-evidence").write_text("must retain\n")

    def test_only_four_verified_successful_work_directories_are_removed(self):
        staged.cleanup(self.root, self.root / "cleanup.json")
        self.assertEqual(len(core.read(self.root / "cleanup.json")["removed_exact_successful_caches"]), 4)
        for stage in ("screen", "genome"):
            self.assertTrue((self.root / stage / "results/science").is_file())
        self.assertTrue((self.outside / "failure-evidence").is_file())

    def test_failed_audit_preserves_every_cache(self):
        core.write(self.root / "genome/results/independent_verification.json", {"verified": False, "files": {}})
        with self.assertRaisesRegex(ValueError, "Both scientific stages"):
            staged.cleanup(self.root, self.root / "cleanup.json")
        self.assertTrue((self.root / "planning/work/disposable").is_file())

    def test_redirected_work_directory_is_never_followed(self):
        work = self.root / "planning/work"
        (work / "disposable").unlink()
        work.rmdir()
        work.symlink_to(self.outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "redirected"):
            staged.cleanup(self.root, self.root / "cleanup.json")
        self.assertTrue((self.outside / "failure-evidence").is_file())


if __name__ == "__main__":
    unittest.main()
