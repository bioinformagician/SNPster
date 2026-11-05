import copy
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import core
import evaluation
from make_smoke import fixture
import models
import population_followup as population


class PopulationContracts(unittest.TestCase):
    def setUp(self):
        self.samples = [
            {"sample_id": "r1", "role": "reference", "fine": "KHV", "source": "SEA"},
            {"sample_id": "r2", "role": "reference", "fine": "CHS", "source": "EAS"},
            {"sample_id": "r3", "role": "reference", "fine": "CDX", "source": "EAS"},
            {"sample_id": "r4", "role": "reference", "fine": "KHV", "source": "SEA"},
            {"sample_id": "heldout", "role": "validation", "fine": "HELDOUT_ONLY", "source": "SEA"},
        ]
        self.sources = ["SEA", "EAS"]

    def test_reference_only_named_mapping_is_deterministic(self):
        params, rows = population.population_design(self.samples, self.sources)
        again, _ = population.population_design(list(reversed(self.samples)), self.sources)
        self.assertEqual(params, again)
        self.assertEqual(params["components"], ["CDX", "CHS", "KHV"])
        self.assertEqual(params["component_map"], [[0, 1], [0, 1], [1, 0]])
        self.assertEqual([row["reference_people"] for row in rows], [1, 1, 2])
        self.assertTrue(all(row["small_reference_warning"] for row in rows))
        meta = {"reference": self.samples[:4], "sources": self.sources}
        self.assertEqual(models.flare_reference_labels(meta, params, [0, 1, 2, 3]),
                         ["KHV", "CHS", "CDX", "KHV"])
        self.assertEqual(models.flare_reference_labels(meta, {}, [0, 1, 2, 3]), ["SEA", "EAS", "EAS", "SEA"])

    def test_ambiguous_missing_unsafe_and_wrong_mapping_fail(self):
        for fine in ("", "country with spaces"):
            bad = copy.deepcopy(self.samples)
            bad[0]["fine"] = fine
            with self.assertRaises(ValueError):
                population.population_design(bad, self.sources)
        bad = copy.deepcopy(self.samples)
        bad[1]["fine"] = "KHV"
        with self.assertRaisesRegex(ValueError, "multiple"):
            population.population_design(bad, self.sources)
        params, _ = population.population_design(self.samples, self.sources)
        meta = {"reference": self.samples[:4], "sources": self.sources}
        params["component_map"][0] = [1, 0]
        with self.assertRaisesRegex(ValueError, "one-hot"):
            models.flare_reference_labels(meta, params, list(range(4)))

    def test_population_axis_reordering_preserves_named_source_sums(self):
        params, _ = population.population_design(self.samples, self.sources)
        value = np.asarray([[.2, .3, .5]])
        reordered = copy.deepcopy(params)
        reordered["components"].reverse()
        reordered["component_map"].reverse()
        meta = {"reference": self.samples[:4], "sources": self.sources}
        models.flare_reference_labels(meta, reordered, list(range(4)))
        np.testing.assert_array_equal(value @ np.asarray(params["component_map"]),
                                      value[:, ::-1] @ np.asarray(reordered["component_map"]))
        self.assertNotEqual(models.base_policy({"family": "flare", "params": params}),
                            models.base_policy({"family": "flare", "params": reordered}))

    def test_validation_requires_explicit_components_mapping_and_no_model_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = fixture(Path(temporary) / "inputs", {"flare_jar": "/absolute/flare.jar"})
            config["models"] = [{"id": "fine", "family": "flare", "params": {"reference_label": "fine"}}]
            with self.assertRaisesRegex(ValueError, "explicit components/mapping"):
                core.validate_config(config)
            config["models"][0]["params"].update(components=["rB", "rA"], component_map=[[1, 0], [0, 1]])
            core.validate_config(config)
            fine_base = models.base_policy(config["models"][0])
            self.assertIn("reference_label", fine_base["params"])
            self.assertNotEqual(fine_base, models.base_policy({"family": "flare", "params": {}}))
            for value in ("country", "", False):
                bad = copy.deepcopy(config)
                bad["models"][0]["params"]["reference_label"] = value
                with self.assertRaisesRegex(ValueError, "reference_label"):
                    core.validate_config(bad)
            config["models"][0]["params"]["model_file"] = "/absolute/old.model"
            with self.assertRaisesRegex(ValueError, "no supplied model"):
                core.validate_config(config)

    def test_native_mapping_sums_probabilities_and_does_not_fit_from_labels(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            aggregate = root / "aggregate"
            aggregate.mkdir()
            params, _ = population.population_design(self.samples, self.sources)
            values = np.asarray([[.3, .25, .45], [.1, .2, .7], [.2, .2, .6]])
            truth = np.asarray([[1, 0], [0, 1], [.5, .5]])
            np.savez(aggregate / "aggregate.npz", values=values, truth=truth)
            metadata = {
                "id": "fine", "sources": self.sources, "components": params["components"],
                "config": {"seed": 1, "models": [{"id": "fine", "params": params}],
                           "interpretations": [{"id": "native", "kind": "native", "params": {}}]},
                "samples": [{"sample_id": str(i), "role": role}
                            for i, role in enumerate(("calibration", "calibration", "validation"))],
                "sha256": core.digest(aggregate / "aggregate.npz"),
            }
            core.write(aggregate / "aggregate.json", metadata)
            evaluation.interpret(aggregate, "native", root / "mapped")
            saved = np.load(root / "mapped/interpreted.npz")
            np.testing.assert_allclose(saved["values"], values @ np.asarray(params["component_map"]))
            np.testing.assert_allclose(saved["values"][0], [.45, .55])
            np.testing.assert_array_equal(saved["truth"], truth)
            note = core.read(root / "mapped/interpretation.json")
            self.assertEqual(note["calibration_people"], [])
            trained = models.load(root / "mapped/interpreter.pkl")
            self.assertEqual(trained["fit_role"], "reference_metadata_only")
            np.savez(aggregate / "aggregate.npz", values=values, truth=truth[:, ::-1])
            metadata["sha256"] = core.digest(aggregate / "aggregate.npz")
            core.write(aggregate / "aggregate.json", metadata)
            evaluation.interpret(aggregate, "native", root / "changed_truth")
            np.testing.assert_array_equal(np.load(root / "changed_truth/interpreted.npz")["values"], saved["values"])


class Readiness(unittest.TestCase):
    def pilot_spec(self):
        return core.read(Path(__file__).parent / "comparisons/population_chr22_parallel_20261031/specification.json")

    def test_parallel_scope_is_one_chr22_candidate_and_fixed_baseline(self):
        spec = self.pilot_spec()
        population.validate_specification(spec)
        self.assertNotIn("predecessor", spec)
        for field, value in (("chromosomes", [20, 22]), ("maximum_new_singleton_calls", 680),
                             ("baseline_id", "full_em280"), ("baseline_policy", "read_current_test")):
            bad = dict(spec, **{field: value})
            with self.subTest(field=field), self.assertRaises(ValueError):
                population.validate_specification(bad)

    def test_parallel_baseline_never_accesses_running_final_or_test(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "java-profile").write_text("c1\n")
            config = root / "screen/results/frozen/frozen/config.json"
            config.parent.mkdir(parents=True)
            model = {"id": "full_em100", "family": "flare", "params": {
                "em": True, "gen": 10, "states": 100, "update-p": False}}
            core.write(config, {"version": 2, "models": [model]})
            summary = root / "screen/results/summary/summary/summary.json"
            summary.parent.mkdir(parents=True)
            core.write(summary, {"test_evaluated": False})
            spec = dict(self.pilot_spec(), guided_run=str(root))
            with patch.object(population, "completed_predecessor", side_effect=AssertionError("Final must not be read")):
                with patch.object(population, "verify_manifest", return_value={"verified": True}):
                    source, results, _, chosen, _ = population.baseline_source(spec)
                    self.assertEqual(source, root)
                    self.assertEqual(results, root / "screen/results")
                    self.assertEqual(chosen, model)
                    core.write(summary, {"test_evaluated": True})
                    with self.assertRaisesRegex(ValueError, "development-only"):
                        population.baseline_source(spec)
                    core.write(summary, {"test_evaluated": False})
                    model["params"]["states"] = 280
                    core.write(config, {"version": 2, "models": [model]})
                    with self.assertRaisesRegex(ValueError, "predeclared"):
                        population.baseline_source(spec)

    def test_running_predecessor_refused_without_loading_test_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(population, "verify_manifest") as verify:
                with self.assertRaisesRegex(ValueError, "Wait for the current final sweep"):
                    population.completed_predecessor(root)
                verify.assert_not_called()

    def test_completed_predecessor_uses_only_development_scores_and_lock(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "reporting/report").mkdir(parents=True)
            (root / "reporting/report/complete.json").write_text("test contents must not be read")
            (root / "reports").mkdir()
            core.write(root / "reports/scoped-cache-cleanup.json", {"durable_outputs_reverified": True})
            (root / "java-profile").write_text("c1\n")
            (root / "testplanning/design").mkdir(parents=True)
            model = {"id": "winner", "family": "flare", "params": {"em": True, "states": 100}}
            core.write(root / "testplanning/design/winner-lock.json", {
                "development_summary_sha256": "sha", "test_seen_before_lock": False, "model": model})
            result = root / "genome/results"
            (result / "summary/summary").mkdir(parents=True)
            core.write(result / "summary/summary/summary.json", {
                "test_evaluated": False, "combine_eas_sea_for_metrics": True,
                "validation_leaderboard": [{"model": "winner", "validation_mean_person_tv_pp": 1}]})
            (result / "frozen/frozen").mkdir(parents=True)
            core.write(result / "frozen/frozen/config.json", {"models": [model]})
            with patch.object(population, "verify_manifest", return_value={"verified": True, "summary_sha256": "sha"}):
                _, chosen, _ = population.completed_predecessor(root)
            self.assertEqual(chosen, model)
            self.assertFalse((root / "test/results").exists())

    def test_launcher_refuses_unready_predecessor_before_creating_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "bin"
            tools.mkdir()
            docker = tools / "docker"
            docker.write_text("#!/bin/sh\nexit 1\n")
            docker.chmod(0o755)
            spec = root / "spec.json"
            spec.write_text("{}")
            run = root / "run"
            env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"])
            script = Path(__file__).with_name("launch_population.sh")
            result = subprocess.run(["bash", str(script), str(spec), str(run)], env=env,
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(run.exists())

    def test_parallel_launcher_passes_eleven_task_cpu_budget_and_freezes_runtime(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            spec = root / "spec.json"
            core.write(spec, self.pilot_spec())
            tools = root / "bin"
            tools.mkdir()
            docker = tools / "docker"
            calls = root / "calls.jsonl"
            docker.write_text(
                "#!/usr/bin/env python3\n"
                "import json,sys\n"
                f"with open({str(calls)!r},'a') as f: f.write(json.dumps(sys.argv[1:])+'\\n')\n"
                "if sys.argv[1]=='run': raise SystemExit(0 if 'check-ready' in sys.argv else 1)\n"
                "raise SystemExit(0)\n")
            docker.chmod(0o755)
            env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"],
                       SNPSTER_CPUS="11", SNPSTER_MEMORY_GB="80")
            env.pop("SNPSTER_INFERENCE_CHUNK_SIZE", None)
            run = root / "run"
            result = subprocess.run(["bash", str(Path(__file__).with_name("launch_population.sh")), str(spec), str(run)],
                                    env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            import json
            commands = [json.loads(line) for line in calls.read_text().splitlines()]
            controller = next(row for row in commands if row[0] == "run" and "check-ready" not in row)
            self.assertIn("BUDGET_CPUS=11", controller)
            self.assertIn("BUDGET_MEMORY=80", controller)
            self.assertIn("INFERENCE_CHUNK_SIZE=1", controller)
            self.assertIn("-t", controller)
            self.assertIn('--budget_cpus "$BUDGET_CPUS"', controller[-1])
            self.assertIn('--inference_chunk_size "$INFERENCE_CHUNK_SIZE"', controller[-1])
            self.assertIn("-ansi-log true", controller[-1])
            self.assertEqual((run / "java-profile").read_text().strip(), "c1")
            self.assertTrue((run / "runtime.sha256").is_file())
            self.assertEqual(core.read(run / "experiment.json")["maximum_new_singleton_calls"], 340)
            frozen_source = (run / "source.sha256").read_bytes()
            (run / "planning/launch/.nextflow").mkdir(parents=True)
            env["SNPSTER_INFERENCE_CHUNK_SIZE"] = "3"
            resumed = subprocess.run(
                ["bash", str(Path(__file__).with_name("launch_population.sh")), str(spec), str(run), "--resume"],
                env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
            commands = [json.loads(line) for line in calls.read_text().splitlines()]
            controller = [row for row in commands if row[0] == "run" and "check-ready" not in row][-1]
            self.assertIn("INFERENCE_CHUNK_SIZE=3", controller)
            self.assertEqual((run / "source.sha256").read_bytes(), frozen_source)

    def test_invalid_inference_chunk_size_is_rejected_before_docker_or_freezing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            spec = root / "spec.json"
            core.write(spec, self.pilot_spec())
            tools = root / "bin"
            tools.mkdir()
            docker = tools / "docker"
            docker.write_text(f"#!/bin/sh\nprintf 'called\\n' >> '{root}/calls'\nexit 1\n")
            docker.chmod(0o755)
            for value in ("0", "-1", "1.5", "invalid"):
                with self.subTest(value=value):
                    env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"],
                               SNPSTER_INFERENCE_CHUNK_SIZE=value)
                    result = subprocess.run(
                        ["bash", str(Path(__file__).with_name("launch_population.sh")), str(spec), str(root / "run")],
                        env=env, capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                    self.assertIn("SNPSTER_INFERENCE_CHUNK_SIZE must be a positive integer", result.stderr)
                    self.assertFalse((root / "calls").exists())
                    self.assertFalse((root / "run").exists())

    def test_resume_rejects_changed_frozen_runtime_before_controller_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run = root / "run"
            (run / "planning/launch/.nextflow").mkdir(parents=True)
            (run / "runtime").mkdir()
            (run / "source").mkdir()
            (run / "population-followup").write_text("true\n")
            (run / "java-profile").write_text("c1\n")
            overlay = run / "runtime/java-c1-pinned.config"
            overlay.write_text("pinned runtime\n")
            (run / "runtime.sha256").write_text(f"{core.digest(overlay)}  runtime/java-c1-pinned.config\n")
            (run / "experiment.json").write_text("{}")
            tools = root / "bin"
            tools.mkdir()
            docker = tools / "docker"
            docker.write_text(f"#!/bin/sh\nprintf 'called\\n' >> '{root}/calls'\nexit 1\n")
            docker.chmod(0o755)
            env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"])
            command = ["bash", str(Path(__file__).with_name("launch_population.sh")),
                       str(run / "experiment.json"), str(run), "--resume"]
            original = subprocess.run(command, env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(original.returncode, 1)
            self.assertEqual((root / "calls").read_text(), "called\n")
            overlay.write_text("different runtime\n")
            changed = subprocess.run(command, env=env, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("FAILED", changed.stdout)
            self.assertEqual((root / "calls").read_text(), "called\n")

    def test_cleanup_rejects_changed_report_and_preserves_durable_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "population-followup").write_text("true\n")
            (root / "source.sha256").write_text("frozen source\n")
            report = root / "reporting/report/population-comparison.json"
            report.parent.mkdir(parents=True)
            core.write(report, {"test_evaluated": False})
            complete = report.with_name("complete.json")
            core.write(complete, {"summary_sha256": "sha", "report_sha256": core.digest(report)})
            caches = [root / stage / "work" for stage in ("planning", "screen", "reporting")]
            for path in caches:
                path.mkdir(parents=True)
            with patch.object(population, "verify_manifest", return_value={"summary_sha256": "wrong"}):
                with self.assertRaisesRegex(ValueError, "changed"):
                    population.cleanup(root, root / "cleanup.json")
            self.assertTrue(all(path.is_dir() for path in caches))
            with patch.object(population, "verify_manifest", return_value={"summary_sha256": "sha"}):
                population.cleanup(root, root / "cleanup.json")
            self.assertFalse(any(path.exists() for path in caches))
            self.assertTrue(report.is_file() and complete.is_file())
            self.assertTrue(core.read(root / "cleanup.json")["durable_outputs_reverified"])

    def test_cleanup_refuses_redirected_work_cache(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "population-followup").write_text("true\n")
            (root / "source.sha256").write_text("frozen source\n")
            report = root / "reporting/report/population-comparison.json"
            report.parent.mkdir(parents=True)
            core.write(report, {})
            core.write(report.with_name("complete.json"),
                       {"summary_sha256": "sha", "report_sha256": core.digest(report)})
            other = root / "unrelated"
            other.mkdir()
            (root / "planning").mkdir()
            (root / "planning/work").symlink_to(other, target_is_directory=True)
            with patch.object(population, "verify_manifest", return_value={"summary_sha256": "sha"}):
                with self.assertRaisesRegex(ValueError, "redirected"):
                    population.cleanup(root, root / "cleanup.json")
            self.assertTrue(other.is_dir())


if __name__ == "__main__":
    unittest.main()
