import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class LaunchRuntime(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix="snpster-launch-test-")
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.source = Path(__file__).resolve().parent
        self.config = self.root / "config.json"
        self.config.write_text("{}\n")
        self.run = self.root / "run"
        self.calls = self.root / "docker-calls.jsonl"
        tools = self.root / "bin"
        tools.mkdir()
        docker = tools / "docker"
        docker.write_text(
            "#!/usr/bin/env python3\n"
            "import json,os,sys\n"
            "with open(os.environ['DOCKER_TEST_CALLS'],'a') as handle:\n"
            "    handle.write(json.dumps(sys.argv[1:])+'\\n')\n"
            "raise SystemExit(1 if sys.argv[1]=='run' else 0)\n")
        docker.chmod(0o755)
        self.environment = {key: value for key, value in os.environ.items() if not key.startswith("SNPSTER_")}
        self.environment.update(PATH=str(tools) + os.pathsep + os.environ["PATH"],
                                DOCKER_TEST_CALLS=str(self.calls))

    def launch(self, *, resume=False, profile=None, reuse_manifest=None, staged=False, guided=False):
        environment = dict(self.environment)
        if profile is not None:
            environment["SNPSTER_JAVA_PROFILE"] = profile
        if reuse_manifest is not None:
            environment["SNPSTER_REUSE_MANIFEST"] = str(reuse_manifest)
        if staged:
            environment["SNPSTER_STAGE_SOURCE_RUN"] = str(self.root)
        command = ["bash", str(self.source / "launch.sh"), str(self.config), str(self.run),
                   "--guided-followup" if guided else "--staged-development" if staged else "--smoke-first"]
        if resume:
            command.append("--resume")
        return subprocess.run(command, env=environment, capture_output=True, text=True, timeout=30)

    def test_runtime_profile_persists_without_changing_frozen_source(self):
        result = self.launch(profile="c1")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertEqual((self.run / "java-profile").read_text().strip(), "c1")
        overlay = self.run / "runtime/java-c1-pinned.config"
        self.assertEqual(overlay.read_bytes(), (self.source / "java_c1.config").read_bytes())
        manifest = (self.run / "source.sha256").read_bytes()
        (self.run / "smoke/launch/.nextflow").mkdir(parents=True)
        result = self.launch(resume=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("JVM runtime profile: c1", result.stdout)
        self.assertEqual(manifest, (self.run / "source.sha256").read_bytes())
        calls = [json.loads(line) for line in self.calls.read_text().splitlines()]
        runs = [args for args in calls if args[0] == "run"]
        self.assertEqual(len(runs), 2)
        self.assertTrue(all("JAVA_PROFILE=c1" in args for args in runs))
        self.assertIn('-C "$runtime_config"', runs[-1][-1])
        self.assertIn('runtime_config="$runtime_config,$EXPERIMENT/runtime/java-c1-pinned.config"', runs[-1][-1])

    def test_runtime_downgrade_is_rejected(self):
        self.launch(profile="c1")
        result = self.launch(resume=True, profile="default")
        self.assertEqual(result.returncode, 2)
        self.assertIn("Do not mix", result.stderr)

    def test_invalid_profile_is_rejected_before_launch(self):
        result = self.launch(profile="unknown")
        self.assertEqual(result.returncode, 2)
        self.assertFalse(self.run.exists())
        self.assertFalse(self.calls.exists())

    def test_changed_runtime_overlay_is_rejected(self):
        self.launch(profile="c1")
        (self.run / "smoke/launch/.nextflow").mkdir(parents=True)
        (self.run / "runtime/java-c1-pinned.config").write_text("changed\n")
        result = self.launch(resume=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("profile changed", result.stderr)

    def test_reuse_manifest_is_frozen_and_cannot_change_on_resume(self):
        reuse_manifest = self.root / "reuse.json"
        reuse_manifest.write_text('{"java_profile":"c1"}\n')
        result = self.launch(profile="c1", reuse_manifest=reuse_manifest)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        saved = self.run / "runtime/reuse-manifest.json"
        self.assertEqual(saved.read_bytes(), reuse_manifest.read_bytes())
        (self.run / "smoke/launch/.nextflow").mkdir(parents=True)
        result = self.launch(resume=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        reuse_manifest.write_text('{"java_profile":"c1","changed":true}\n')
        result = self.launch(resume=True, reuse_manifest=reuse_manifest)
        self.assertEqual(result.returncode, 1)
        self.assertIn("cannot add or change", result.stderr)

    def test_staged_policy_and_source_are_frozen_and_graphical(self):
        result = self.launch(profile="c1", staged=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertEqual((self.run / "staged-development").read_text().strip(), "true")
        self.assertEqual((self.run / "runtime/staged-source-run").read_text().strip(), str(self.root))
        self.assertTrue((self.run / "source/staged.nf").is_file())
        calls = [json.loads(line) for line in self.calls.read_text().splitlines()]
        controller = next(call for call in calls if call[0] == "run")
        self.assertIn("STAGED_DEVELOPMENT=true", controller)
        self.assertIn("-t", controller)
        self.assertIn("-ansi-log true", controller[-1])
        self.assertIn('metadata planning plan', controller[-1])
        self.assertIn('metadata advancement advance', controller[-1])
        self.assertIn('audit genome', controller[-1])
        self.assertIn('cleanup --source-run', controller[-1])
        (self.run / "planning/launch/.nextflow").mkdir(parents=True)
        resumed = self.launch(resume=True, staged=True)
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        (self.run / "reports/scoped-cache-cleanup.json").write_text("{}\n")
        completed = self.launch(resume=True, staged=True)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("already passed both audits", completed.stderr)

    def test_staged_recovery_requires_matching_c1_profile_before_freezing(self):
        result = self.launch(staged=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires SNPSTER_JAVA_PROFILE=c1", result.stderr)
        self.assertFalse(self.run.exists())

    def test_guided_policy_is_frozen_and_has_no_genome_expansion(self):
        result = self.launch(profile="c1", guided=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertEqual((self.run / "guided-followup").read_text().strip(), "true")
        self.assertTrue((self.run / "source/guided_followup.nf").is_file())
        self.assertIn("No test or genome expansion", result.stdout)
        calls = [json.loads(line) for line in self.calls.read_text().splitlines()]
        controller = next(call for call in calls if call[0] == "run")
        self.assertIn("GUIDED_FOLLOWUP=true", controller)
        script = controller[-1]
        self.assertIn('guided_metadata planning plan', script)
        self.assertIn('pipeline screen "$EXPERIMENT/planning/design/config.json"', script)
        self.assertIn('guided_metadata reporting report', script)
        self.assertIn('guided_followup.py cleanup', script)
        (self.run / "planning/launch/.nextflow").mkdir(parents=True)
        resumed = self.launch(resume=True, guided=True)
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        wrong_policy = self.launch(resume=True, staged=True)
        self.assertEqual(wrong_policy.returncode, 2)

    def test_guided_runtime_must_match_verified_phasing_runtime(self):
        result = self.launch(guided=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires SNPSTER_JAVA_PROFILE=c1", result.stderr)
        self.assertFalse(self.run.exists())


if __name__ == "__main__":
    unittest.main()
