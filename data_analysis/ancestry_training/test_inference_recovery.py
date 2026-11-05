from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import core
import inference_recovery as recovery


class SingletonRecovery(unittest.TestCase):
    def test_incomplete_singleton_is_not_success_shaped(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(recovery.projection, "capture", return_value=None):
                self.assertIsNone(recovery.checked_singleton(Path(temporary), {}, {}, "s", np.zeros((2, 2))))

    def test_completed_singleton_requires_exact_genotypes_keys_and_axis(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            keys = [("22", 10, "A", "G"), ("22", 20, "A", "G")]
            calls = np.asarray([[[0, 1], [1, 0]]], dtype=np.int8)
            core.write_vcf(root / "s.vcf.gz", keys, ["s"], calls)
            (root / "s.model").write_text("# ancestry axis\nB A\n")
            meta = {"chromosome": 22, "keys": keys, "span": [0, 30], "sources": ["A", "B"]}
            with patch.object(recovery.projection, "capture", return_value={"files": {}}):
                with patch.object(recovery.models, "flare_posterior") as posterior:
                    recovery.checked_singleton(root, meta, {}, "s", calls[0])
                    self.assertEqual(posterior.call_args.args[2], ["B", "A"])
                    with self.assertRaisesRegex(ValueError, "genotype"):
                        recovery.checked_singleton(root, meta, {}, "s", 1 - calls[0])
                    (root / "s.model").write_text("A C\n")
                    with self.assertRaisesRegex(ValueError, "axis"):
                        recovery.checked_singleton(root, meta, {}, "s", calls[0])

    def test_snapshot_requires_all_hashes_and_becomes_work_independent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = root / "work"
            original.mkdir()
            (original / "s.model").write_text("A B\n")
            record = {"directory": str(original), "files": {"s.model": core.digest(original / "s.model")}}
            saved = recovery.snapshot(record, root / "reports/saved")
            (original / "s.model").write_text("changed\n")
            self.assertEqual(core.digest(Path(saved["directory"]) / "s.model"), saved["files"]["s.model"])
            with self.assertRaisesRegex(ValueError, "changed"):
                recovery.snapshot(record, root / "reports/wrong")

    def test_recovery_rejects_running_or_wrong_root_stop_proof(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = root / "stop.json"
            core.write(receipt, {"root": str(root), "controller_after": {"State": {"Running": True}},
                                 "owned_workers_after": []})
            with self.assertRaisesRegex(ValueError, "stopped"):
                recovery.build(root, "screen", receipt, root / "reports/recovered")
            core.write(receipt, {"root": str(root / "other"), "controller_after": {"State": {"Running": False}},
                                 "owned_workers_after": []})
            with self.assertRaisesRegex(ValueError, "stopped"):
                recovery.build(root, "screen", receipt, root / "reports/recovered")


class LauncherCapacity(unittest.TestCase):
    def test_every_launch_command_sets_launch_pool_to_task_budget(self):
        source = Path(__file__).parent
        for name in ("launch.sh", "launch_final.sh", "launch_population.sh"):
            text = (source / name).read_text()
            commands = [part for part in text.split("nextflow -log")[1:]]
            self.assertTrue(commands)
            for command in commands:
                self.assertIn('-pool-size "$BUDGET_CPUS"', command.split("\n        }\n", 1)[0])


if __name__ == "__main__":
    unittest.main()
