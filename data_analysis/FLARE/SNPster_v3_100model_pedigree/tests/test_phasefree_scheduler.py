import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import schedule_phasefree_arrays as scheduler


class PhasefreeSchedulerTests(unittest.TestCase):
    def test_incomplete_preparation_is_archived_without_changing_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            directory = out / "chr1"
            directory.mkdir()
            (directory / "partial.vcf.gz").write_bytes(b"interrupted")
            result = scheduler.archive_incomplete_preparation(out, 1)
            self.assertFalse(directory.exists())
            self.assertEqual(Path(result["archive"]).joinpath("partial.vcf.gz").read_bytes(),
                             b"interrupted")
            self.assertIsNone(scheduler.archive_incomplete_preparation(out, 1))

    def test_completed_preparation_is_kept_and_inconsistent_cases_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            directory = out / "chr1"
            directory.mkdir()
            (directory / "prepared.json").write_text("{}")
            self.assertIsNone(scheduler.archive_incomplete_preparation(out, 1))
            (directory / "prepared.json").unlink()
            case = directory / "person"
            case.mkdir()
            (case / "completion.json").write_text("{}")
            with self.assertRaisesRegex(ValueError, "completed work"):
                scheduler.archive_incomplete_preparation(out, 1)
            self.assertTrue((case / "completion.json").exists())

    def test_child_failure_is_persistent_and_never_success_shaped(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            state = {"status": "running"}
            with self.assertRaisesRegex(RuntimeError, "failed"):
                scheduler.run_batch([
                    ("good", [sys.executable, "-c", "print('completed')"]),
                    ("bad", [sys.executable, "-c", "raise SystemExit(7)"]),
                ], out, state)
            report = json.loads((out / "durable_scheduler.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["failed_jobs"], {"bad": 7})
            self.assertIn("completed", (out / "durable_good.log").read_text())


if __name__ == "__main__":
    unittest.main()
