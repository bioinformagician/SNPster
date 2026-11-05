import copy
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import core
from make_smoke import fixture
import models
import reuse


class InferenceChunks(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.panel = self.root / "panel"
        self.panel.mkdir()
        self.samples = [{"sample_id": f"person_{index}"} for index in range(5)]
        self.truth = np.asarray([[1., 0.], [0., 1.], [.5, .5], [.25, .75], [.75, .25]])
        self.meta = {
            "chromosome": 22, "route": {"id": "observed"}, "span": [0, 100],
            "sources": ["B", "A"], "config": {}, "training": [{"sample_id": "train"}],
            "queries": self.samples,
        }
        core.write(self.panel / "panel.json", self.meta)
        np.savez(self.panel / "data.npz", truth=np.concatenate([[[1., 0.]], self.truth]))
        self.inputs = {"panel.json": core.digest(self.panel / "panel.json"),
                       "data.npz": core.digest(self.panel / "data.npz"), "query.npy": "frozen-query"}

    def chunk(self, name, start, end):
        folder = self.root / name
        folder.mkdir()
        for sample in self.samples[start:end]:
            for suffix in (".model", ".anc.vcf.gz"):
                (folder / f"{sample['sample_id']}{suffix}").write_bytes(b"native evidence")
        artifacts = {path.name: core.digest(path) for path in folder.iterdir()}
        np.savez(folder / "predictions.npz", values=self.truth[start:end], truth=self.truth[start:end])
        core.write(folder / "prediction.json", {
            **{key: self.meta[key] for key in ("chromosome", "route", "span", "sources", "config")},
            "id": "flare", "family": "flare", "license_track": "permissive_software",
            "components": ["B", "A"], "samples": self.samples[start:end],
            "inference_input_sha256": self.inputs,
            "inference_chunk": {"start": start, "end": end, "total": len(self.samples)},
            "artifact_sha256": artifacts, "seconds": 1., "max_native_rounding_correction": .01,
            "native_probability_policy": {}, "sha256": core.digest(folder / "predictions.npz"),
        })
        return folder

    def test_shuffled_chunks_restore_complete_person_order_and_native_evidence(self):
        chunks = [self.chunk("last", 4, 5), self.chunk("first", 0, 2), self.chunk("middle", 2, 4)]
        output = self.root / "complete"
        models.collect_inference(self.panel, chunks, output)
        with np.load(output / "predictions.npz") as data:
            np.testing.assert_array_equal(data["values"], self.truth)
            np.testing.assert_array_equal(data["truth"], self.truth)
        record = core.read(output / "prediction.json")
        self.assertEqual(record["samples"], self.samples)
        self.assertNotIn("inference_chunk", record)
        self.assertEqual([chunk["start"] for chunk in record["inference_chunks"]], [0, 2, 4])
        self.assertEqual(record["seconds"], 3.)
        for name, expected in record["artifact_sha256"].items():
            self.assertEqual(core.digest(output / name), expected)

    def test_missing_overlapping_and_duplicate_chunks_fail(self):
        chunks = [self.chunk("first", 0, 2), self.chunk("last", 2, 5)]
        variants = ([chunks[0]], [chunks[1]], [chunks[0], chunks[0], chunks[1]],
                    [chunks[0], self.chunk("overlap", 1, 5)])
        for index, folders in enumerate(variants):
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, "Incomplete|overlapping|duplicate"):
                models.collect_inference(self.panel, folders, self.root / f"bad{index}")

    def test_single_person_chunks_restore_order_and_preserve_native_evidence(self):
        chunks = [self.chunk(f"singleton{index}", index, index + 1) for index in range(len(self.samples))]
        output = self.root / "singletons-complete"
        models.collect_inference(self.panel, list(reversed(chunks)), output)
        with np.load(output / "predictions.npz") as data:
            np.testing.assert_array_equal(data["values"], self.truth)
            np.testing.assert_array_equal(data["truth"], self.truth)
        record = core.read(output / "prediction.json")
        self.assertEqual(record["samples"], self.samples)
        self.assertEqual(len(record["inference_chunks"]), len(self.samples))
        self.assertTrue(all(chunk["end"] - chunk["start"] == 1 for chunk in record["inference_chunks"]))
        for name, expected in record["artifact_sha256"].items():
            self.assertEqual(core.digest(output / name), expected)

    def test_changed_input_model_person_and_axis_fail(self):
        chunks = [self.chunk("first", 0, 2), self.chunk("last", 2, 5)]
        original = core.read(chunks[1] / "prediction.json")
        for index, key in enumerate(("inference_input_sha256", "id", "samples", "components", "config")):
            changed = copy.deepcopy(original)
            if key == "inference_input_sha256":
                changed[key]["query.npy"] = "different-query"
            elif key == "samples":
                changed[key] = list(reversed(changed[key]))
            elif key == "components":
                changed[key] = list(reversed(changed[key]))
            else:
                changed[key] = "different"
            core.write(chunks[1] / "prediction.json", changed)
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "mismatch"):
                models.collect_inference(self.panel, chunks, self.root / f"bad{index}")

    def test_changed_predictions_truth_and_native_evidence_fail(self):
        for index, kind in enumerate(("hash", "truth", "native", "missing", "path")):
            chunks = [self.chunk(f"first{index}", 0, 2), self.chunk(f"last{index}", 2, 5)]
            record = core.read(chunks[1] / "prediction.json")
            if kind in ("hash", "truth"):
                np.savez(chunks[1] / "predictions.npz", values=self.truth[2:], truth=self.truth[2:, ::-1])
                if kind == "truth":
                    record["sha256"] = core.digest(chunks[1] / "predictions.npz")
            elif kind == "native":
                (chunks[1] / "person_2.model").write_text("changed")
            elif kind == "missing":
                del record["artifact_sha256"]["person_2.model"]
            else:
                record["artifact_sha256"]["../escape"] = "unexpected"
            core.write(chunks[1] / "prediction.json", record)
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                models.collect_inference(self.panel, chunks, self.root / f"bad{index}")

    def test_invalid_ranges_fail_before_native_execution(self):
        fitted = self.root / "fitted"
        fitted.mkdir()
        core.write(fitted / "model.json", {"family": "flare"})
        for index, (start, end) in enumerate(((-1, 1), (0, 0), (4, 6), (True, 2), (1, None))):
            with self.subTest(start=start, end=end), self.assertRaisesRegex(ValueError, "range"):
                models.infer(self.panel, self.root, self.root, fitted, self.root / f"bad{index}",
                             query_start=start, query_end=end)
        core.write(fitted / "model.json", {"family": "gnomix"})
        with self.assertRaisesRegex(ValueError, "only for singleton FLARE"):
            models.infer(self.panel, self.root, self.root, fitted, self.root / "nonflare",
                         query_start=0, query_end=2)

    def test_completed_stage_adoption_checks_exact_scope_and_all_hashes(self):
        frozen = self.root / "frozen"
        frozen.mkdir()
        core.write(frozen / "config.json", self.meta["config"])
        reuse_path = self.root / "reuse.json"
        record = {
            "version": 1, "config": self.meta["config"], "source_run": "original",
            "assets": {"prepare/observed/22": {"directory": str(self.panel), "files": reuse.inventory(self.panel)}},
        }
        core.write(reuse_path, record)
        args = SimpleNamespace(stage="prepare", reuse_manifest=reuse_path, frozen=frozen, route="observed",
                               chromosome=22, output=self.root / "adopted")
        self.assertTrue(reuse.adopt_stage(args))
        self.assertEqual(core.digest(args.output / "data.npz"), core.digest(self.panel / "data.npz"))
        self.assertEqual(core.read(args.output / "execution_reuse.json")["status"], "adopted_not_refitted")
        args.output = self.root / "changed"
        (self.panel / "data.npz").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "Reusable output changed"):
            reuse.adopt_stage(args)
        record["config"] = {"different": True}
        core.write(reuse_path, record)
        with self.assertRaisesRegex(ValueError, "scientific scope"):
            reuse.adopt_stage(args)


@unittest.skipUnless(os.environ.get("SNPSTER_TEST_FLARE_JAR"), "Set SNPSTER_TEST_FLARE_JAR for real FLARE equivalence")
class NativeInferenceChunks(unittest.TestCase):
    def test_real_fixed_and_em_flare_are_identical_to_serial_singletons(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            specs = [{"id": name, "family": "flare", "params": {"em": em, "gen": 10, "states": 20}}
                     for name, em in (("fixed", False), ("em", True))]
            config = fixture(root / "inputs", {"flare_jar": os.environ["SNPSTER_TEST_FLARE_JAR"]}, models=specs)
            core.freeze(root / "inputs/experiment.json", root / "frozen")
            core.prepare(root / "frozen", "imputed", 22, root / "panel")
            core.phase(root / "panel", root / "phased")
            base_id = next(iter(models.plan(config)["bases"]))
            models.fit_base(root / "panel", base_id, root / "base")
            meta = core.read(root / "panel/panel.json")
            for spec in specs:
                fitted = root / spec["id"]
                models.fit_model(root / "panel", root / "base", spec["id"], fitted)
                serial = root / f"{spec['id']}_serial"
                models.infer(root / "panel", root / "phased", root / "base", fitted, serial)
                commands = []
                native_run = models.run

                def checked_run(command, log):
                    commands.append(command)
                    native_run(command, log)

                chunks = []
                with patch("models.run", side_effect=checked_run):
                    for start in range(0, len(meta["queries"]), 5):
                        chunk = root / f"{spec['id']}_{start}"
                        models.infer(root / "panel", root / "phased", root / "base", fitted, chunk,
                                     query_start=start, query_end=min(start + 5, len(meta["queries"])))
                        chunks.append(chunk)
                merged = root / f"{spec['id']}_merged"
                models.collect_inference(root / "panel", chunks[::-1], merged)
                with np.load(serial / "predictions.npz") as left, np.load(merged / "predictions.npz") as right:
                    np.testing.assert_array_equal(left["values"], right["values"])
                    np.testing.assert_array_equal(left["truth"], right["truth"])
                self.assertEqual(len(commands), len(meta["queries"]))
                self.assertTrue(all("nthreads=1" in command and f"seed={config['seed']}" in command
                                    for command in commands))
                self.assertEqual(core.read(serial / "prediction.json")["samples"],
                                 core.read(merged / "prediction.json")["samples"])
                native = {}
                for sample in meta["queries"][:3]:
                    identity = sample["sample_id"]
                    native[identity] = {
                        "directory": str(serial),
                        "files": {f"{identity}{suffix}": core.digest(serial / f"{identity}{suffix}")
                                  for suffix in reuse.NATIVE_SUFFIXES if (serial / f"{identity}{suffix}").is_file()},
                    }
                key = f"{meta['route']['id']}/{meta['chromosome']}/{spec['id']}"
                reuse_record = {
                    "version": 1, "config": meta["config"],
                    "native": {key: {"samples": native, "input_sha256": reuse.inference_inputs(
                        root / "panel", root / "phased", root / "base", fitted)}},
                }
                reuse_path = root / f"{spec['id']}_reuse.json"
                core.write(reuse_path, reuse_record)
                resumed = root / f"{spec['id']}_resumed"
                commands.clear()
                with patch("models.run", side_effect=checked_run):
                    models.infer(root / "panel", root / "phased", root / "base", fitted, resumed,
                                 query_start=0, query_end=len(meta["queries"]), reuse_manifest=reuse_path)
                self.assertEqual(len(commands), len(meta["queries"]) - 3)
                with np.load(serial / "predictions.npz") as left, np.load(resumed / "predictions.npz") as right:
                    np.testing.assert_array_equal(left["values"], right["values"])
                self.assertEqual(core.read(resumed / "prediction.json")["reused_singletons"], list(native))
                reuse_record["native"][key]["input_sha256"]["query.npy"] = "different"
                core.write(reuse_path, reuse_record)
                with self.assertRaisesRegex(ValueError, "different panel, phase"):
                    models.infer(root / "panel", root / "phased", root / "base", fitted,
                                 root / f"{spec['id']}_invalid", query_start=0, query_end=1, reuse_manifest=reuse_path)


if __name__ == "__main__":
    unittest.main()
