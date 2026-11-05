import csv
import json
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch
import uuid

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import model_family_gnomix as adapter


class GnomixAdapterTests(unittest.TestCase):
    def setUp(self):
        # Explicit project-relative test scratch, never the system temp directory.
        self.work = Path.cwd() / (".gnomix-test-" + uuid.uuid4().hex)
        self.work.mkdir()
        self.addCleanup(shutil.rmtree, self.work)

    def tsv(self, name, columns, rows):
        path = self.work / name
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)
        return path

    def layout(self, count=640, candidate="gnomix_default_020"):
        keys = [("22", 100 + i * 10, "A", "G") for i in range(count)]
        points = np.array([[100., 0.], [100 + (count - 1) * 10, 64.]])
        return keys, adapter.window_layout(keys, points, [90, 110 + (count - 1) * 10],
                                           adapter.config_for(candidate))

    def indexed_vcf(self, name="input", order=("extra", "a", "b"), variants=None, chroms=(22,)):
        path = self.work / (name + ".vcf")
        rows = variants if variants is not None else [
            (100, "A", "G", ["1|1", "0|1", "1|0"]),
            (200, "C", "T", ["0|0", "1|0", "0|1"])]
        with path.open("w") as handle:
            handle.write('##fileformat=VCFv4.2\n')
            for chrom in chroms:
                handle.write(f'##contig=<ID={chrom},length=10000>\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                         + "\t".join(order) + "\n")
            for chrom in chroms:
                for pos, ref, alt, genotypes in rows:
                    handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
                                 + "\t".join(genotypes) + "\n")
        gz = self.work / (name + ".vcf.gz")
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(gz), str(path)], check=True,
                       capture_output=True)
        subprocess.run(["bcftools", "index", str(gz)], check=True, capture_output=True)
        return gz

    def test_exact_six_configs_and_smoother_genetic_span(self):
        for mode in ("fast", "default", "best"):
            for suffix, size, smoother in (("020", .2, 75), ("050", .5, 31)):
                config = adapter.config_for(f"gnomix_{mode}_{suffix}")
                self.assertEqual(config["window_size_cM"], size)
                self.assertEqual(config["smooth_size"], smoother)
                self.assertEqual(config["seed"], 20261013)
                self.assertFalse(config["calibrate"])
                self.assertFalse(config["retrain_base"])
                self.assertLessEqual(abs(size * smoother - 15), .5)
        with self.assertRaises(ValueError):
            adapter.config_for("gnomix_large_020")

    def test_upstream_window_conversion_and_remainder(self):
        keys = [("22", i + 1, "A", "G") for i in range(651)]
        layout = adapter.window_layout(keys, np.array([[1, 2], [651, 64]]), [0, 652],
                                       adapter.config_for("gnomix_default_050"))
        self.assertEqual(layout["M"], round(.5 * 651 / 64))
        self.assertEqual(layout["ends"][-1], 651)
        self.assertEqual(layout["marker_weights"].sum(), 651)
        self.assertEqual(layout["physical_weights"].sum(), 652)
        self.assertGreater(layout["marker_weights"][-1], layout["M"])

    def test_sparse_short_and_out_of_span_windows_fail(self):
        keys, _ = self.layout()
        with self.assertRaisesRegex(ValueError, "sparse/short"):
            adapter.window_layout(keys, np.array([[100, 0], [6490, 1.]]), [90, 6500],
                                  adapter.config_for("gnomix_default_020"))
        with self.assertRaisesRegex(ValueError, "density"):
            adapter.window_layout(keys, np.array([[100, 0], [6490, 999.]]), [90, 6500],
                                  adapter.config_for("gnomix_fast_020"))
        with self.assertRaisesRegex(ValueError, "span"):
            adapter.window_layout(keys, np.array([[100, 0], [6490, 64.]]), [100, 6490],
                                  adapter.config_for("gnomix_default_020"))
        with self.assertRaisesRegex(ValueError, "bracket"):
            adapter.window_layout(keys, np.array([[101, 0], [6490, 64.]]), [90, 6500],
                                  adapter.config_for("gnomix_default_020"))

    def test_phasing_and_requested_sample_order_are_preserved(self):
        path = self.indexed_vcf()
        keys = [("22", 100, "A", "G"), ("22", 200, "C", "T")]
        x = adapter.load_haplotypes(path, ["b", "a"], 22, keys)
        np.testing.assert_array_equal(x, [[1, 0], [0, 1], [0, 1], [1, 0]])
        with self.assertRaisesRegex(ValueError, "Missing"):
            adapter.load_haplotypes(path, ["absent"], 22, keys)
        with self.assertRaisesRegex(ValueError, "mismatch"):
            adapter.load_haplotypes(path, ["a"], 22, [("22", 100, "G", "A")])
        with self.assertRaisesRegex(ValueError, "missing 1"):
            adapter.load_haplotypes(path, ["a"], 22, keys + [("22", 300, "A", "T")])

    def test_unphased_and_multiallelic_GT_are_not_silently_repaired(self):
        self.assertEqual(adapter.phased_gt(".|1"), [2, 1])
        for value in ("0/1", "0/0", "0", "0|2", "0|1|1", "./."):
            with self.assertRaises(ValueError):
                adapter.phased_gt(value)

    def test_reference_filter_excludes_unlabelled_samples(self):
        path = self.indexed_vcf(variants=[
            (100, "A", "G", ["1|1", "0|0", "0|0"]),
            (200, "C", "T", ["0|0", "1|0", "0|1"])])
        keys = adapter.reference_markers(path, ["a", "b"], 22, dict(min_maf=.1, min_mac=1))
        self.assertEqual(keys, [("22", 200, "C", "T")])
        # Evaluation is not refiltered, even when all requested alleles are identical.
        x = adapter.load_haplotypes(path, ["extra"], 22, keys)
        np.testing.assert_array_equal(x, [[0], [0]])

    def test_exact_shared_marker_check(self):
        path = self.work / "markers.tsv"
        path.write_text("22\t100\tA\tG\n22\t200\tC\tT\n")
        adapter.check_shared_markers(path, 22, [("22", 100, "A", "G"), ("22", 200, "C", "T")])
        with self.assertRaises(ValueError):
            adapter.check_shared_markers(path, 22, [("22", 100, "A", "G")])

    def test_segments_use_direct_one_based_position_and_half_open_end(self):
        path = self.tsv("segments.tsv",
                        ["sample_id", "chrom", "haplotype", "start_bp", "end_bp", "source_group"],
                        [["a", 22, 1, 0, 200, "EUR"], ["a", 22, 1, 200, 301, "SSA"],
                         ["a", 22, 2, 0, 301, "EAS_SEA"],
                         ["b", 22, 1, 0, 301, "OCE"], ["b", 22, 2, 0, 301, "SAS"]])
        labels = adapter.segment_labels(path, ["b", "a"], 22, np.array([100, 200, 300]))
        np.testing.assert_array_equal(labels, [[5, 5, 5], [2, 2, 2], [0, 6, 6], [1, 1, 1]])

    def test_missing_overlapping_and_invalid_haplotype_segments_fail(self):
        columns = ["sample_id", "chrom", "haplotype", "start_bp", "end_bp", "source_group"]
        for rows, text in [
                ([["a", 22, 1, 0, 301, "EUR"]], "Incomplete"),
                ([["a", 22, 1, 0, 301, "EUR"], ["a", 22, 1, 100, 200, "SAS"]], "Overlapping"),
                ([["a", 22, 0, 0, 301, "EUR"]], "Invalid")]:
            path = self.tsv("segments.tsv", columns, rows)
            with self.assertRaisesRegex(ValueError, text):
                adapter.segment_labels(path, ["a"], 22, np.array([100, 200, 300]))

    def test_window_truth_ties_match_upstream_mode_and_include_remainder(self):
        layout = dict(starts=np.array([0, 2]), ends=np.array([2, 5]))
        labels = np.array([[6, 0, 1, 1, 6], [2, 2, 3, 4, 4]])
        np.testing.assert_array_equal(adapter.window_labels(labels, layout), [[0, 1], [2, 4]])

    def test_primary_marker_counts_and_secondary_physical_weights_differ(self):
        p = np.zeros((2, 2, 7))
        p[:, 0, 0], p[:, 1, 6] = 1, 1
        layout = dict(W=2, marker_weights=np.array([2, 3]), physical_weights=np.array([90, 10]))
        primary, physical = adapter.aggregate_probabilities(p, layout)
        np.testing.assert_array_equal(primary[0], [2, 0, 0, 0, 0, 0, 3])
        np.testing.assert_array_equal(physical[0], [90, 0, 0, 0, 0, 0, 10])
        with self.assertRaisesRegex(ValueError, "Unpaired"):
            adapter.aggregate_probabilities(p[:1], layout)

    def test_crf_probability_axis_is_reordered(self):
        class CRF:
            classes_ = ["6", "0", "1", "2", "3", "4", "5"]
        class Smooth:
            model = CRF()
            def predict_proba(self, base):
                return np.tile(np.arange(1, 8) / 28, (2, 2, 1))
        class Model:
            smooth = Smooth()
            W = 2
        result = adapter.ordered_probabilities(Model(), np.empty((2, 2, 7)))
        np.testing.assert_allclose(result[0, 0], np.array([2, 3, 4, 5, 6, 7, 1]) / 28)

    def test_output_schema_and_probability_normalization(self):
        path = self.work / "predictions.tsv"
        values = {("pure", "test"): np.arange(1, 8, dtype=float)}
        adapter.write_predictions(path, values, 28)
        rows = adapter.records(path)
        self.assertEqual(list(rows[0]), ["dataset", "sample_id"] + ["pred_" + s for s in adapter.SOURCES])
        self.assertAlmostEqual(sum(float(rows[0]["pred_" + s]) for s in adapter.SOURCES), 1.)

    def test_manifest_rejects_evaluation_training_overlap(self):
        ref = self.tsv("ref.tsv", ["sample_id", "source_group"],
                       [[f"ref_{s}", s] for s in adapter.SOURCES])
        pure = self.tsv("train.tsv", ["sample_id", "source_group"],
                        [[f"train_{s}", s] for s in adapter.SOURCES])
        empty = dict(batches=[])
        manifest = dict(sources=list(adapter.SOURCES), reference=dict(labels_tsv=str(ref)),
                        pure_train=dict(labels_tsv=str(pure)), pedigree_train=empty,
                        evaluation={key: empty for key in ("pure", "pedigree", "real")})
        counts = adapter.validate_splits(manifest)
        self.assertEqual(counts["reference"], 7)
        manifest["evaluation"]["pure"] = dict(batches=[dict(samples_tsv=str(ref))])
        with self.assertRaisesRegex(ValueError, "overlap"):
            adapter.validate_splits(manifest)

    def test_all_six_real_upstream_models_fit_and_roundtrip(self):
        """Tiny pure-only fixtures test APIs, not synthetic benchmark mosaics."""
        gnomix_class, windowed_class = adapter.load_upstream()
        rng = np.random.default_rng(17)
        x = rng.integers(0, 2, size=(28, 640), dtype=np.int8)
        pure_classes = np.repeat(np.arange(7), 4)
        x_path, y_path = self.work / "x.npy", self.work / "y.npy"
        np.save(x_path, x)
        np.save(y_path, np.repeat(pure_classes[:, None], 640, axis=1))
        for mode in ("fast", "default", "best"):
            for suffix in ("020", "050"):
                with self.subTest(mode=mode, suffix=suffix):
                    candidate = f"gnomix_{mode}_{suffix}"
                    keys, layout = self.layout(candidate=candidate)
                    config = adapter.config_for(candidate)
                    model = adapter.initialize_model(config, keys, layout, 1, gnomix_class)
                    self.assertEqual(type(model).__module__, "src.model")
                    data = windowed_class([str(x_path)], [str(y_path)], 640, layout["M"])
                    model.base.train_windowed(data, verbose=False)
                    y = np.repeat(pure_classes[:, None], layout["W"], axis=1)
                    b = model.base.predict_proba(x)
                    model.smooth.train(b, y)
                    before = adapter.ordered_probabilities(model, b)
                    loaded = pickle.loads(pickle.dumps(model))
                    after = adapter.ordered_probabilities(loaded, loaded.base.predict_proba(x))
                    np.testing.assert_allclose(before, after, rtol=0, atol=1e-12)
                    self.assertFalse(model.calibrate)
                    self.assertEqual(model.base.models[0].random_state, 20261013)
                    self.assertEqual(model.smooth.S, config["smooth_size"])

    def test_evaluation_does_not_read_evaluation_segments_or_truth(self):
        keys, layout = self.layout()
        manifest = dict(evaluation={name: dict(batches=[]) for name in ("pure", "pedigree", "real")})
        manifest["evaluation"]["pure"]["batches"] = [
            dict(vcf="query", samples_tsv="ids", segments_tsv="DO_NOT_READ")]
        class Base:
            def predict_proba(self, x):
                return np.ones((2, layout["W"], 7)) / 7
        class Smooth:
            model = type("Classes", (), dict(classes_=np.arange(7)))()
            def predict_proba(self, b):
                return b
        model = type("Model", (), dict(base=Base(), smooth=Smooth(), W=layout["W"]))()
        bundle = dict(config={}, keys=keys, layout=layout, upstream_sha=adapter.SHA,
                      chromosome=22, model=model)
        with patch.object(adapter, "sample_records", return_value=[dict(sample_id="query")]) as read:
            with patch.object(adapter, "load_haplotypes", return_value=np.zeros((2, 640))):
                primary, physical, delta = adapter.evaluate_chromosome(manifest, bundle, bundle, 64)
        read.assert_called_once_with("ids")
        self.assertEqual(len(primary), 1)
        self.assertEqual(len(physical), 1)
        self.assertEqual(delta, 0)

    def cohort_fixture(self, chroms=(22,)):
        rng = np.random.default_rng(26)
        roles = {role: [role + "_" + source for source in adapter.SOURCES]
                 for role in ("ref", "train", "eval")}
        names = roles["ref"] + roles["train"] + roles["eval"] + ["pedtrain", "dev", "real"]
        rows = [(100 + 10 * i, "A", "G",
                 [f"{a}|{b}" for a, b in rng.integers(0, 2, size=(len(names), 2))])
                for i in range(640)]
        vcf = self.indexed_vcf("cohort", names, rows, chroms)
        paths = {}
        for role in ("ref", "train", "eval"):
            paths[role] = self.tsv(role + ".tsv", ["sample_id", "source_group"],
                                   [[name, source] for name, source in zip(roles[role], adapter.SOURCES)])
        for role in ("pedtrain", "dev", "real"):
            paths[role] = self.tsv(role + ".tsv", ["sample_id"], [[role]])
        segments = self.tsv("truth.tsv",
                            ["sample_id", "chrom", "haplotype", "start_bp", "end_bp", "source_group"],
                            [["pedtrain", chrom, h, 90, 6500, "EUR"]
                             for chrom in chroms for h in (1, 2)])
        map_path = self.work / "map"
        map_path.write_text("".join(
            f"{chrom}\tleft\t0\t100\n{chrom}\tright\t{32 if len(chroms) > 1 and chrom == 22 else 64}\t6490\n"
            for chrom in chroms))
        def batch(role):
            return dict(vcf=str(vcf), samples_tsv=str(paths[role]))
        manifest = dict(
            sources=list(adapter.SOURCES), seed=adapter.SEED,
            reference=dict(vcf=str(vcf), labels_tsv=str(paths["ref"])),
            pure_train=dict(vcf=str(vcf), labels_tsv=str(paths["train"])),
            pedigree_train=dict(batches=[dict(**batch("pedtrain"), segments_tsv=str(segments))]),
            evaluation={name: dict(batches=[batch(role)]) for name, role in
                        (("pure", "eval"), ("pedigree", "dev"), ("real", "real"))},
            reference_filters=dict(min_maf=0, min_mac=1), genetic_map=str(map_path),
            chromosome_spans={str(chrom): [90, 6500] for chrom in chroms})
        (self.work / "cohort.json").write_text(json.dumps(manifest))
        return manifest

    def test_complete_smoke_cli_core_and_separate_evaluation_roundtrip(self):
        manifest = self.cohort_fixture()
        directory = adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8)
        before = (directory / "predictions.tsv").read_text()
        audit = json.loads((directory / "output_audit.json").read_text())
        self.assertTrue(audit["smoke_only"])
        self.assertTrue(audit["model_roundtrip_verified"])
        self.assertEqual(audit["n_predictions"], 9)
        self.assertEqual(audit["n_markers"], 640)
        self.assertEqual(audit["chromosomes"], [22])
        self.assertIn(manifest["reference"]["vcf"] + ".csi", audit["input_sha256"])
        adapter.run(self.work, "gnomix_fast_050", True, stage="evaluate", ncores=1, batch_size=8)
        self.assertEqual(before, (directory / "predictions.tsv").read_text())
        complete_audit = (directory / "audit.json").read_text()
        with patch.object(adapter, "fit_chromosome", side_effect=AssertionError("must not refit")):
            adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8,
                        resume=True, verify_only=True)
        self.assertEqual(complete_audit, (directory / "audit.json").read_text())
        original = self.work / "revisions/model_family_gnomix.before_resume.py"
        original.parent.mkdir()
        source = "# Archived revision\n" + Path(adapter.__file__).read_text()
        original.write_text(source)
        legacy = json.loads(complete_audit)
        legacy.pop("execution_fingerprint")
        legacy["adapter_sha256"] = adapter.digest(original)
        for name in ("audit.json", "output_audit.json"):
            adapter.write_json(directory / name, legacy)
        adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8,
                    resume=True, verify_only=True)
        original.write_text(source.replace("model.smooth.train(b, y)",
                                           "model.smooth.train(b, y[::-1])"))
        legacy["adapter_sha256"] = adapter.digest(original)
        for name in ("audit.json", "output_audit.json"):
            adapter.write_json(directory / name, legacy)
        with self.assertRaisesRegex(ValueError, "implementation changed"):
            adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8,
                        resume=True, verify_only=True)
        for name in ("predictions.tsv", "predictions_physical.tsv", "output_audit.json"):
            (directory / name).unlink()
        paused = json.loads(complete_audit)
        paused["status"] = "paused"
        adapter.write_json(directory / "audit.json", paused)
        with patch.object(adapter, "fit_chromosome", side_effect=AssertionError("must not refit")):
            adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8, resume=True)
        self.assertEqual(before, (directory / "predictions.tsv").read_text())
        metadata_path = directory / "models/chr22.json"
        metadata = json.loads(metadata_path.read_text())
        metadata["model_sha256"] = "0" * 64
        adapter.write_json(metadata_path, metadata)
        with self.assertRaisesRegex(ValueError, "checksum"):
            adapter.run(self.work, "gnomix_fast_050", True, ncores=1, batch_size=8,
                        resume=True, verify_only=True)
        self.assertEqual(before, (directory / "predictions.tsv").read_text())

    def test_two_chromosome_pause_resume_matches_uninterrupted(self):
        self.cohort_fixture((21, 22))
        select = lambda *args: range(21, 23) if args == (1, 23) else range(*args)
        candidate = "gnomix_fast_050"
        with patch.object(adapter, "range", side_effect=select, create=True):
            directory = adapter.run(self.work, candidate, ncores=1, batch_size=8)
            before = {name: (directory / name).read_bytes()
                      for name in ("predictions.tsv", "predictions_physical.tsv")}
            for name in ("chr22.pkl", "chr22.json"):
                (directory / "models" / name).unlink()
            evaluate = adapter.evaluate_chromosome
            def request_after_checkpoint(*args):
                result = evaluate(*args)
                adapter.write_json(self.work / "pause_requested.json", {"status": "pause_requested"})
                return result
            with patch.object(adapter, "evaluate_chromosome", side_effect=request_after_checkpoint):
                with self.assertRaises(SystemExit) as stopped:
                    adapter.run(self.work, candidate, ncores=1, batch_size=8)
            self.assertEqual(stopped.exception.code, 75)
            paused = json.loads((directory / "audit.json").read_text())
            self.assertEqual(paused["status"], "paused")
            self.assertEqual([item["chromosome"] for item in paused["chromosome_audits"]], [21])
            self.assertFalse((directory / "predictions.tsv").exists())
            (self.work / "pause_requested.json").unlink()
            with patch.object(adapter, "fit_chromosome", wraps=adapter.fit_chromosome) as fit:
                adapter.run(self.work, candidate, ncores=1, batch_size=8, resume=True)
            self.assertEqual([call.args[2] for call in fit.call_args_list], [22])
            for name, expected in before.items():
                self.assertEqual((directory / name).read_bytes(), expected)

    def test_resume_fingerprint_ignores_only_orchestration(self):
        source = Path(adapter.__file__).read_text()
        original = self.work / "original.py"
        original.write_text(source)
        expected = adapter.implementation_fingerprint(original)
        original.write_text(source.replace("raise ValueError(\"Verification requires --resume\")",
                                           "raise ValueError(\"Different orchestration message\")"))
        self.assertEqual(expected, adapter.implementation_fingerprint(original))
        original.write_text(source.replace("model.smooth.train(b, y)",
                                           "model.smooth.train(b, y[::-1])"))
        self.assertNotEqual(expected, adapter.implementation_fingerprint(original))


if __name__ == "__main__":
    unittest.main()
