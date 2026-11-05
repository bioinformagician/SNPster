from pathlib import Path
import json
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import validate_phasefree_arrays as validation


class PhasefreeArrayTests(unittest.TestCase):
    def test_pure_selection_is_balanced_and_eligible_only(self):
        rows = [dict(sample_id=f"{source}_{i}", source_group=source)
                for source in validation.SOURCES for i in range(2)]
        eligible = {row["sample_id"] for row in rows}
        selected = validation.pure_controls(rows, eligible)
        self.assertEqual(len(selected), 7)
        np.testing.assert_array_equal(np.sum([row["truth"] for row in selected], axis=0),
                                      np.ones(7))
        self.assertEqual(selected, validation.pure_controls(list(reversed(rows)), eligible))
        with self.assertRaisesRegex(ValueError, "No eligible"):
            validation.pure_controls(rows, eligible - {"OCE_0", "OCE_1"})

    def test_unphased_input_drops_parental_orientation_and_preserves_dosage(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "array.vcf"
            haplotypes = np.array([[[1, 0, 1], [0, 1, 1]]], dtype=np.int8)
            keys = [("22", pos, "A", "G") for pos in (10, 20, 30)]
            validation.write_truth_vcf(path, 22, keys, ["person"], haplotypes, unphased=True)
            rows = [line.split("\t") for line in path.read_text().splitlines()
                    if not line.startswith("#")]
            self.assertEqual([row[-1] for row in rows], ["0/1", "0/1", "1/1"])
            np.testing.assert_array_equal(
                [sum(map(int, row[-1].split("/"))) for row in rows],
                haplotypes.sum(axis=1)[0])
            with self.assertRaisesRegex(ValueError, "axes"):
                validation.write_truth_vcf(path, 22, keys, ["person"], haplotypes[:, :, :2])

    def test_numpy_case_slices_preserve_sample_haplotype_marker_axes(self):
        values = np.arange(4 * 2 * 5).reshape(4, 2, 5)
        selected = values[2:3, :, [0, 3, 4]]
        self.assertEqual(selected.shape, (1, 2, 3))
        np.testing.assert_array_equal(selected[0], values[2][:, [0, 3, 4]])

    def test_checkpoint_reuse_rejects_changed_results_and_experiment(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            data = directory / "dosages.npy"
            np.save(data, np.array([0, 1, 2], dtype=np.int8))
            checkpoint = directory / "completion.json"
            checkpoint.write_text(json.dumps(dict(
                status="complete", experiment_sha256="frozen",
                input_sha256={}, output_sha256={str(data): validation.adapter.digest(data)})))
            self.assertTrue(validation.verified(checkpoint, "frozen"))
            with self.assertRaisesRegex(ValueError, "Incompatible checkpoint"):
                validation.verified(checkpoint, "different")
            np.save(data, np.array([0, 0, 2], dtype=np.int8))
            with self.assertRaisesRegex(ValueError, "provenance changed"):
                validation.verified(checkpoint, "frozen")


if __name__ == "__main__":
    unittest.main()
