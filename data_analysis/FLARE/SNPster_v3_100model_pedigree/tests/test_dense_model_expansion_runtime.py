from pathlib import Path
import pickle
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from audit_production_ancestry import write_json
import dense_model_expansion as expansion


class ExpandedGnomixRuntimeTests(unittest.TestCase):
    def test_shared_reference_base_is_unchanged_and_l2_reaches_real_crf(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel = root / "panel"
            panel.mkdir()
            n_markers = 200
            rng = np.random.default_rng(20261022)
            labels = np.repeat(np.arange(8), 4)
            haplotypes = rng.integers(0, 2, (2 * len(labels), n_markers), dtype=np.int8)
            keys = [("22", 100000 + 5000 * i, "A", "C") for i in range(n_markers)]
            genetic_map = root / "map.txt"
            genetic_map.write_text(f"22 . 0 {keys[0][1]}\n22 . 12 {keys[-1][1]}\n")
            metadata = dict(keys=keys, map=str(genetic_map), chromosome=22,
                            span=[keys[0][1], keys[-1][1] + 1], evaluation=["a", "b"])
            write_json(panel / "panel.json", metadata)
            np.savez_compressed(panel / "reference.npz", labels=labels, haplotypes=haplotypes)
            training_labels = np.repeat(np.repeat(np.arange(8), 2)[:, None], n_markers, axis=1)
            np.savez_compressed(panel / "lonni.training.npz",
                                haplotypes=haplotypes[np.repeat(np.arange(8) * 8, 2) +
                                                     np.tile(np.arange(2), 8)], labels=training_labels)
            np.savez_compressed(panel / "lonni.evaluation.npz", haplotypes=haplotypes[:4])
            expansion.fit_base(panel, .5, root / "base")
            before = expansion.digest(root / "base/base.pkl")
            values = []
            for c2 in (.01, 10):
                output = root / str(c2)
                output.mkdir()
                result = expansion.fit_smoother(panel, metadata,
                    dict(window_size_cM=.5, crf_c2=c2), "lonni", root / "base", output)
                with (output / "lonni.model.pkl").open("rb") as handle:
                    model = pickle.load(handle)["model"]
                self.assertEqual(model.smooth.model.CRF.get_params()["c2"], c2)
                self.assertEqual(result.shape, (2, 8))
                np.testing.assert_allclose(result.sum(axis=1), 1, rtol=0, atol=1e-10)
                values.append(result)
            self.assertEqual(before, expansion.digest(root / "base/base.pkl"))
            self.assertGreater(float(np.max(np.abs(values[0] - values[1]))), 1e-6)


if __name__ == "__main__":
    unittest.main()
