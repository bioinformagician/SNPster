import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from make_plans import generate, SOURCES
from prepare_training_pedigrees import validate_pool
from screen_training_decoders import (
    PEDIGREE_WEIGHT, assemble_training, check_disjoint, fit_matched, predict_decoder,
)


class TrainingScreenTests(unittest.TestCase):
    def test_training_pool_cannot_overlap_reference_or_evaluation(self):
        training = {s: (s, s + "_fine") for s in SOURCES}
        labels = dict(training)
        validate_pool(training, {}, labels, set(training), set(), set())
        for reference, evaluation in (({"EUR"}, set()), (set(), {"SEA"})):
            with self.assertRaisesRegex(ValueError, "overlap"):
                validate_pool(training, {}, labels, set(training), reference, evaluation)
        with self.assertRaisesRegex(ValueError, "exactly"):
            validate_pool(training, {}, labels, {"EUR"}, set(), set())

    def test_explicit_cohort_identity_and_donor_checks(self):
        check_disjoint({"train"}, {"train"}, {"eval"}, {"train_node"}, {"eval_node"})
        with self.assertRaisesRegex(ValueError, "founder donors"):
            check_disjoint({"same"}, {"same"}, {"same"}, {"train_node"}, {"eval_node"})
        with self.assertRaisesRegex(ValueError, "IDs collide"):
            check_disjoint({"train"}, {"train"}, {"eval"}, {"same"}, {"same"})

    def test_training_namespace_preserves_default_plans(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            meta = base / "meta.tsv"
            meta.write_text("sample\tfine\tsource_group\n" + "".join(
                f"{s}_{i}\t{s}_fine\t{s}\n" for s in SOURCES for i in range(8)))
            generate(meta, [], base / "default", 9, 10, 1, True)
            generate(meta, [], base / "train", 9, 10, 1, True, family_prefix="trainfam")
            original = pd.read_csv(base / "default/expected_ancestry.tsv", sep="\t")
            training = pd.read_csv(base / "train/expected_ancestry.tsv", sep="\t")
            self.assertFalse(set(original.sample_id) & set(training.sample_id))
            self.assertEqual(original.founder_donor.fillna("").tolist(),
                             training.founder_donor.fillna("").tolist())
            manifest = json.loads((base / "train/manifest.json").read_text())
            self.assertEqual(manifest["family_prefix"], "trainfam")
            with self.assertRaisesRegex(ValueError, "ASCII"):
                generate(meta, [], base / "bad", 9, 10, 1, True, family_prefix="../bad")

    def test_modest_pedigree_weight_and_no_old_mosaic_input(self):
        pure = pd.DataFrame({"sample_id": ["a", "b"], "source_group": ["EUR", "SEA"],
                             "anc_0": [1., 0.], "anc_1": [0., 1.]})
        ped = pd.DataFrame({"sample_id": ["c", "d", "e"], "family_id": ["one", "two", "two"],
                            "anc_0": [.5, .25, .25], "anc_1": [.5, .75, .75]})
        for s in SOURCES:
            pure[s] = (pure.source_group == s).astype(float)
            ped["realized_" + s] = ped.anc_0 if s == "EUR" else ped.anc_1 if s == "SEA" else 0.
        x, y, weights = assemble_training(pure, ped, ["anc_0", "anc_1"], True)
        self.assertEqual(len(x), 5)
        self.assertAlmostEqual(weights[:2].sum(), 1 - PEDIGREE_WEIGHT)
        self.assertAlmostEqual(weights[2:].sum(), PEDIGREE_WEIGHT)
        self.assertAlmostEqual(weights[2], weights[3] + weights[4])
        self.assertEqual(len(assemble_training(pure, ped, ["anc_0", "anc_1"], False)[0]), 2)

    def test_matched_architectures_produce_valid_proportions(self):
        x = np.eye(8)
        y = np.eye(8)
        for kind in ("ridge", "stochastic", "extratrees"):
            model = fit_matched(kind, x, y, np.ones(8) / 8)
            pred = predict_decoder(model, x)
            self.assertTrue(np.isfinite(pred).all())
            self.assertTrue((pred >= -1e-8).all())
            np.testing.assert_allclose(pred.sum(axis=1), 1)


if __name__ == "__main__":
    unittest.main()
