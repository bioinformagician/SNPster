from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import train_matched_global as matched
from screen_v1_decoders import predict_decoder


class MatchedGlobalTests(unittest.TestCase):
    def test_selection_exactly_balances_sources_and_does_not_reuse_donors(self):
        candidates = []
        for i in range(7):
            y = np.zeros(7)
            y[i] = y[(i + 1) % 7] = .5
            candidates.append(dict(sample_id=f"train{i}", family_id=f"family{i}",
                                   ancestors=[f"donor{i}a", f"donor{i}b"], truth=y.tolist()))
        selected = matched.select_training(candidates)
        self.assertEqual(len(selected), 7)
        np.testing.assert_array_equal(np.sum([row["truth"] for row in selected], axis=0), np.ones(7))
        self.assertEqual(len({d for row in selected for d in row["ancestors"]}), 14)
        conflicting = [dict(row) for row in candidates]
        conflicting[0]["ancestors"] = conflicting[1]["ancestors"]
        with self.assertRaisesRegex(ValueError, "infeasible"):
            matched.select_training(conflicting)

    def test_truth_merging_is_individual_level_not_haplotype_specific(self):
        row = {"realized_" + source: 0 for source in ("EUR", "EAS", "SEA", "SAS", "GME",
                                                      "NAM", "OCE", "SSA")}
        row.update(realized_EUR=.5, realized_EAS=.2, realized_SEA=.3)
        np.testing.assert_allclose(matched.truth_vector(row), [.5, .5, 0, 0, 0, 0, 0])

    def training_fixture(self):
        ids = [f"person{i}" for i in range(799)]
        values = np.eye(7)[np.arange(799) % 7]
        metadata = pd.DataFrame(dict(sample_id=ids, family_id=ids,
                                     source_group=np.asarray(matched.adapter.SOURCES)[np.arange(799) % 7]))
        for i, source in enumerate(matched.adapter.SOURCES):
            metadata["truth_" + source] = values[:, i]
        pure, ped = metadata.iloc[:307].copy(), metadata.iloc[307:].copy()
        features = pd.DataFrame(np.tile(values, (1, 2)), columns=matched.FEATURES)
        features["sample_id"] = ids
        imputed = features.iloc[307:314].copy()
        return features, imputed, pure, ped

    def test_fits_two_global_heads_and_rejects_nontraining_individuals(self):
        clean, imputed, pure, ped = self.training_fixture()
        models = matched.fit_global_heads(clean, imputed, pure, ped)
        self.assertEqual(set(models), set(matched.ARMS))
        for model in models.values():
            predictions = predict_decoder(model, clean[matched.FEATURES].to_numpy())
            self.assertEqual(predictions.shape, (799, 7))
            np.testing.assert_allclose(predictions.sum(axis=1), 1, atol=1e-12)
            self.assertTrue(np.all(np.isfinite(predictions)))
        imputed.loc[imputed.index[0], "sample_id"] = "evaluation_person"
        with self.assertRaisesRegex(ValueError, "Only selected training descendants"):
            matched.fit_global_heads(clean, imputed, pure, ped)

    def test_clean_control_does_not_depend_on_imputed_training_features(self):
        clean, imputed, pure, ped = self.training_fixture()
        before = matched.fit_global_heads(clean, imputed, pure, ped)
        imputed[matched.FEATURES] = np.roll(imputed[matched.FEATURES].to_numpy(), 1, axis=1)
        after = matched.fit_global_heads(clean, imputed, pure, ped)
        np.testing.assert_array_equal(before["global_ridge_clean"].coef_,
                                      after["global_ridge_clean"].coef_)
        self.assertGreater(np.max(np.abs(before["global_ridge_matched"].coef_
                                        - after["global_ridge_matched"].coef_)), 1e-6)

    def test_full_genome_summary_fits_and_roundtrips_heads_with_separate_array_track(self):
        clean, imputed, pure, ped = self.training_fixture()
        clean["group"], imputed["group"] = "clean_training", "imputed_training"
        groups = [clean, imputed]
        for name, count in (("pure", 159), ("pedigree", 123), ("real", 1), ("array_imputed", 5)):
            frame = clean.iloc[:count].copy()
            frame["sample_id"] = [f"{name}{i}" for i in range(count)]
            frame["group"] = name
            groups.append(frame)
        frame = pd.concat(groups, ignore_index=True)
        self.assertEqual(len(frame), 1094)
        config = dict(selected_training_ids=imputed.sample_id.tolist(), limitations=[])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            screen, out = root / "screen", root / "out"
            screen.mkdir()
            out.mkdir()
            matched.adapter.write_json(out / "experiment.json", config)
            matched.adapter.write_json(out / "pilot.json", {"input_sha256": {}})
            matched.adapter.write_json(screen / "cohort.json", {})
            original = frame[frame.group.isin(["pure", "pedigree", "real"])][
                ["group", "sample_id", *["smooth_" + s for s in matched.adapter.SOURCES]]].rename(
                    columns={"group": "dataset",
                             **{f"smooth_{s}": f"pred_{s}" for s in matched.adapter.SOURCES}})
            folder = screen / "candidates" / matched.pilot.CANDIDATE
            folder.mkdir(parents=True)
            original.to_csv(folder / "predictions.tsv", sep="\t", index=False)
            for chrom in range(1, 23):
                folder = out / f"features/chr{chrom}"
                folder.mkdir(parents=True)
                part = frame.copy()
                part["chromosome"] = chrom
                part["n_markers"] = 100000 if chrom == 1 else 4500 if chrom < 22 else 7260
                part.to_csv(folder / "features.tsv", sep="\t", index=False)

            def fake_score(directory, arm):
                folder = directory / "candidates" / arm
                predictions = pd.read_csv(folder / "predictions.tsv", sep="\t")
                self.assertEqual(len(predictions), 283)
                self.assertEqual(set(predictions.dataset), {"pure", "pedigree", "real"})
                matched.adapter.write_json(folder / "metrics.json", dict(
                    pure={"mean_tv": 0.}, pedigree={"mean_tv": 0.}, real=[],
                    all_gates_pass=True, gates={}))

            with (patch.object(matched, "frozen_config", return_value=config),
                  patch.object(matched, "checkpoint", return_value=True),
                  patch.object(matched, "training_table", return_value=(pure, ped)),
                  patch.object(matched, "score", side_effect=fake_score)):
                matched.summarize(screen, out)
            result = json.loads((out / "completion.json").read_text())
            self.assertEqual(result["status"], "complete")
            self.assertEqual(set(result["results"]), {"frozen_crf", *matched.ARMS})
            self.assertLess(result["max_original_clean_validation_error"], 1e-12)
            for arm in matched.ARMS:
                self.assertEqual(result["results"][arm]["max_decoder_roundtrip_error"], 0.)
                self.assertTrue((out / f"candidates/{arm}/global_decoder.pkl").is_file())


if __name__ == "__main__":
    unittest.main()
