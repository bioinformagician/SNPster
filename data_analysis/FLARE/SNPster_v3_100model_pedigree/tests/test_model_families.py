from pathlib import Path
import gzip
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from screen_model_families import (
    MODES, SOURCES, ancestry_priors, canonical_genotypes, collapse, component_source, copying_profiles,
    empirical_centers, midpoint_weights, physical_predictions, read_model, unified_seven, write_model,
)
from screen_v1_decoders import SOURCES as ORIGINAL_SOURCES


class ModelFamilyTests(unittest.TestCase):
    def test_exact_twelve_configs_and_only_asian_merge(self):
        self.assertEqual(len(MODES), 12)
        self.assertEqual(len(set(MODES)), 12)
        original = np.eye(len(ORIGINAL_SOURCES))
        result = collapse(original)
        self.assertEqual(result.shape, (8, 7))
        np.testing.assert_array_equal(result[1], result[3])
        self.assertEqual(len(np.unique(result, axis=0)), 7)
        with self.assertRaises(ValueError):
            collapse(np.ones((1, 8)))

    def test_physical_weights_cover_span_and_reject_bad_positions(self):
        np.testing.assert_allclose(midpoint_weights([10, 20, 50], [0, 60]), [15, 20, 25])
        np.testing.assert_allclose(midpoint_weights([20], [0, 60]), [60])
        for positions, span in (([10, 10], [0, 60]), ([20, 10], [0, 60]),
                                ([10, 20], [15, 60]), ([10, 20], [0, 20])):
            with self.assertRaises(ValueError):
                midpoint_weights(positions, span)

    def test_source_components_and_training_priors_cannot_cross_source(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            frame = pd.DataFrame(np.eye(7), columns=SOURCES)
            frame.insert(0, "sample_id", list(SOURCES))
            frame.to_csv(path / "truth.tsv", sep="\t", index=False)
            names = [f"{s}__{i}" for s in SOURCES for i in (1, 2)]
            ancestry_priors(path / "truth.tsv", {"ancestries": names}, path / "priors.tsv")
            result = pd.read_csv(path / "priors.tsv", sep="\t").set_index("SAMPLE")
            for source in SOURCES:
                self.assertEqual(result.loc[source, f"{source}__1"], .5)
                self.assertEqual(result.loc[source, f"{source}__2"], .5)
                self.assertEqual(result.loc[source].sum(), 1.)
                for name in names:
                    if component_source(name) != source:
                        self.assertEqual(result.loc[source, name], 0.)
            with self.assertRaises(ValueError):
                component_source("unanchored")

    def test_physical_posterior_aggregation_and_exact_marker_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "local.vcf.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson\n")
                stream.write("22\t10\t.\tA\tC\t.\t.\t.\tANP1:ANP2\t1,0:1,0\n")
                stream.write("22\t20\t.\tA\tC\t.\t.\t.\tANP1:ANP2\t0,1:0,1\n")
            markers = pd.DataFrame({"chrom": [22, 22], "pos": [10, 20],
                                    "ref": ["A", "A"], "alt": ["C", "C"]})
            ids, prediction = physical_predictions(
                path, ["EUR", "SSA"], {"22": [0, 30]}, markers)
            self.assertEqual(ids, ["person"])
            self.assertEqual(prediction[0, SOURCES.index("EUR")], .5)
            self.assertEqual(prediction[0, SOURCES.index("SSA")], .5)
            markers.loc[1, "alt"] = "G"
            with self.assertRaises(ValueError):
                physical_predictions(path, ["EUR", "SSA"], {"22": [0, 30]}, markers)

    def test_model_roundtrip_uniform_priors_and_fixed_generation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "model"
            model = {"ancestries": list(SOURCES), "panels": ["a", "b"],
                     "p": np.tile([.1, .9], (7, 1)),
                     "theta": np.full((7, 2), .001), "rho": np.ones(7)}
            write_model(path, model)
            restored = read_model(path)
            np.testing.assert_array_equal(restored["p"], model["p"])
            np.testing.assert_allclose(restored["mu"], np.full(7, 1 / 7))
            self.assertEqual(restored["generation"], 10.)
            model["p"][0] = [-1e-18, 1.]
            with self.assertRaises(ValueError):
                write_model(path, model)

    def test_empirical_cluster_centers_preserve_exact_zero_coordinates(self):
        values = np.array([[0., .1, .9, 2.], [0., .3, .7, 4.],
                           [0., .8, .2, 6.], [0., 1., 0., 8.]])
        labels = np.array([0, 0, 1, 1])
        centers = empirical_centers(values, labels, 2)
        np.testing.assert_allclose(centers, [[0., .2, .8, 3.], [0., .9, .1, 7.]])
        self.assertTrue((centers >= 0).all())
        self.assertTrue((centers[:, 0] == 0).all())
        with self.assertRaises(ValueError):
            empirical_centers(values, np.zeros(4, int), 2)

    def test_canonicalization_keeps_matching_second_alt_without_changing_phase(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            header = ("##fileformat=VCFv4.2\n##contig=<ID=22,length=100>\n"
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
            for name, records in (
                ("ref", "22\t10\t.\tC\tT\t.\t.\t.\tGT\t0|1\n"
                        "22\t20\t.\tA\tC\t.\t.\t.\tGT\t0|1\n"),
                ("query", "22\t10\t.\tC\tG\t.\t.\t.\tGT\t0|0\n"
                          "22\t10\t.\tC\tT\t.\t.\t.\tGT\t0|1\n"
                          "22\t20\t.\tA\tC\t.\t.\t.\tGT\t1|1\n"),
            ):
                (root / f"{name}.vcf").write_text(header + records)
                subprocess.run(["bcftools", "view", "-Oz", "-o", str(root / f"{name}.vcf.gz"),
                                str(root / f"{name}.vcf")], check=True)
                subprocess.run(["bcftools", "index", "-t", str(root / f"{name}.vcf.gz")],
                                check=True)
            markers = pd.DataFrame({"chrom": [22, 22], "pos": [10, 20],
                                    "ref": ["C", "A"], "alt": ["T", "C"]})
            result, audit = canonical_genotypes(
                root / "query.vcf.gz", root / "ref.vcf.gz", markers, root, "test")
            genotypes = subprocess.run(["bcftools", "query", "-f", "[%GT\\n]", result],
                                       check=True, capture_output=True, text=True).stdout.splitlines()
            self.assertEqual(genotypes, ["0|1", "1|1"])
            self.assertEqual(audit["original_duplicate_position_records"], 1)
            self.assertEqual(audit["canonical_records"], 2)

    def test_copying_profiles_verify_full_haplotype_order_and_source_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "profiles"
            header = "##nRefHaps=20\n#WIND COORD GT_HAP SWITCH p1 p2\n"
            rows = [f"0 22:10-20 {i} {2 if i % 2 == 0 else 4} "
                    f"{'.2 .8' if i % 2 == 0 else '.4 .6'}\n" for i in range(14)]
            path.write_text(header + "".join(rows))
            result = copying_profiles(path, list(SOURCES), dict(zip(SOURCES, SOURCES)), 1)
            np.testing.assert_allclose(result[2][0], [.3, .7, 3.])
            path.write_text(header + rows[0] + rows[0])
            with self.assertRaises(ValueError):
                copying_profiles(path, list(SOURCES), dict(zip(SOURCES, SOURCES)), 1)

    def test_balancing_families_and_combining_minor_asian_components(self):
        cfg = {"absent_ancestry_max": .01, "minor_min": .02, "minor_max": .05,
               "minor_detection_threshold": .01}
        frame = pd.DataFrame({"sample_id": ["a", "b", "c"],
                              "family_id": ["f1", "f2", "f2"],
                              "scenario": ["f1"] * 3})
        for source in SOURCES:
            frame[f"truth_{source}"] = 1. if source == "EUR" else 0.
            frame[f"pred_{source}"] = ([0., 1., 1.] if source == "EUR"
                                      else [1., 0., 0.] if source == "SSA"
                                      else [0., 0., 0.])
        result = unified_seven(frame, cfg)
        self.assertEqual(result["balanced_mean_tv"], .5)
        self.assertEqual(result["balanced_p95_tv"], 1.)
        self.assertEqual(result["minor_2_5_n"], 0)
        self.assertIsNone(result["minor_2_5_recall"])


if __name__ == "__main__":
    unittest.main()
