import copy
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))

from audit_production_ancestry import SOURCES, digest, write_json
import dense_model_expansion as expansion
from production_matched_preprocessing import write_vcf


class ExpandedModelTests(unittest.TestCase):
    def test_98_unique_effective_configurations_include_unchanged_baselines(self):
        baseline = expansion.read(PROJECT / "config/production_matched_protocol.json")["candidates"]
        plan = expansion.read(PROJECT / "config/dense_model_expansion.json")
        candidates = expansion.candidate_grid(baseline, plan)
        self.assertEqual(len(candidates), 98)
        self.assertEqual(len({expansion.signature(policy) for policy in candidates.values()}), 98)
        self.assertEqual(sum(policy["family"] == "flare" for policy in candidates.values()), 64)
        self.assertEqual(sum(policy["family"] == "gnomix" for policy in candidates.values()), 32)
        self.assertEqual(sum(policy["family"] == "phasefree" for policy in candidates.values()), 2)
        for name, policy in baseline.items():
            if name in candidates:
                self.assertEqual(candidates[name], policy)
            else:
                self.assertIn(expansion.signature(policy), {expansion.signature(row) for row in candidates.values()})
        self.assertEqual(len(set(baseline).intersection(candidates)), 8)
        self.assertEqual(expansion.signature(dict(family="flare", reference="full", em=False, array=True)),
                         expansion.signature(dict(family="flare", reference="full", em=False, array=False, gen=10)))
        plan["maximum_configurations"] = 97
        with self.assertRaisesRegex(ValueError, "count contract"):
            expansion.candidate_grid(baseline, plan)

    def test_flare_grid_settings_reach_exact_command(self):
        policy = dict(family="flare", reference="finecap8", array=False, em=True, gen=50)
        command = list(map(str, expansion.flare_command(
            Path("/panel"), dict(map="/map"), policy, "lonni", Path("/model/lonni"))))
        for value in ("gen=50", "array=false", "em=true", "update-p=false",
                      "ref=/panel/finecap8.vcf.gz", "ref-panel=/panel/finecap8.panel", "seed=20261022"):
            self.assertIn(value, command)
        self.assertFalse(any(value.startswith("model=") for value in command))
        policy["reference"] = "full"
        self.assertIn("ref=/panel/reference.vcf.gz",
                      expansion.flare_command(Path("/panel"), dict(map="/map"), policy, "lonni", Path("/out")))

    def test_expansion_cannot_change_donors_truth_or_production(self):
        original = dict(protocol=dict(candidates={}, production=dict(unchanged=True)),
                        people=[dict(sample_id="heldout")], truth_by_chromosome={"22": [1, 0]})
        expanded = copy.deepcopy(original)
        expanded["candidate_expansion"] = {}
        expanded["protocol"]["candidates"] = {"new": {}}
        expansion.verify_same_population(original, expanded)
        expanded["truth_by_chromosome"]["22"][0] = .5
        with self.assertRaisesRegex(ValueError, "changed"):
            expansion.verify_same_population(original, expanded)

    def test_fine_caps_are_deterministic_without_prediction_based_selection(self):
        records = [dict(sample_id=identity, fine=fine) for identity, fine in
                   (("c", "x"), ("a", "x"), ("b", "x"), ("e", "y"), ("d", "y"))]
        self.assertEqual([row["sample_id"] for row in expansion.fine_subset(records, 2)], ["a", "b", "d", "e"])
        self.assertEqual(expansion.fine_subset(records, 2), expansion.fine_subset(records[::-1], 2))

    def test_partial_source_scope_is_rejected_before_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scope = root / "full"
            (scope / "summary/summary").mkdir(parents=True)
            original = dict(protocol=dict(candidates={}, masks=["lonni", "other"]))
            original_path = root / "original.json"
            write_json(original_path, original)
            cohort = copy.deepcopy(original)
            cohort["candidate_expansion"] = dict(source_cohort=str(original_path),
                source_cohort_sha256=digest(original_path), baseline_aliases={})
            cohort_path = root / "expanded.json"
            write_json(cohort_path, cohort)
            write_json(scope / "published_verification.json", dict(published_sha256={}, source_sha256={}))
            write_json(scope / "post_cleanup_verification.json",
                       dict(manifest_sha256=digest(scope / "published_verification.json")))
            write_json(scope / "summary/summary/independent_verification.json",
                       dict(status="complete_verified_development_comparison", model_configurations=10,
                            cohort_sha256=digest(original_path), chromosomes=list(range(1, 22)),
                            masks=["lonni", "other"], n_evaluation_people=15))
            with self.assertRaisesRegex(ValueError, "complete, independently verified"):
                expansion.verify_source(cohort_path, scope, root / "gate.json")
            self.assertFalse((root / "gate.json").exists())

    def test_capped_vcfs_preserve_every_allele_and_all_eight_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scope = root / "original"
            panel = scope / "panels/chr22/panel"
            panel.mkdir(parents=True)
            refs = [dict(sample_id=f"ref{source}_{i:02d}", source=source, fine=f"fine_{source}")
                    for source in SOURCES for i in range(40)]
            ids = [row["sample_id"] for row in refs]
            keys = [("22", 1000 + 10 * i, "A", "C") for i in range(120)]
            calls = np.random.default_rng(1).integers(0, 2, (len(ids), 2, len(keys)), dtype=np.int8)
            write_vcf(panel / "reference.vcf", keys, ids, calls, phased=True)
            (panel / "reference.ids").write_text("\n".join(ids) + "\n")
            np.savez_compressed(panel / "reference.npz", haplotypes=calls.reshape(2 * len(ids), -1))
            write_json(panel / "panel.json", dict(chromosome=22, keys=keys, finecap_people=128))
            cohort_path = root / "cohort.json"
            write_json(cohort_path, dict(reference=refs))
            gate = root / "source_verified.json"
            write_json(gate, dict(status="verified_complete_original_scope_reusable",
                                 cohort_sha256=digest(cohort_path), source_scope=str(scope)))
            expansion.panels(cohort_path, panel, gate, root / "expanded")
            metadata = expansion.read(root / "expanded/panel.json")
            self.assertEqual(metadata["expanded_reference_people"],
                             dict(full=320, finecap=128, finecap8=64, finecap32=256))
            self.assertTrue((root / "expanded/source_gate.json").is_file())

    def test_baseline_reuse_records_original_provenance_and_rejects_policy_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            policy = dict(family="phasefree", reference="full", input="raw_harmonized")
            original_path = root / "original_cohort.json"
            write_json(original_path, dict(protocol=dict(candidates={"baseline": policy})))
            source_scope = root / "original_scope"
            model = source_scope / "models/baseline/chr22/model"
            model.mkdir(parents=True)
            (model / "predictions.tsv").write_text("unchanged scientific predictions\n")
            write_json(model / "complete.json", dict(cohort_sha256=digest(original_path),
                       panel_sha256="original_panel", predictions_sha256=digest(model / "predictions.tsv")))
            cohort_path = root / "expanded_cohort.json"
            cohort = dict(protocol=dict(candidates={"baseline": policy}),
                          candidate_expansion=dict(source_cohort=str(original_path),
                                                   source_cohort_sha256=digest(original_path)))
            write_json(cohort_path, cohort)
            panel = root / "panel"
            panel.mkdir()
            write_json(panel / "panel.json", dict(chromosome=22, source_panel_sha256="original_panel"))
            write_json(panel / "source_gate.json", dict(source_scope=str(source_scope)))
            expansion.adopt_baseline(cohort_path, panel, "baseline", root / "adopted")
            record = expansion.read(root / "adopted/complete.json")
            self.assertEqual(record["original_cohort_sha256"], digest(original_path))
            self.assertEqual(record["execution_mode"], "reused_independently_verified_identical_baseline")
            self.assertEqual(digest(root / "adopted/predictions.tsv"), digest(model / "predictions.tsv"))
            cohort["protocol"]["candidates"]["baseline"]["input"] = "different"
            write_json(cohort_path, cohort)
            with self.assertRaisesRegex(ValueError, "candidate policy"):
                expansion.adopt_baseline(cohort_path, panel, "baseline", root / "wrong")

    def test_capacity_gate_measures_only_new_assets_and_rejects_insufficient_disk(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            smoke = root / "smoke"
            (smoke / "summary/summary").mkdir(parents=True)
            write_json(smoke / "summary/summary/independent_verification.json",
                       dict(status="smoke_passed_not_a_whole_genome_leaderboard", model_configurations=98))
            write_json(smoke / "published_verification.json", dict(published_sha256={}, source_sha256={}))
            write_json(smoke / "post_cleanup_verification.json",
                       dict(manifest_sha256=digest(smoke / "published_verification.json")))
            for name in ("panels", "models", "reference_bases"):
                (smoke / name).mkdir()
                (smoke / name / "asset").write_bytes(b"12345")
            for chromosome in range(1, 23):
                folder = root / f"full/panels/chr{chromosome}/panel"
                folder.mkdir(parents=True)
                write_json(folder / "panel.json", dict(marker_count=100))
            cohort = root / "cohort.json"
            write_json(cohort, dict(chromosome_spans={str(i): [1, 101] for i in range(1, 23)}))
            with patch.object(expansion.shutil, "disk_usage", return_value=SimpleNamespace(free=1)):
                with self.assertRaisesRegex(ValueError, "capacity gate failed"):
                    expansion.capacity(cohort, root / "full", smoke, root / "capacity.json")
            result = expansion.read(root / "capacity.json")
            self.assertEqual(result["measured_smoke_model_panel_and_base_bytes"], 15)
            self.assertEqual(result["scaling_ratio"], 22)
            self.assertEqual(result["required_bytes"], 4 * 15 * 22 + 20 * 1024**3)

    def test_pinned_flare_uses_nondefault_generation_parameters(self):
        import gzip
        from screen_model_families import read_model
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ids = [f"ref{source}_{i}" for source in SOURCES for i in range(16)]
            keys = [("22", 100000 + 10000 * i, "A", "C") for i in range(120)]
            calls = np.random.default_rng(20261022).integers(0, 2, (len(ids), 2, len(keys)), dtype=np.int8)
            write_vcf(root / "reference.vcf", keys, ids, calls, phased=True)
            write_vcf(root / "lonni.evaluation.vcf", keys, ["target_a", "target_b"], calls[:2], phased=True)
            (root / "full.panel").write_text("".join(f"{identity}\t{source}\n" for source in SOURCES
                                                   for identity in ids if identity.startswith(f"ref{source}_")))
            genetic_map = root / "map.txt"
            genetic_map.write_text(f"22 . 0 {keys[0][1]}\n22 . 12 {keys[-1][1]}\n")
            local = {}
            for generation, array in ((5, True), (50, True), (50, False)):
                policy = dict(family="flare", reference="full", array=array, em=False, gen=generation)
                prefix = root / f"gen{generation}_{array}"
                expansion.run(expansion.flare_command(root, dict(map=str(genetic_map)),
                              policy, "lonni", prefix), root / f"run{generation}_{array}")
                self.assertEqual(read_model(Path(str(prefix) + ".model"))["generation"], generation)
                with gzip.open(str(prefix) + ".anc.vcf.gz", "rt") as handle:
                    local[generation, array] = [line for line in handle if not line.startswith("#")]
            self.assertTrue(local[50, True])
            self.assertEqual(local[50, True], local[50, False])


if __name__ == "__main__":
    unittest.main()
