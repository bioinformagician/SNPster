import copy
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

import core
import guided_followup as guided
from make_smoke import fixture
import models


class FounderTruth(unittest.TestCase):
    def setUp(self):
        self.samples = [
            {"sample_id": "d1", "role": "train", "donors": "d1", "source": "EAS", "track": "pure"},
            {"sample_id": "d2", "role": "train", "donors": "d2", "source": "SEA", "track": "pure"},
            {"sample_id": "child", "role": "train", "donors": "d1,d2", "source": "", "track": "mixed_g1"},
        ]
        self.authority = {"d1": "EAS", "d2": "SEA"}
        self.records = []
        for sample in self.samples:
            for hap in ("0", "1"):
                pieces = [(0, 100, sample["sample_id"])] if sample["track"] == "pure" else (
                    [(0, 30, "d1"), (30, 100, "d2")] if hap == "0" else [(0, 100, "d2")])
                for start, end, donor in pieces:
                    self.records.append({"sample_id": sample["sample_id"], "chromosome": "22",
                        "haplotype": hap, "start": str(start), "end": str(end), "source": "EAS_SEA",
                        "donor_id": donor, "donor_haplotype": hap})

    def test_exact_founder_haplotype_recovery_not_random_split_of_merged_truth(self):
        result = guided.expand_inheritance(self.records, self.samples, self.samples, self.authority, 22, [0, 100])
        child = [row for row in result if row["sample_id"] == "child"]
        counts = {source: sum(int(row["end"]) - int(row["start"]) for row in child if row["source"] == source)
                  for source in ("EAS", "SEA")}
        self.assertEqual(counts, {"EAS": 30, "SEA": 170})
        self.assertEqual([{**row, "source": guided.merged(row["source"])} for row in result],
                         [{key: row[key] for key in ("sample_id", "chromosome", "haplotype", "start", "end", "source")}
                          for row in self.records])

    def test_invalid_founder_role_label_haplotype_gap_and_duplicate_fail(self):
        for field, value in (("donor_id", "other"), ("source", "EUR"), ("donor_haplotype", "2"),
                             ("start", "3")):
            changed = copy.deepcopy(self.records)
            changed[-1][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                guided.expand_inheritance(changed, self.samples, self.samples, self.authority, 22, [0, 100])
        with self.assertRaises(ValueError):
            guided.expand_inheritance(self.records + [self.records[-1]], self.samples,
                                      self.samples, self.authority, 22, [0, 100])
        original = copy.deepcopy(self.samples)
        original[1]["role"] = "test"
        with self.assertRaisesRegex(ValueError, "cross-role"):
            guided.expand_inheritance(self.records, self.samples, original, self.authority, 22, [0, 100])

    def test_enrichment_preserves_roles_and_does_not_import_test(self):
        sources = ("EUR", "EAS_SEA", "GME")
        original = [{"sample_id": f"{role}_{index}", "role": role, "donors": f"{role}_{index}",
                     "family_id": f"{role}_{index}", "source": source, "track": "pure", "fine": source}
                    for role in core.ROLES for index, source in enumerate(sources)]
        original.append({**original[1], "sample_id": "calibration_extraSEA", "donors": "calibration_extraSEA",
                         "role": "calibration"})
        authority = {row["sample_id"]: ("EAS" if row["source"] == "EAS_SEA" else row["source"]) for row in original}
        authority["calibration_extraSEA"] = "SEA"
        previous = [row for row in original if row["role"] in ("reference", "train")
                    or (row["role"] in ("calibration", "validation") and row["source"] == "EUR")]
        selected = guided.enriched_samples(original, previous, authority)
        self.assertFalse(any(row["role"] == "test" for row in selected))
        self.assertIn("calibration_extraSEA", {row["sample_id"] for row in selected})
        self.assertEqual(sum(row["role"] == "reference" for row in selected), 3)
        self.assertEqual(sum(row["role"] == "train" for row in selected), 3)
        changed = copy.deepcopy(previous)
        changed[0]["role"] = "test"
        with self.assertRaisesRegex(ValueError, "ownership"):
            guided.enriched_samples(original, changed, authority)


class ReferenceBalance(unittest.TestCase):
    def test_balanced_cap_keeps_rare_sources_and_rotates_fine_populations(self):
        reference = [{"sample_id": f"{source}_{fine}_{i}", "source": source, "fine": fine}
                     for source, fine, size in (("GME", "a", 3), ("GME", "b", 2),
                                                ("EUR", "a", 12), ("EUR", "b", 3),
                                                ("OCE", "a", 2)) for i in range(size)]
        meta = {"reference": reference, "sources": ["EUR", "GME", "OCE"]}
        selected = models.reference_subset(meta, {"reference_source_cap": 5})
        rows = [reference[i] for i in selected]
        self.assertEqual(len([row for row in rows if row["source"] == "GME"]), 5)
        self.assertEqual(len([row for row in rows if row["source"] == "OCE"]), 2)
        self.assertEqual(len([row for row in rows if row["source"] == "EUR"]), 5)
        self.assertEqual(len([row for row in rows if row["source"] == "EUR" and row["fine"] == "b"]), 2)
        self.assertEqual(len(selected), len(set(selected)))
        self.assertEqual(selected, sorted(selected))
        self.assertNotEqual(models.base_policy({"family": "flare", "params": {"reference_source_cap": 5}}),
                            models.base_policy({"family": "flare", "params": {}}))
        for cap in (0, True, 1.5):
            with self.assertRaises(ValueError):
                models.reference_subset(meta, {"reference_source_cap": cap})
        with self.assertRaises(ValueError):
            models.reference_subset(meta, {"reference_source_cap": 5, "reference_cap": 2})


class PhaseChunks(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.panel = self.root / "panel"
        self.panel.mkdir()
        self.samples = [{"sample_id": f"s{i}"} for i in range(5)]
        self.calls = np.array([[[0, 1, 0], [1, 0, 1]]] * 5, dtype=np.int8)
        self.meta = {"queries": self.samples}
        core.write(self.panel / "panel.json", self.meta)
        np.savez(self.panel / "data.npz", query=self.calls)

    def chunk(self, name, start, end):
        output = self.root / name
        output.mkdir()
        np.save(output / "query.npy", self.calls[start:end])
        core.write(output / "chunk.json", {
            "start": start, "end": end, "total": 5, "samples": self.samples[start:end],
            "panel_sha256": core.digest(self.panel / "panel.json"),
            "data_sha256": core.digest(self.panel / "data.npz"),
            "reuse_manifest_sha256": "frozen", "query_sha256": core.digest(output / "query.npy"),
            "artifacts": {}, "reused_singletons": [row["sample_id"] for row in self.samples[start:end]],
        })
        return output

    def test_collection_restores_order_and_exact_genotypes(self):
        folders = [self.chunk("last", 3, 5), self.chunk("first", 0, 3)]
        guided.collect_phase(self.panel, folders, self.root / "phased")
        np.testing.assert_array_equal(np.load(self.root / "phased/query.npy"), self.calls)
        record = core.read(self.root / "phased/phase.json")
        self.assertEqual(record["reused_singletons"], [row["sample_id"] for row in self.samples])
        self.assertEqual(record["changed_diploid_calls"], 0)

    def test_missing_duplicate_changed_and_overlapping_chunks_fail(self):
        a, b = self.chunk("a", 0, 3), self.chunk("b", 3, 5)
        for index, folders in enumerate(([a], [b], [a, a, b], [a, self.chunk("overlap", 2, 5)])):
            with self.subTest(index=index), self.assertRaises(ValueError):
                guided.collect_phase(self.panel, folders, self.root / f"bad{index}")
        np.save(b / "query.npy", self.calls[3:] * 0)
        record = core.read(b / "chunk.json")
        record["query_sha256"] = core.digest(b / "query.npy")
        core.write(b / "chunk.json", record)
        with self.assertRaisesRegex(ValueError, "diploid"):
            guided.collect_phase(self.panel, [a, b], self.root / "changed")

    @unittest.skipUnless(os.environ.get("SNPSTER_TEST_BEAGLE_JAR"), "requires pinned Beagle")
    def test_real_partial_phase_reuse_equals_fresh_singleton_phasing(self):
        tools = {"beagle_jar": os.environ["SNPSTER_TEST_BEAGLE_JAR"], "flare_jar": os.environ["SNPSTER_TEST_FLARE_JAR"]}
        config = fixture(self.root / "inputs", tools,
                         models=[{"id": "flare", "family": "flare", "params": {"em": False, "gen": 2}}])
        config["routes"] = [config["routes"][0]]
        core.write(self.root / "config.json", config)
        core.freeze(self.root / "config.json", self.root / "frozen")
        panel = self.root / "real-panel"
        core.prepare(self.root / "frozen", "array", 22, panel)
        core.phase(panel, self.root / "fresh-phase")
        meta = core.read(panel / "panel.json")
        source_panel = self.root / "source-panel"
        source_panel.mkdir()
        source = self.root / "source-phase"
        source.mkdir()
        old_meta = {**meta, "queries": meta["queries"][:3]}
        core.write(source_panel / "panel.json", old_meta)
        np.save(source / "query.npy", np.load(self.root / "fresh-phase/query.npy")[:3])
        old_phase = core.read(self.root / "fresh-phase/phase.json")
        core.write(source / "phase.json", {**old_phase, "query_sha256": core.digest(source / "query.npy")})
        reuse_path = self.root / "reuse.json"
        core.write(reuse_path, {"seed": config["seed"], "reference_filters": config["reference_filters"],
            "assets": {"array/22": {"panel": str(source_panel), "phased": str(source),
                "hashes": {str(path): core.digest(path) for path in (
                    source_panel / "panel.json", source / "phase.json", source / "query.npy")}}}})
        guided.phase_chunk(panel, reuse_path, 0, 5, self.root / "phase-chunk")
        np.testing.assert_array_equal(np.load(self.root / "phase-chunk/query.npy"),
                                      np.load(self.root / "fresh-phase/query.npy")[:5])
        self.assertEqual(core.read(self.root / "phase-chunk/chunk.json")["reused_singletons"],
                         [row["sample_id"] for row in meta["queries"][:3]])
        (source / "phase.json").write_text("{}\n")
        with self.assertRaisesRegex(ValueError, "reuse evidence changed"):
            guided.phase_chunk(panel, reuse_path, 0, 5, self.root / "changed-phase")


if __name__ == "__main__":
    unittest.main()
