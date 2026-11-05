import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from audit_production_ancestry import digest
from sgdp_genetic_qc import collect, kinship
from verify_sgdp_full_manifest import require_hashes


class GeneticQCTests(unittest.TestCase):
    def test_all22_actual_concat_pruning_and_king_include_every_candidate_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "config").mkdir()
            protocol = json.loads((PROJECT / "config/sgdp_normalization.json").read_text())
            policy = protocol["relatedness"]
            policy.update(minimum_pruned_markers=20, minimum_pair_called_markers=20)
            (root / "config/sgdp_normalization.json").write_text(json.dumps(protocol))
            ids = [f"panel{i}" for i in range(60)] + ["candidate1", "candidate2"]
            randomizer, panels = random.Random(20261022), []
            for chromosome in range(1, 23):
                folder = root / f"chr{chromosome}"
                folder.mkdir()
                panels.append(folder)
                path = folder / "panel.vcf"
                with path.open("w") as handle:
                    handle.write(f"##fileformat=VCFv4.2\n##contig=<ID={chromosome}>\n"
                                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                                 "\t".join(ids) + "\n")
                    for position in range(1, 31):
                        calls = [randomizer.choice(("0/0", "0/1", "1/1")) for _ in ids]
                        calls[-2] = calls[0]
                        handle.write(f"{chromosome}\t{position * 1000000}\t{chromosome}:{position}:A:C"
                                     "\tA\tC\t.\t.\t.\tGT\t" + "\t".join(calls) + "\n")
                subprocess.run(["bgzip", str(path)], check=True)
                subprocess.run(["bcftools", "index", str(path) + ".gz"], check=True)
                (folder / "panel.json").write_text(json.dumps(dict(
                    chromosome=chromosome, sites=30, sample_ids=ids, candidate_ids=ids[-2:],
                    called={identity: 30 for identity in ids},
                )))
            kinship(root, panels, root / "king")
            record = json.loads((root / "king/genetic_eligibility.json").read_text())
            self.assertEqual(record["pairs_checked"], 121)
            self.assertIn("candidate1", record["excluded_candidate_ids"])
            self.assertGreaterEqual(record["pruned_markers"], 20)
            self.assertFalse(record["scientific_ancestry_launch_authorized"])

    def test_normalized_collector_uses_only_durable_paths_and_detects_changed_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            project, run_dir, staged = root / "project", root / "published", root / "task"
            (project / "config").mkdir(parents=True)
            (project / "config/sgdp_normalization.json").write_text(json.dumps({"remaining_gates": ["cohort"]}))
            run_dir.mkdir()
            (run_dir / "source_sha256.txt").write_text("")
            folders, reports, panels = [], [], []
            for chromosome in range(1, 23):
                folder = staged / f"normalized/chr{chromosome}"
                folder.mkdir(parents=True)
                folders.append(folder)
                calls = folder / "normalized.vcf.gz"
                calls.write_bytes(b"calls")
                (folder / "normalization.json").write_text(json.dumps({"chromosome": chromosome}))
                report = root / f"chr{chromosome}.verified.json"
                report.write_text(json.dumps(dict(
                    chromosome=chromosome, sites=1, genotype_cells_verified=2, sample_ids=["a", "b"],
                    normalized_sha256=digest(calls),
                )))
                reports.append(report)
                panel = staged / f"panels/chr{chromosome}"
                panel.mkdir(parents=True)
                (panel / "panel.json").write_text(json.dumps({"chromosome": chromosome}))
                panels.append(panel)
            king = staged / "kinship"
            king.mkdir()
            (king / "genetic_eligibility.json").write_text(json.dumps({
                "pairs_checked": 1, "scientific_ancestry_launch_authorized": False,
            }))
            output = root / "normalization_complete.json"
            with self.assertRaisesRegex(ValueError, "per autosome"):
                collect(project, run_dir, folders, reports[:-1], panels, king, output)
            collect(project, run_dir, folders, reports, panels, king, output)
            record = json.loads(output.read_text())
            self.assertEqual(record["genotype_cells_verified"], 44)
            self.assertFalse(record["scientific_ancestry_launch_authorized"])
            self.assertFalse(any(str(staged) in path for path in record["published_sha256"]))
            shutil.copytree(staged / "normalized", run_dir / "normalized")
            shutil.copytree(staged / "panels", run_dir / "identity_panels")
            shutil.copytree(king, run_dir / "kinship")
            (run_dir / "verification").mkdir()
            for report in reports:
                shutil.copyfile(report, run_dir / "verification" / report.name)
            require_hashes(record["published_sha256"])
            (run_dir / "normalized/chr1/normalized.vcf.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Durable evidence changed"):
                require_hashes(record["published_sha256"])


if __name__ == "__main__":
    unittest.main()
