import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from audit_production_ancestry import digest
from recover_sgdp_duplicates import collect, duplicates, repair, select_record
from verify_sgdp_full_manifest import require_hashes


class DuplicateRecoveryTests(unittest.TestCase):
    def record(self, value, position=100):
        return f"8\t{position}\t.\tG\tA\t.\t.\t.\tGT:GQ:DP:PL:AD\t{value}\n"

    def test_supported_record_wins_in_either_order_without_changing_evidence(self):
        supported = self.record("0/0:120:40:0,120,255:40,0")
        empty = self.record("./.:.:.:.:.")
        self.assertEqual(select_record([supported, empty]), supported)
        self.assertEqual(select_record([empty, supported]), supported)
        self.assertEqual(select_record([empty, supported, supported]), supported)

    def test_identical_and_wholly_missing_duplicates_collapse(self):
        for value in ("0/1:60:40:60,0,60:20,20", "./.:.:.:.:."):
            line = self.record(value)
            self.assertEqual(select_record([line, line]), line)

    def test_conflicting_genotypes_or_raw_evidence_fail(self):
        supported = self.record("0/0:120:40:0,120,255:40,0")
        for value in ("0/1:60:40:60,0,60:20,20", "0/0:123:41:0,123,255:41,0",
                      "./.:1:40:0,1,255:40,0"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "Conflicting"):
                select_record([supported, self.record(value)])

    def test_real_vcf_repair_retains_every_other_record_and_index_matches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "original"
            smoke = Path(directory) / "smoke"
            output = Path(directory) / "recovery"
            for folder in (root / "basecalls", root / "source", smoke / "metadata",
                           smoke / "raw_control", output):
                folder.mkdir(parents=True)
            header = (
                "##fileformat=VCFv4.2\n##contig=<ID=8>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson\n"
            )
            supported = self.record("0/0:120:40:0,120,255:40,0")
            empty = self.record("./.:.:.:.:.")
            other = self.record("0/1:60:40:60,0,60:20,20", 200)
            import subprocess
            source = root / "basecalls/basecalls.vcf"
            source.write_text(header + empty + supported + other)
            subprocess.run(["bgzip", str(source)], check=True)
            calls = root / "basecalls/basecalls.vcf.gz"
            subprocess.run(["bcftools", "index", str(calls)], check=True)
            sites, repeated = duplicates(calls, 8)
            self.assertEqual(sites, 3)
            self.assertEqual(sum(repeated.values()), 1)
            (smoke / "metadata/selected.json").write_text(json.dumps({"people": [{"sample_id": "person"}]}))
            (root / "basecalls/basecalls.json").write_text(json.dumps({
                "sites": 3, "output_sha256": digest(calls), "called": {"person": 2},
            }))
            (root / "basecalls/likelihood_query.command.json").write_text("[]")
            (root / "source/download.json").write_text("{}")
            (smoke / "raw_control/raw_control.vcf.gz").write_bytes(b"raw")
            real_run = subprocess.run

            def run(command, **kwargs):
                if command[0] == "python3":
                    control = output / "raw_control.vcf.gz.csi"
                    control.write_bytes(b"index")
                    return subprocess.CompletedProcess(command, 0)
                return real_run(command, **kwargs)

            with patch("recover_sgdp_duplicates.subprocess.run", side_effect=run):
                normalized, evidence = repair(root, smoke, 8, repeated, output)
            with gzip.open(normalized / "basecalls/basecalls.vcf.gz", "rt") as handle:
                self.assertEqual(handle.read(), header + supported + other)
            result = json.loads((normalized / "basecalls/basecalls.json").read_text())
            self.assertEqual(result["sites"], 2)
            self.assertEqual(result["called"], {"person": 2})
            self.assertEqual(evidence[0]["removed_records"], 1)
            self.assertEqual(subprocess.check_output([
                "bcftools", "index", "-n", str(normalized / "basecalls/basecalls.vcf.gz"),
            ], text=True).strip(), "2")

    def test_collector_binds_repaired_staged_files_to_durable_published_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            import_run, smoke, previous, recovery_run, staging = (
                base / name for name in ("import", "smoke", "previous", "recovery", "staging"))
            for run in (import_run, previous, recovery_run):
                run.mkdir()
                (run / "source_sha256.txt").write_text("")
            (smoke / "metadata").mkdir(parents=True)
            (smoke / "metadata/selected.json").write_text(json.dumps({"people": [{"sample_id": "person"}]}))
            (import_run / "summary").mkdir()
            (import_run / "summary/full_import.json").write_text(json.dumps({
                "raw_call_qc": [{"chromosome": chromosome, "concordance": 0.999, "n_compared": 1000}
                               for chromosome in range(1, 23)],
            }))
            folders = []
            for chromosome in range(1, 23):
                folder = staging / f"chr{chromosome}"
                folder.mkdir(parents=True)
                folders.append(folder)
                original = smoke if chromosome == 22 else import_run / f"chromosomes/chr{chromosome}"
                record = dict(
                    chromosome=chromosome, sites=1, genotype_cells_independently_rederived=1,
                    sample_ids=["person"], original_raw_control_concordance=0.999,
                    original_raw_control_compared=1000, sha256={}, shared_sha256={},
                    status="published_genotypes_verified_not_ancestry_execution_ready",
                    scientific_ancestry_launch_authorized=False,
                )
                details = dict(chromosome=chromosome, original_root=str(original),
                               original_sites=1, removed_records=0, status="no_duplicates")
                if chromosome == 8:
                    path = folder / "normalized/basecalls/basecalls.vcf.gz"
                    path.parent.mkdir(parents=True)
                    path.write_bytes(b"retained evidence")
                    record["sha256"] = {"basecalls/basecalls.vcf.gz": digest(path)}
                    details.update(original_sites=2, removed_records=1,
                                   status="corrected_redundant_duplicates")
                (folder / "independent.json").write_text(json.dumps(record))
                (folder / "recovery.json").write_text(json.dumps(details))
            linked = base / "task"
            linked.mkdir()
            for folder in folders:
                (linked / folder.name).symlink_to(folder, target_is_directory=True)
            folders = [linked / folder.name for folder in folders]
            manifest = recovery_run / "published_verification.json"
            with self.assertRaisesRegex(ValueError, "per autosome"):
                collect(import_run, smoke, previous, recovery_run, folders[:-1], manifest)
            collect(import_run, smoke, previous, recovery_run, folders, manifest)
            result = json.loads(manifest.read_text())
            self.assertEqual(result["removed_redundant_records"], 1)
            self.assertEqual(result["data_roots"]["8"],
                             str(recovery_run / "chromosomes/chr8/normalized"))
            self.assertFalse(result["scientific_ancestry_launch_authorized"])
            self.assertFalse(any(str(staging) in path for path in result["inputs_sha256"]))
            import shutil
            shutil.copytree(staging, recovery_run / "chromosomes")
            require_hashes(result["inputs_sha256"])
            (recovery_run / "chromosomes/chr8/normalized/basecalls/basecalls.vcf.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Durable evidence changed"):
                require_hashes(result["inputs_sha256"])


if __name__ == "__main__":
    unittest.main()
