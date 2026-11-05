import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from audit_production_ancestry import digest
from verify_public_sgdp import verify
from verify_sgdp_full_manifest import collect, recheck


class PublicGenotypeVerificationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.project = Path(self.directory.name)
        self.smoke = self.project / "smoke"
        self.run = self.project / "import"
        self.checks = self.project / "independent"
        for root, folders in ((self.project, ("config", "results_v3_production_preflight_20261021_193442/summary")),
                              (self.smoke, ("metadata", "raw_control")),
                              (self.run, ("basecalls", "source", "verification", "summary")),
                              (self.checks, ("summary",))):
            for folder in folders:
                (root / folder).mkdir(parents=True)
        self.protocol = json.loads((PROJECT / "config/sgdp_public_import.json").read_text())
        self.protocol["expected_metadata_disjoint_counts"] = {"NAM": 1}
        self.write(self.project / "config/sgdp_public_import.json", self.protocol)
        preflight = self.project / "results_v3_production_preflight_20261021_193442/summary/preflight.json"
        self.write(preflight, {"references": [{"sample_ids": ["forbidden"]}], "cohorts": {}})
        metadata = self.smoke / "metadata/metadata.original.txt"
        metadata.write_text("SGDP_ID\tIllumina_ID\tSample_ID\tSample_ID(Aliases)\tPopulation_ID\tEmbargo\tcategory\n"
                            "S_Chane-1\tDNA1\tperson\talias\tChane\tFullyPublic\tX\n")
        self.write(self.smoke / "metadata/selected.json", {
            "source_metadata_sha256": digest(metadata), "production_preflight_sha256": digest(preflight),
            "people": [{"sample_id": "S_Chane-1.DG", "source": "NAM",
                        "aliases": sorted(["S_Chane-1", "DNA1", "person", "alias"])}],
        })
        self.calls = self.run / "basecalls/basecalls.vcf.gz"
        self.calls.write_bytes(b"calls")
        (self.run / "basecalls/basecalls.vcf.gz.csi").write_bytes(b"index")
        raw = self.smoke / "raw_control/raw_control.vcf.gz"
        raw.write_bytes(b"raw")
        self.write(self.run / "basecalls/basecalls.json", {
            "output_sha256": digest(self.calls), "source_sha256": "source", "chromosome": 1,
            "genotypes_from": "mpileup_PL_and_AD_only", "input_imputed_GT_DS_GP_used": False,
            "sites": 1, "called": {"S_Chane-1.DG": 1},
        })
        self.write(self.run / "source/download.json", {"sha256": "source"})
        self.write(self.run / "verification/verification.json", {
            "basecalls_sha256": digest(self.calls), "raw_control_sha256": digest(raw),
            "selected_sha256": digest(self.smoke / "metadata/selected.json"), "chromosome": 1,
            "concordance": 0.999, "n_compared": 1000, "scientific_ancestry_launch_authorized": False,
        })
        self.write(self.run / "basecalls/likelihood_query.command.json",
                   ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%AD]\n"])

    def write(self, path, record):
        path.write_text(json.dumps(record))

    def run_verify(self, line):
        process = Mock(stdout=io.StringIO(line))
        process.wait.return_value = 0
        process.poll.return_value = 0
        with patch("verify_public_sgdp.subprocess.check_output", side_effect=["S_Chane-1.DG\n", "1\n"]), \
                patch("verify_public_sgdp.subprocess.Popen", return_value=process):
            verify(self.run, self.project, 1, self.smoke, self.checks / "chr1.json")

    def test_non22_chromosome_and_separate_metadata_are_verified(self):
        self.run_verify("1\t100\tA\tC\t0/1\t60\t40\t60,0,60\t20,20\n")
        result = json.loads((self.checks / "chr1.json").read_text())
        self.assertEqual(result["chromosome"], 1)
        self.assertEqual(result["genotype_cells_independently_rederived"], 1)
        self.assertFalse(result["scientific_ancestry_launch_authorized"])
        self.assertIn("metadata/selected.json", result["shared_sha256"])
        self.assertNotIn("metadata/selected.json", result["sha256"])

    def test_wrong_chromosome_and_tampered_genotypes_fail(self):
        for line in ("22\t100\tA\tC\t0/1\t60\t40\t60,0,60\t20,20\n",
                     "1\t100\tA\tC\t0/0\t60\t40\t60,0,60\t20,20\n"):
            with self.subTest(line=line), self.assertRaises(ValueError):
                self.run_verify(line)

    def test_missing_index_or_changed_callability_cannot_pass(self):
        process = Mock(stdout=io.StringIO("1\t100\tA\tC\t0/1\t60\t40\t60,0,60\t20,20\n"))
        process.wait.return_value = 0
        process.poll.return_value = 0
        with patch("verify_public_sgdp.subprocess.check_output", side_effect=["S_Chane-1.DG\n", "2\n"]), \
                patch("verify_public_sgdp.subprocess.Popen", return_value=process), \
                self.assertRaisesRegex(ValueError, "index disagrees"):
            verify(self.run, self.project, 1, self.smoke, self.checks / "chr1.json")

    def full_reports(self):
        reports = []
        qc = []
        for chromosome in range(1, 23):
            report = self.checks / f"chr{chromosome}.json"
            self.write(report, {
                "status": "published_genotypes_verified_not_ancestry_execution_ready",
                "chromosome": chromosome, "sample_ids": ["S_Chane-1.DG"], "sites": 1,
                "genotype_cells_independently_rederived": 1, "original_raw_control_concordance": 0.999,
                "original_raw_control_compared": 1000, "scientific_ancestry_launch_authorized": False,
                "sha256": {}, "shared_sha256": {"metadata/selected.json": digest(self.smoke / "metadata/selected.json")},
            })
            reports.append(report)
            qc.append({"chromosome": chromosome, "concordance": 0.999, "n_compared": 1000})
        self.write(self.run / "summary/full_import.json", {
            "chromosomes": list(range(1, 23)), "reused_chr22": str(self.smoke),
            "new_chromosome_root": str(self.run / "chromosomes"),
            "scientific_ancestry_launch_authorized": False, "raw_call_qc": qc,
        })
        for root in (self.run, self.checks):
            (root / "source_sha256.txt").write_text("")
        return reports

    def test_all22_coverage_and_post_cleanup_durable_hashes(self):
        reports = self.full_reports()
        output = self.checks / "summary/published_verification.json"
        with self.assertRaisesRegex(ValueError, "per autosome"):
            collect(self.run, self.smoke, self.checks, reports[:-1], output)
        collect(self.run, self.smoke, self.checks, reports, output)
        result = json.loads(output.read_text())
        self.assertEqual(result["genotype_cells_independently_rederived"], 22)
        self.assertFalse(result["scientific_ancestry_launch_authorized"])
        recheck(output, self.checks / "summary/recheck.json")
        (self.smoke / "metadata/selected.json").write_text("changed")
        with self.assertRaisesRegex(ValueError, "Durable evidence changed"):
            recheck(output, self.checks / "summary/recheck.json")

    def test_path_escape_in_report_is_rejected(self):
        reports = self.full_reports()
        record = json.loads(reports[0].read_text())
        record["sha256"]["../../../outside"] = "bad"
        self.write(reports[0], record)
        with self.assertRaisesRegex(ValueError, "escapes"):
            collect(self.run, self.smoke, self.checks, reports, self.checks / "bad.json")


if __name__ == "__main__":
    unittest.main()
