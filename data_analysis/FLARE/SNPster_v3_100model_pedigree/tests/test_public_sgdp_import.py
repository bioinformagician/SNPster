import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from import_sgdp_public import call_from_likelihoods, select_public
from verify_public_sgdp import expected_call
from prepare_sgdp_full_import import collect, require_verified


class PublicSGDPImportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = json.loads((PROJECT / "config/sgdp_public_import.json").read_text())

    def row(self, access="FullyPublic", category="X", alias="person"):
        return dict(Illumina_ID="DNA1", Sample_ID=alias, **{
            "Sample_ID(Aliases)": alias, "SGDP_ID": "S_Chane-1",
            "Population_ID": "Chane", "Embargo": access, "category": category,
        })

    def test_only_public_and_all_alias_disjoint_people_are_eligible(self):
        self.assertEqual(len(select_public([self.row()], set(), self.protocol)), 1)
        for forbidden in ({"person"}, {"DNA1"}, {"S_Chane-1"}):
            self.assertEqual(select_public([self.row()], forbidden, self.protocol), [])
        for access, category in (("SignedLetter", "Y"), ("FullyPublic", "DO_NOT_USE")):
            self.assertEqual(select_public([self.row(access, category)], set(), self.protocol), [])

    def test_calls_use_only_likelihoods_and_read_depth(self):
        for pl, expected in (("0,60,255", "0/0"), ("60,0,60", "0/1"), ("255,60,0", "1/1")):
            self.assertEqual(call_from_likelihoods(pl, "20,20", self.protocol)[0], expected)
        for pl, ad in (("0,0,100", "20,20"), ("0,29,100", "20,20"),
                       ("0,60,255", "4,5"), ("0,60,255", "60,41"), (".", ".")):
            self.assertEqual(call_from_likelihoods(pl, ad, self.protocol)[0], "./.")
        self.assertEqual(call_from_likelihoods("0,30,100", "5,5", self.protocol)[0], "0/0")

    def test_invalid_likelihood_shapes_and_values_fail(self):
        for pl, ad in (("0,60", "20,20"), ("-1,60,100", "20,20"), ("0,60,100", "20")):
            with self.assertRaises(ValueError):
                call_from_likelihoods(pl, ad, self.protocol)

    def test_independent_verifier_reconstructs_calls_and_missingness(self):
        for pl, ad in (("0,60,255", "40,0"), ("60,0,60", "20,20"),
                       ("255,60,0", "0,40"), ("0,0,100", "20,20"),
                       ("0,29,100", "20,20"), ("0,30,100", "5,5"),
                       ("0,60,255", "60,41"), (".", ".")):
            with self.subTest(pl=pl, ad=ad):
                self.assertEqual(expected_call(pl, ad, self.protocol),
                                 call_from_likelihoods(pl, ad, self.protocol))

    def test_geographic_oceania_is_not_automatically_oceanian_ancestry(self):
        sources = self.protocol["candidate_population_sources"]
        self.assertEqual(sources["Australian"], "OCE")
        for population in ("Dusun", "Igorot", "Hawaiian", "Maori", "Saharawi"):
            self.assertNotIn(population, sources)
        self.assertEqual(self.protocol["required_ancestry_resolutions"],
                         ["eight_category", "seven_category"])
        self.assertFalse(self.protocol["scientific_ancestry_launch_authorized"])

    def test_full_import_requires_a_verified_smoke(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "published_verification.json").write_text(json.dumps({
                "status": "not_verified", "chromosome": 22,
                "scientific_ancestry_launch_authorized": False,
            }))
            with self.assertRaisesRegex(ValueError, "verified public"):
                require_verified(root, PROJECT)

    def test_full_import_covers_all_autosomes_and_preserves_unresolved_gates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "verification").mkdir()
            status = "raw_control_concordance_passed_not_ancestry_execution_ready"
            (root / "verification/verification.json").write_text(json.dumps({
                "status": status, "chromosome": 22,
            }))
            reports = []
            for chromosome in range(1, 22):
                path = root / f"{chromosome}.json"
                path.write_text(json.dumps(dict(status=status, chromosome=chromosome)))
                reports.append(path)
            with self.assertRaisesRegex(ValueError, "exactly one"):
                collect(root, reports[:-1], root)
            with patch("prepare_sgdp_full_import.write_json") as write:
                collect(root, reports, root)
            result = write.call_args.args[1]
            self.assertEqual(result["chromosomes"], list(range(1, 23)))
            self.assertFalse(result["scientific_ancestry_launch_authorized"])
            self.assertFalse(result["independent_all_autosome_genotype_verification_complete"])


if __name__ == "__main__":
    unittest.main()
