import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from audit_production_ancestry import digest
from dense_donor_qc import classify, merge_candidates, verify_merge, verify_panel
from sgdp_genetic_qc import genotype_rows


def vcf(path, identities, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.with_suffix("")
    with plain.open("w") as handle:
        handle.write('##fileformat=VCFv4.2\n##contig=<ID=22>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                     "\t".join(identities) + "\n")
        for position, ref, alt, genotypes in rows:
            handle.write(f"22\t{position}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t" +
                         "\t".join(genotypes) + "\n")
    with path.open("wb") as handle:
        subprocess.run(["bgzip", "-c", str(plain)], stdout=handle, check=True)
    subprocess.run(["bcftools", "index", "-f", str(path)], check=True)
    return path


def fixture(root):
    new = vcf(root / "new.vcf.gz", ["new", "control"], [
        (100, "A", "C", ["0/1", "1/1"]), (200, "G", "T", ["./.", "0/0"]),
        (300, "A", "C", ["1/1", "0/1"])])
    old = vcf(root / "old.vcf.gz", ["old", "control"], [
        (100, "A", "C", ["0/0", "1/1"]), (200, "G", "T", ["1/1", "0/0"]),
        (300, "A", "G", ["0/1", "0/1"])])
    markers = vcf(root / "markers.vcf.gz", ["marker"], [
        (100, "A", "C", ["0/1"]), (200, "G", "T", ["0/0"]), (300, "A", "C", ["1/1"])])
    proof = root / "normalized.json"
    proof.write_text(json.dumps(dict(status="coordinate_allele_and_genotype_preservation_verified",
                                    normalized_sha256=digest(new), chromosome=22,
                                    sample_ids=["new", "control"])))
    paths = [new, old, Path(str(new) + ".csi"), Path(str(old) + ".csi"), proof]
    plan = dict(new_ids=["new"], historical_ids=["old", "control"], public_ids=["new", "old", "control"],
                recalled_control_id="control", minimum_control_shared_calls=2, minimum_control_concordance=0.995,
                marker_vcf=str(markers), jobs={"22": dict(
                    new=str(new), old=str(old), verification=str(proof),
                    input_sha256={str(path): digest(path) for path in paths})})
    return plan, new, old, proof


class DenseDonorQCTests(unittest.TestCase):
    def test_exact_merge_preserves_missingness_and_excludes_wrong_alleles(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan, new, old, proof = fixture(root)
            folder = root / "merged"
            merge_candidates(plan, 22, new, old, proof, folder)
            output = root / "verified.json"
            verify_merge(plan, 22, new, old, proof, folder, output)
            record = json.loads(output.read_text())
            self.assertEqual(record["sites"], 2)
            self.assertEqual(record["genotype_cells_verified"], 6)
            self.assertEqual(record["control_concordance"], 1)
            rows = dict(genotype_rows(folder / "candidates.vcf.gz", folder / "positions.tsv", plan["public_ids"]))
            self.assertEqual(rows[("22", 200, "G", "T")], ["./.", "1/1", "0/0"])
            self.assertNotIn(("22", 300, "A", "C"), rows)
            vcf(folder / "candidates.vcf.gz", plan["public_ids"], [
                (100, "A", "C", ["1/1", "0/0", "1/1"]),
                (200, "G", "T", ["./.", "1/1", "0/0"])])
            with self.assertRaisesRegex(ValueError, "genotypes differ"):
                verify_merge(plan, 22, new, old, proof, folder, output)

    def test_merge_rejects_sample_axis_changes_and_dropped_loci(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan, new, old, proof = fixture(root)
            folder = root / "merged"
            merge_candidates(plan, 22, new, old, proof, folder)
            vcf(folder / "candidates.vcf.gz", ["old", "new", "control"], [
                (100, "A", "C", ["0/1", "0/0", "1/1"])])
            with self.assertRaisesRegex(ValueError, "sample axis"):
                verify_merge(plan, 22, new, old, proof, folder, root / "proof.json")
            vcf(folder / "candidates.vcf.gz", plan["public_ids"], [
                (100, "A", "C", ["0/1", "0/0", "1/1"])])
            with self.assertRaisesRegex(ValueError, "loci differ"):
                verify_merge(plan, 22, new, old, proof, folder, root / "proof.json")

    def test_panel_verification_checks_reference_historical_and_public_genotypes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan, new, old, proof = fixture(root)
            merged = root / "merged"
            merge_candidates(plan, 22, new, old, proof, merged)
            fixed = vcf(root / "production/beagle_references/fixed.vcf.gz", ["fixed"], [
                (100, "A", "C", ["0|0"]), (200, "G", "T", ["1|1"])])
            historical = vcf(root / "historical.vcf.gz", ["historical"], [
                (100, "A", "C", ["0/1"]), (200, "G", "T", ["0/0"])])
            plan.update(all_ids=["fixed", "historical", *plan["public_ids"]], fixed_ids=["fixed"],
                        groups=[dict(path=str(historical), ids=["historical"])])
            folder = root / "panel"
            calls = vcf(folder / "panel.vcf.gz", plan["all_ids"], [
                (100, "A", "C", ["0/0", "0/1", "0/1", "0/0", "1/1"]),
                (200, "G", "T", ["1/1", "0/0", "./.", "1/1", "0/0"])])
            report = dict(chromosome=22, sample_ids=plan["all_ids"], candidate_ids=plan["public_ids"],
                          panel_sha256=digest(calls), sites=2,
                          called=dict(fixed=2, historical=2, new=1, old=2, control=2))
            (folder / "panel.json").write_text(json.dumps(report))
            preflight = dict(references=[dict(chromosome=22, reference=str(fixed),
                                             files={str(fixed): {"sha256": digest(fixed)}})])
            output = root / "panel_verified.json"
            verify_panel(plan, 22, folder, merged, root / "production", preflight, output)
            self.assertEqual(json.loads(output.read_text())["genotype_cells_verified"], 10)
            vcf(calls, plan["all_ids"], [
                (100, "A", "C", ["1/1", "0/1", "0/1", "0/0", "1/1"]),
                (200, "G", "T", ["1/1", "0/0", "./.", "1/1", "0/0"])])
            report["panel_sha256"] = digest(calls)
            (folder / "panel.json").write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, "genotype changed"):
                verify_panel(plan, 22, folder, merged, root / "production", preflight, output)

    def test_fixed_reference_exclusion_identity_bridges_and_family_capacity(self):
        candidates = ["evaluation", "training", "alias1", "alias2", "oce1", "oce2"]
        plan = dict(candidate_ids=candidates, coverage_eligible_ids=candidates, fixed_ids=["reference"],
                    sources={identity: "OCE" if identity.startswith("oce") else "EUR" for identity in candidates},
                    historical_roles={"legacy": ["pure_training"]})
        def pair(first, second, concordance=1):
            return dict(first=first, second=second, genotype_concordance=concordance)
        related = [pair("evaluation", "reference"), pair("training", "legacy"),
                   pair("alias1", "old_development"), pair("alias2", "old_development"),
                   pair("oce1", "oce2", 0.85)]
        result = classify(plan, related, {"duplicate_concordance_threshold": 0.995})
        self.assertEqual(result["fixed_reference_direct_related_exclusions"], ["evaluation"])
        self.assertIn("training", result["eligible_candidate_ids"])
        self.assertEqual(result["independent_component_capacity_by_source"]["EUR"], 2)
        self.assertEqual(result["independent_component_capacity_by_source"]["OCE"], 1)
        self.assertTrue(result["historical_identity_bridges"][0]["must_be_absent_from_new_fitting_roles"])
        self.assertFalse(result["ancestry_models_authorized"])
        self.assertFalse(result["donor_roles_frozen"])
        connected = classify(plan, related + [pair("alias1", "evaluation", 0.85)],
                             {"duplicate_concordance_threshold": 0.995})
        self.assertEqual(connected["fixed_reference_connected_exclusions"],
                         ["alias1", "alias2", "evaluation"])
        self.assertEqual(connected["independent_component_capacity_by_source"]["EUR"], 1)

    def test_conflicting_identity_labels_are_not_reported_as_reference_overlap(self):
        plan = dict(candidate_ids=["a", "b"], coverage_eligible_ids=["a", "b"], fixed_ids=[],
                    sources={"a": "EUR", "b": "EAS"}, historical_roles={})
        result = classify(plan, [dict(first="a", second="b", genotype_concordance=1)],
                          {"duplicate_concordance_threshold": 0.995})
        self.assertEqual(result["fixed_reference_direct_related_exclusions"], [])
        self.assertEqual(result["identity_label_conflict_exclusions"], ["a", "b"])
        self.assertEqual(result["eligible_candidate_ids"], [])


if __name__ == "__main__":
    unittest.main()
