from collections import Counter
import sys
from pathlib import Path
import unittest
import gzip
import json
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from prepare_dense_public_support import finalize, select_candidates, verify_raw
from audit_production_ancestry import digest


class DensePublicSupportTests(unittest.TestCase):
    def test_public_alias_exclusions_and_all_eight_source_counts(self):
        sources = ("EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA")
        protocol = dict(
            allowed_access={"Embargo": "FullyPublic", "category": "X"},
            candidate_population_sources={source: source for source in sources},
            minimum_metadata_candidates_per_source=2)
        records = []
        for source in sources:
            for number in range(3):
                sample = f"{source}_{number}"
                records.append(dict(
                    Illumina_ID=f"instrument_{sample}", Sample_ID=sample,
                    **{"Sample_ID(Aliases)": f"alias_{sample}"},
                    SGDP_ID=sample, Population_ID=source,
                    Embargo="FullyPublic", category="X"))
        people, counts = select_candidates(records, {"alias_EUR_0"}, protocol)
        self.assertEqual(counts["EUR"], 2)
        self.assertEqual(len(people), 23)
        self.assertNotIn("EUR_0.DG", {person["sample_id"] for person in people})
        records[-1]["Embargo"] = "SignedLetter"
        _, counts = select_candidates(records, {"alias_EUR_0"}, protocol)
        self.assertEqual(counts["SSA"], 2)
        records[-2]["category"] = "Y"
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            select_candidates(records, {"alias_EUR_0"}, protocol)

    def test_dai_is_eas_and_admixed_oceania_not_assumed_pure(self):
        protocol = json.loads((Path(__file__).resolve().parents[1] /
                               "config/dense_public_support.json").read_text())
        mapping = protocol["candidate_population_sources"]
        self.assertEqual(mapping["Dai"], "EAS")
        self.assertEqual(mapping["Australian"], "OCE")
        self.assertNotIn("Hawaiian", mapping)
        self.assertNotIn("Maori", mapping)
        self.assertEqual(set(Counter(mapping.values())), {"EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA"})

    def test_raw_evidence_is_rederived_and_only_identical_duplicates_removed(self):
        project = Path(__file__).resolve().parents[1]
        protocol = project / "config/sgdp_public_import.json"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = root / "input.vcf.gz"
            selection = root / "selected.json"
            summary = root / "summary.json"
            selection.write_text(json.dumps({"people": [{"sample_id": "a"}]}))
            header = (
                "##fileformat=VCFv4.2\n##contig=<ID=22>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n")
            first = "22\t100\t.\tA\tC\t.\t.\t.\tGT:GQ:DP:PL:AD\t0/1:60:20:60,0,80:10,10\n"
            second = "22\t200\t.\tG\tT\t.\t.\t.\tGT:GQ:DP:PL:AD\t0/0:60:20:0,60,80:20,0\n"
            with gzip.open(calls, "wt") as handle:
                handle.write(header + first + first + second)
            summary.write_text(json.dumps(dict(output_sha256=digest(calls), chromosome=22, sites=3)))
            verify_raw(calls, summary, selection, protocol, root / "verified")
            record = json.loads((root / "verified/verified.json").read_text())
            self.assertEqual(record["retained_sites"], 2)
            self.assertEqual(record["removed_strict_redundant_records"], 1)
            self.assertEqual(record["genotype_cells_independently_rederived"], 3)
            conflict = first.replace("0/1:60:20:60,0,80:10,10", "1/1:60:20:80,60,0:0,20")
            with gzip.open(calls, "wt") as handle:
                handle.write(header + first + conflict + second)
            summary.write_text(json.dumps(dict(output_sha256=digest(calls), chromosome=22, sites=3)))
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                verify_raw(calls, summary, selection, protocol, root / "conflicting")

    def test_finalizer_binds_indexes_and_requires_same_eligible_people_across_autosomes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "policy.json"
            source.write_text("{}")
            (root / "source_sha256.txt").write_text(f"{digest(source)}  {source}\n")
            people = [dict(sample_id=f"{source}_{number}", source=source)
                      for source in ("EUR", "EAS", "SAS", "SEA", "GME", "NAM", "OCE", "SSA")
                      for number in range(3)]
            for chromosome in range(1, 23):
                run = root / ("smoke" if chromosome == 22 else "full")
                metadata = run / "metadata"
                metadata.mkdir(parents=True, exist_ok=True)
                (metadata / "selected.json").write_text(json.dumps({
                    "people": people, "input_sha256": {str(source): digest(source)}}))
                normalized = run / f"normalized/chr{chromosome}/normalized.vcf.gz"
                observed = run / f"raw_observed/chr{chromosome}/observed/basecalls.vcf.gz"
                for path in (normalized, observed):
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text("genotypes")
                index = Path(str(normalized) + ".csi")
                index.write_text("index")
                summary = normalized.parent / "normalization.json"
                summary.write_text("{}")
                (observed.parent / "verified.json").write_text(json.dumps(dict(
                    status="every_observed_GT_GQ_DP_independently_rederived_from_PL_AD",
                    output_sha256=digest(observed))))
                verification = run / f"verification/chr{chromosome}.verified.json"
                verification.parent.mkdir(exist_ok=True)
                verification.write_text(json.dumps(dict(
                    normalized_sha256=digest(normalized), source_vcf_sha256=digest(observed),
                    normalized_index_sha256=digest(index),
                    normalization_summary_sha256=digest(summary))))
                report = run / f"density/chr{chromosome}/density/density.json"
                report.parent.mkdir(parents=True)
                report.write_text(json.dumps(dict(
                    chromosome=chromosome, normalized_sha256=digest(normalized),
                    minimum_real_mask_coverage=0.9,
                    people=[dict(person=person, mask_support={
                        mask: {"fraction": 0.9} for mask in ("lonni", "frederik")})
                        for person in people])))
            finalize(root / "smoke", False)
            finalize(root / "full", False)
            record = json.loads((root / "full/published_verification.json").read_text())
            self.assertEqual(record["chromosomes"], list(range(1, 23)))
            self.assertIn(str(root / "smoke/normalized/chr22/normalized.vcf.gz"),
                          record["published_sha256"])
            finalize(root / "full", True)
            index = root / "full/normalized/chr1/normalized.vcf.gz.csi"
            index.write_text("changed")
            with self.assertRaisesRegex(ValueError, "independent checks"):
                finalize(root / "full", False)
            index.write_text("index")
            for chromosome, number in ((1, 0), (2, 1)):
                report = root / f"full/density/chr{chromosome}/density/density.json"
                data = json.loads(report.read_text())
                for row in data["people"]:
                    if row["person"]["sample_id"] == f"EUR_{number}":
                        row["mask_support"]["lonni"]["fraction"] = 0.89
                report.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "entire scope"):
                finalize(root / "full", False)


if __name__ == "__main__":
    unittest.main()
