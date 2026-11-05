import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from pyliftover import LiftOver

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "pedigree_benchmark/scripts"))
from normalize_sgdp import Fasta, IndependentChain, index_fasta, mapping, normalize, swap_call, verify_normalized
from sgdp_genetic_qc import audit_pairs


class NormalizationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.fasta("source", "ACGT" * 5)
        self.target = self.fasta("target", "NNNNN" + "ACGT" * 5 + "NNNNN")
        self.chain = self.root / "chain.gz"
        self.write_chain("+")
        self.protocol = dict(minimum_retained_fraction=0.90, minimum_called_fraction=0.90)

    def fasta(self, name, sequence):
        path = self.root / f"{name}.fa"
        path.write_text(">1\n" + "\n".join(sequence[i:i + 7] for i in range(0, len(sequence), 7)) + "\n")
        index_fasta(path)
        return path

    def write_chain(self, strand):
        with gzip.open(self.chain, "wt") as handle:
            handle.write(f"chain 100 chr1 20 + 0 20 chr1 30 {strand} 5 25 1\n20\n")

    def test_fasta_index_line_boundaries_and_coordinate_bounds(self):
        fasta = Fasta(self.source)
        try:
            self.assertEqual("".join(fasta.base("1", pos) for pos in range(1, 21)), "ACGT" * 5)
            for pos in (0, 21):
                with self.assertRaises(ValueError):
                    fasta.base("1", pos)
        finally:
            fasta.close()

    def test_independent_chain_matches_forward_reverse_endpoints_and_gaps(self):
        for strand in ("+", "-"):
            self.write_chain(strand)
            direct, library = IndependentChain(self.chain, "1"), LiftOver(str(self.chain))
            for zero in range(-1, 22):
                expected = [hit[:3] for hit in library.convert_coordinate("chr1", zero) or []]
                self.assertEqual(direct.convert_coordinate("chr1", zero), expected)
            self.assertEqual(direct.convert_coordinate("chr1", 0)[0][1], 5 if strand == "+" else 24)

    def test_mapping_uses_one_based_vcf_positions_and_complements_reverse_alleles(self):
        for strand, expected in (("+", (6, "A", "C", "+", False)),
                                 ("-", (25, "T", "G", "-", False))):
            self.write_chain(strand)
            source, target = Fasta(self.source), Fasta(self.target)
            try:
                self.assertEqual(mapping("1", 1, "A", "C", LiftOver(str(self.chain)), source, target),
                                 (expected, None))
                self.assertEqual(mapping("1", 1, "C", "A", LiftOver(str(self.chain)), source, target),
                                 (None, "source_reference_mismatch"))
            finally:
                source.close()
                target.close()

    def test_ref_alt_swap_reorders_raw_evidence_and_preserves_missingness(self):
        self.assertEqual(swap_call("0/0:60:20:0,60,255:20,0"), "1/1:60:20:255,60,0:0,20")
        self.assertEqual(swap_call("0/1:60:40:60,0,60:20,20"), "0/1:60:40:60,0,60:20,20")
        self.assertEqual(swap_call("./.:.:.:.:."), "./.:.:.:.:.")

    def source_vcf(self):
        path = self.root / "original.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tperson\n"
            "1\t1\t.\tA\tC\t.\t.\t.\tGT:GQ:DP:PL:AD\t0/0:60:20:0,60,255:20,0\n"
            "1\t2\t.\tC\tA\t.\t.\t.\tGT:GQ:DP:PL:AD\t0/1:60:40:60,0,60:20,20\n"
            "1\t20\t.\tT\tA\t.\t.\t.\tGT:GQ:DP:PL:AD\t1/1:60:20:255,60,0:0,20\n")
        subprocess.run(["bgzip", str(path)], check=True)
        return path.with_suffix(".vcf.gz")

    def test_real_normalization_and_independent_verification_with_ref_swap(self):
        original = self.source_vcf()
        self.target = self.fasta("target_swap", "NNNNN" + "C" + ("ACGT" * 5)[1:] + "NNNNN")
        output = self.root / "normalized"
        normalize(original, "1", self.source, self.target, self.chain, self.protocol, output)
        previous = Path.cwd()
        import os
        try:
            os.chdir(self.root)
            verify_normalized(original, output, "1", self.source, self.target, self.chain,
                              self.root / "verified.json")
        finally:
            os.chdir(previous)
        record = json.loads((self.root / "verified.json").read_text())
        self.assertEqual(record["sites"], 3)
        self.assertEqual(record["genotype_cells_verified"], 3)
        self.assertEqual(json.loads((output / "normalization.json").read_text())["counts"]["ref_alt_swaps"], 1)
        self.assertFalse(record["scientific_ancestry_launch_authorized"])

    def test_king_real_table_column_contract_and_pair_audit(self):
        path = self.root / "king.vcf"
        header = ("##fileformat=VCFv4.2\n##contig=<ID=1>\n"
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcandidate\tpanel\tcandidate2\n")
        values = ("0/0", "0/1", "1/1")
        path.write_text(header + "".join(
            f"1\t{i + 1}\tS{i}\tA\tC\t.\t.\t.\tGT\t{values[i % 3]}\t{values[(i // 3) % 3]}\t{values[(i // 7) % 3]}\n"
            for i in range(200)))
        prefix = self.root / "king"
        ids = self.root / "require.id"
        ids.write_text("0\tcandidate\n0\tcandidate2\n")
        subprocess.run([
            "plink2", "--vcf", str(path), "--const-fid", "0", "--make-king-table",
            "cols=id,nsnp,hethet,ibs0,ibs1,kinship", "--king-table-require", str(ids),
            "--threads", "1", "--memory", "768", "--out", str(prefix),
        ], check=True, stdout=subprocess.DEVNULL)
        with prefix.with_suffix(".kin0").open() as handle:
            names = handle.readline().lstrip("#").split()
            records = [dict(zip(names, line.split())) for line in handle]
        policy = dict(minimum_pair_called_markers=100, king_kinship_exclusion_threshold=0.0442,
                      duplicate_concordance_threshold=0.995)
        result = audit_pairs(records, ["candidate", "panel", "candidate2"],
                             ["candidate", "candidate2"], policy)
        self.assertEqual(result["pairs_checked"], 3)
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            audit_pairs(records[:-1], ["candidate", "panel", "candidate2"], ["candidate", "candidate2"], policy)
        records[0]["NSNP"] = "1"
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            audit_pairs(records, ["candidate", "panel", "candidate2"], ["candidate", "candidate2"], policy)


if __name__ == "__main__":
    unittest.main()
