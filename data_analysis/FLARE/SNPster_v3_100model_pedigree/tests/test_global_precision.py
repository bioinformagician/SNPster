import gzip
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_global_precision
from screen_global_precision import digest, global_frame, local_hash, precision_source, verify_inputs


class GlobalPrecisionTests(unittest.TestCase):
    def test_source_transform_is_export_only_and_rejects_unexpected_source(self):
        source = (
            "unchanged accumulator\nString[] df3 = df3();\nfloat sum = 0f;\n"
            "float factor = 1000f/sum;\n"
            "sb.append(df3[ (int) Math.rint(factor*fa.get(j)) ]);\n"
        )
        result = precision_source(source)
        self.assertTrue(result.startswith("unchanged accumulator\n"))
        self.assertIn("Double.toString(factor*fa.get(j))", result)
        self.assertIn("double factor = 1.0/sum;", result)
        with self.assertRaises(ValueError):
            precision_source(result)
        with self.assertRaises(ValueError):
            precision_source(source + "\nfloat sum = 0f;")

    def test_local_hash_ignores_metadata_but_checks_calls_and_sample_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.gz", Path(tmp) / "b.gz"
            with gzip.open(a, "wb") as out:
                out.write(b"##date=old\n#CHROM\tA\tB\n22\t0|1\t1|1\n")
            with gzip.open(b, "wb") as out:
                out.write(b"##date=new\n#CHROM\tA\tB\n22\t0|1\t1|1\n")
            self.assertEqual(local_hash(a), local_hash(b))
            with gzip.open(b, "wb") as out:
                out.write(b"#CHROM\tA\tB\n22\t1|1\t1|1\n")
            self.assertNotEqual(local_hash(a), local_hash(b))

    def test_global_parser_retains_small_probabilities_and_rejects_invalid_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "global.tsv"
            path.write_text("SAMPLE\tanc_0\tanc_1\nA\t0.9997\t0.0003\n")
            np.testing.assert_allclose(global_frame(path).to_numpy(), [[.9997, .0003]])
            compressed = Path(tmp) / "global.anc.gz"
            with gzip.open(compressed, "wt") as out:
                out.write(path.read_text())
            np.testing.assert_array_equal(global_frame(path), global_frame(compressed))
            for row in ("A\t0\t0", "A\t1\t-1", "A\tNaN\t1", "A\t1\t0\nA\t1\t0"):
                path.write_text("SAMPLE\tanc_0\tanc_1\n" + row + "\n")
                with self.assertRaises(ValueError):
                    global_frame(path)

    def test_archived_script_override_verifies_original_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / "inference_screen.py"
            script.write_text("original inference implementation")
            provenance = {"input_sha256": {screen_global_precision.__file__: digest(script)}}
            verify_inputs(provenance, script)
            script.write_text("changed implementation")
            with self.assertRaises(ValueError):
                verify_inputs(provenance, script)
            with self.assertRaises(ValueError):
                verify_inputs({"input_sha256": {}}, script)


if __name__ == "__main__":
    unittest.main()
