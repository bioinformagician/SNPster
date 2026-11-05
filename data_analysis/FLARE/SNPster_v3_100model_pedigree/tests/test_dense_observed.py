import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import benchmark_dense_observed as benchmark


class DenseObservedTests(unittest.TestCase):
    def test_vectorized_gt_parser_preserves_order_and_pairing(self):
        actual = benchmark.decode_calls("0|1\t1|0\t1|1\t0|0", 4)
        np.testing.assert_array_equal(actual, [[0, 1], [1, 0], [1, 1], [0, 0]])
        for bad in ("0/1\t1|0", "0|.\t1|0", "0|2\t1|0", "0|1\t1", "0|1\t1|0\t0|0"):
            with self.assertRaises(ValueError):
                benchmark.decode_calls(bad, 2)

    def test_whole_families_are_purged_when_any_founder_is_unavailable(self):
        donors = {"a": ["x", "y"], "b": ["x", "absent"], "c": ["y"]}
        self.assertEqual(benchmark.safe_families(donors, ["x", "y"]), {"a", "c"})
        self.assertEqual(benchmark.safe_families(donors, ["x"]), set())


if __name__ == "__main__":
    unittest.main()
