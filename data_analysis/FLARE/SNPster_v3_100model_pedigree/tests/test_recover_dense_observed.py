import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from recover_dense_observed import FAILED, output_prefix


class RecoveryTests(unittest.TestCase):
    def test_only_known_failed_call_uses_a_new_output_prefix(self):
        work = Path("/benchmark/task/candidate")
        candidate, sample, stage = FAILED
        self.assertEqual(output_prefix(work, candidate, sample, stage),
                         work / f"{sample}_{stage}_native_retry")
        for policy, person, phase in (
                ("flare_array_fixed", sample, stage),
                (candidate, "different", stage),
                (candidate, sample, "inherited_phase")):
            self.assertEqual(output_prefix(work, policy, person, phase), work / f"{person}_{phase}")


if __name__ == "__main__":
    unittest.main()
