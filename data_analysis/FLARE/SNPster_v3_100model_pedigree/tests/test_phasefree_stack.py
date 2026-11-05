from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
import screen_phasefree_stack as stack
import test_supervised_genotypes as fixtures


class PhasefreeStackTests(unittest.TestCase):
    def test_decoders_train_only_on_disjoint_clean_or_imputed_features(self):
        with tempfile.TemporaryDirectory() as tmp:
            axes, cohort, pure, ped = fixtures.SupervisedGenotypeTests().fixture(Path(tmp))
            full = np.eye(7)[np.arange(len(axes)) % 7]
            feature_rows = np.flatnonzero(~axes.group.eq("reference"))
            values = np.tile(full[feature_rows], (1, 2))
            import screen_supervised_genotypes as shared
            with patch.object(shared, "training_table", return_value=(pure, ped)):
                before, _ = stack.fit_decoders(values, feature_rows, axes, cohort)
                changed = values.copy()
                changed[-1] = 1e5
                after, _ = stack.fit_decoders(changed, feature_rows, axes, cohort)
                for mode in before:
                    np.testing.assert_array_equal(before[mode].predict(values[:3]),
                                                  after[mode].predict(values[:3]))
                with self.assertRaisesRegex(ValueError, "nontraining"):
                    stack.fit_decoders(values[1:], feature_rows[1:], axes, cohort)

    def test_feature_construction_combines_two_separate_seven_source_views(self):
        classifier = np.array([[.5, .5, 0, 0, 0, 0, 0]])
        frequency = np.array([[.4, .6, 0, 0, 0, 0, 0]])
        with (patch.object(stack, "pca_predict", return_value=classifier),
              patch.object(stack, "frequency_predict", return_value=(frequency, None, None))):
            values = stack.features({"classifier": {}, "frequency": {}}, np.ones((1, 4)))
        np.testing.assert_array_equal(values, np.c_[classifier, frequency])
        self.assertEqual(values.shape, (1, 14))


if __name__ == "__main__":
    unittest.main()
