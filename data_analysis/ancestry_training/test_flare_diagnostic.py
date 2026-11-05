import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np

from diagnose_flare import crossed_parameters, kinship_coefficients, native_vector


class FlareDiagnostics(unittest.TestCase):
    def test_identical_sample_has_half_kinship(self):
        values, zero = kinship_coefficients([[0, 1, 1, 2]], [0, 1, 1, 2])
        np.testing.assert_array_equal(values, [.5])
        np.testing.assert_array_equal(zero, [0])

    def test_kinship_rejects_invalid_dosages_and_undefined_mass(self):
        for reference, target in [([[0, 3]], [0, 1]), ([[0, 2]], [0, 2]), ([[1, 1]], [1])]:
            with self.assertRaises(ValueError):
                kinship_coefficients(reference, target)

    def test_parent_like_sharing_has_quarter_kinship_without_ibs0(self):
        values, zero = kinship_coefficients([[0, 1, 0, 1]], [0, 1, 1, 2])
        np.testing.assert_array_equal(values, [.25])
        np.testing.assert_array_equal(zero, [0])

    def test_crossing_aligns_rows_and_columns_by_name(self):
        destination = {"ancestries": ["A", "B"], "panels": ["A", "B"],
                       "p": np.eye(2), "rho": np.array([1., 2.]), "theta": np.zeros((2, 2))}
        donor = {"ancestries": ["B", "A"], "panels": ["B", "A"],
                 "p": np.array([[.8, .2], [.3, .7]]), "rho": np.array([20., 10.]),
                 "theta": np.full((2, 2), .01)}
        result = crossed_parameters(destination, donor, ["p", "rho"])
        np.testing.assert_array_equal(result["p"], [[.7, .3], [.2, .8]])
        np.testing.assert_array_equal(result["rho"], [10., 20.])
        np.testing.assert_array_equal(result["theta"], destination["theta"])
        np.testing.assert_array_equal(destination["rho"], [1., 2.])
        with self.assertRaises(ValueError):
            crossed_parameters(destination, donor, ["unknown"])

    def test_native_global_axes_are_named_not_positional(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "global.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("SAMPLE\tB\tA\nperson\t0.75\t0.25\n")
            np.testing.assert_array_equal(native_vector(path, ["A", "B"]), [.25, .75])
            with self.assertRaises(ValueError):
                native_vector(path, ["A", "missing"])


if __name__ == "__main__":
    unittest.main()
