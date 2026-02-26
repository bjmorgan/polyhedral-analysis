import unittest

import numpy as np

from polyhedral_analysis.csm import SymmetryMeasureResult, continuous_symmetry_measure


class TestContinuousSymmetryMeasurePerfectGeometry(unittest.TestCase):

    def test_perfect_octahedron_has_zero_csm(self):
        perfect = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        result = continuous_symmetry_measure(perfect, perfect)
        self.assertAlmostEqual(result.symmetry_measure, 0.0)

    def test_scaled_perfect_octahedron_has_zero_csm(self):
        perfect = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        result = continuous_symmetry_measure(perfect * 2.0, perfect)
        self.assertAlmostEqual(result.symmetry_measure, 0.0)


class TestContinuousSymmetryMeasureDistortedGeometry(unittest.TestCase):

    def test_distorted_octahedron_matches_expected_csm(self):
        """CSM for this distorted octahedron was verified against
        pymatgen's chemenv implementation."""
        distorted = np.array([
            [0.0, 0.0, 1.1],
            [0.0, 0.1, -1.0],
            [1.0, 0.0, 0.1],
            [-1.0, 0.1, 0.0],
            [0.0, 1.0, -0.1],
            [0.1, -1.0, 0.0],
        ])
        perfect = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        result = continuous_symmetry_measure(distorted, perfect)
        self.assertAlmostEqual(result.symmetry_measure, 0.892184544299399, places=10)
        self.assertAlmostEqual(result.scaling_factor, 0.9746359136170526, places=10)

    def test_rotation_matrix_is_orthogonal(self):
        distorted = np.array([
            [0.0, 0.0, 1.1],
            [0.0, 0.1, -1.0],
            [1.0, 0.0, 0.1],
            [-1.0, 0.1, 0.0],
            [0.0, 1.0, -0.1],
            [0.1, -1.0, 0.0],
        ])
        perfect = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        result = continuous_symmetry_measure(distorted, perfect)
        R = result.rotation_matrix
        np.testing.assert_array_almost_equal(R.T @ R, np.eye(3))


class TestContinuousSymmetryMeasureSinglePoint(unittest.TestCase):

    def test_single_point_returns_zero_csm(self):
        result = continuous_symmetry_measure(
            np.array([[1.0, 2.0, 3.0]]),
            np.array([[4.0, 5.0, 6.0]]),
        )
        self.assertEqual(result.symmetry_measure, 0.0)
        np.testing.assert_array_equal(result.rotation_matrix, np.eye(3))
        self.assertEqual(result.scaling_factor, 1.0)


class TestContinuousSymmetryMeasureDegenerateInput(unittest.TestCase):

    def test_all_zero_distorted_points_raises_value_error(self):
        zeros = np.zeros((6, 3))
        perfect = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        with self.assertRaises(ValueError):
            continuous_symmetry_measure(zeros, perfect)


class TestSymmetryMeasureResult(unittest.TestCase):

    def test_result_is_named_tuple(self):
        result = SymmetryMeasureResult(
            symmetry_measure=1.0,
            rotation_matrix=np.eye(3),
            scaling_factor=1.0,
        )
        self.assertEqual(result.symmetry_measure, 1.0)
        self.assertEqual(result.scaling_factor, 1.0)


if __name__ == '__main__':
    unittest.main()
