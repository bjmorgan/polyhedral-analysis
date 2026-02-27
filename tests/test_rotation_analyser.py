import unittest
import numpy as np
from polyhedral_analysis.rotation_analyser import RotationAnalyser
from unittest.mock import Mock
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron

class TestRotationAnalyserInit(unittest.TestCase):

    def test_rotation_analyser_is_initialised(self):
        reference_points = np.array([[1.0, -1.0, 1.0],
                                     [-1.0, -1.0, -1.0],
                                     [1.0, 1.0, -1.0],
                                     [-1.0, 1.0, 1.0]])
        ra = RotationAnalyser(reference_points=reference_points)
        self.assertEqual(ra.reference_points.shape, (1, 4, 3))
        np.testing.assert_array_equal(
            ra.reference_points, np.array([reference_points]))

    def test_rotation_analyser_is_initialised_with_multiple_reference_geometries(self):
        reference_points = np.array([[[1.0, -1.0, 1.0],
                                      [-1.0, -1.0, -1.0]],
                                     [[1.0, 1.0, -1.0],
                                      [-1.0, 1.0, 1.0]]])
        ra = RotationAnalyser(reference_points=reference_points)
        self.assertEqual(ra.reference_points.shape, (2, 2, 3))
        np.testing.assert_array_equal(ra.reference_points, reference_points)

    def test_rotation_analyser_init_raises_error_for_misformed_reference_points(self):
        reference_points = np.array([[[1.0, -1.0, 1.0]],
                                     [[-1.0, -1.0, -1.0],
                                      [1.0, 1.0, -1.0],
                                      [-1.0, 1.0, 1.0]]], dtype='object')
        with self.assertRaises(ValueError):
            ra = RotationAnalyser(reference_points=reference_points)

class TestRotationAnalyser(unittest.TestCase):

    def setUp(self):
        reference_points = np.array([[1.0, -1.0, 1.0],
                                     [-1.0, -1.0, -1.0],
                                     [1.0, 1.0, -1.0],
                                     [-1.0, 1.0, 1.0]])
        all_points = np.array([reference_points, reference_points * (-1.0)])
        self.ra = RotationAnalyser(reference_points=all_points)

    def test_polyhedron_orientation(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.vertex_vectors = Mock(return_value='foo')
        self.ra.discrete_orientation = Mock(return_value='bar')
        orientation = self.ra.polyhedron_orientation(mock_polyhedron)
        mock_polyhedron.vertex_vectors.assert_called_once_with(reference='central_atom')
        self.ra.discrete_orientation.assert_called_once_with('foo')
        self.assertEqual(orientation, 'bar')

    def test_discrete_orientation_identifies_correct_reference(self):
        """A slightly distorted copy of reference orientation 0 should be identified as such."""
        rng = np.random.RandomState(42)
        distorted = self.ra.reference_points[0] + rng.normal(0, 0.05, self.ra.reference_points[0].shape)
        result = self.ra.discrete_orientation(distorted)
        self.assertEqual(result['reference_geometry_index'], 0)
        self.assertAlmostEqual(result['rotational_distance'], 0.0, places=1)
        self.assertLess(result['symmetry_measure'], 1.0)

    def test_discrete_orientation_identifies_inverted_reference(self):
        """A slightly distorted copy of reference orientation 1 should be identified as such."""
        rng = np.random.RandomState(123)
        distorted = self.ra.reference_points[1] + rng.normal(0, 0.05, self.ra.reference_points[1].shape)
        result = self.ra.discrete_orientation(distorted)
        self.assertEqual(result['reference_geometry_index'], 1)
        self.assertAlmostEqual(result['rotational_distance'], 0.0, places=1)

    def test_discrete_orientation_perfect_geometry_has_zero_csm(self):
        """A perfect (undistorted) reference geometry should have CSM close to zero."""
        points = self.ra.reference_points[0].copy()
        result = self.ra.discrete_orientation(points)
        self.assertAlmostEqual(result['symmetry_measure'], 0.0, places=5)
        self.assertAlmostEqual(result['rotational_distance'], 0.0, places=5)


if __name__ == '__main__':
    unittest.main()
