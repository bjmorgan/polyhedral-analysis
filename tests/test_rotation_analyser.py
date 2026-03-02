import itertools
import unittest
import numpy as np
from polyhedral_analysis.rotation_analyser import RotationAnalyser, _analyse_point_group
from polyhedral_analysis.csm import continuous_symmetry_measure
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

TETRAHEDRON_VERTICES = np.array([[1.0, -1.0, 1.0],
                                  [-1.0, -1.0, -1.0],
                                  [1.0, 1.0, -1.0],
                                  [-1.0, 1.0, 1.0]])


class TestAnalysePointGroup(unittest.TestCase):

    def test_tetrahedron_reduced_permutation_count(self):
        """Td has |G|=24, so 4!/24 = 1 reduced permutation."""
        reduced_perms, _ = _analyse_point_group(TETRAHEDRON_VERTICES)
        self.assertEqual(len(reduced_perms), 1)
        self.assertEqual(reduced_perms.shape[1], 4)

    def test_tetrahedron_proper_rotation_count(self):
        """Td has 12 proper rotations."""
        _, proper_rots = _analyse_point_group(TETRAHEDRON_VERTICES)
        self.assertEqual(len(proper_rots), 12)

    def test_proper_rotations_are_orthogonal_with_det_plus_one(self):
        """All proper rotation matrices should be orthogonal with det = +1."""
        _, proper_rots = _analyse_point_group(TETRAHEDRON_VERTICES)
        for R in proper_rots:
            np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-10)
            self.assertAlmostEqual(np.linalg.det(R), 1.0, places=10)


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
        rng = np.random.default_rng(42)
        distorted = self.ra.reference_points[0] + rng.normal(0, 0.05, self.ra.reference_points[0].shape)
        result = self.ra.discrete_orientation(distorted)
        self.assertEqual(result['reference_geometry_index'], 0)
        self.assertAlmostEqual(result['rotational_distance'], 0.0, places=1)
        self.assertLess(result['symmetry_measure'], 1.0)

    def test_discrete_orientation_identifies_inverted_reference(self):
        """A slightly distorted copy of reference orientation 1 should be identified as such."""
        rng = np.random.default_rng(123)
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

    def test_discrete_orientation_matches_brute_force(self):
        """Reduced approach should give the same minimum CSM as brute-force N! permutations."""
        rng = np.random.default_rng(42)
        distorted = self.ra.reference_points[0] + rng.normal(0, 0.05, self.ra.reference_points[0].shape)
        result = self.ra.discrete_orientation(distorted)
        # Brute force over all 4! = 24 permutations for reference 0
        points = distorted - np.mean(distorted, axis=0)
        best_csm = float('inf')
        for perm in itertools.permutations(range(4)):
            csm_result = continuous_symmetry_measure(
                points, self.ra.reference_points[0][list(perm)])
            if csm_result.symmetry_measure < best_csm:
                best_csm = csm_result.symmetry_measure
        self.assertAlmostEqual(result['symmetry_measure'], best_csm, places=10)

    def test_discrete_orientation_rejects_wrong_vertex_count(self):
        """discrete_orientation should raise ValueError for mismatched vertex count."""
        wrong_shape = np.zeros((3, 3))  # 3 vertices instead of 4
        with self.assertRaises(ValueError):
            self.ra.discrete_orientation(wrong_shape)

    def test_all_rotational_distances_shape(self):
        """all_rotational_distances should have length n_refs * n_proper."""
        points = self.ra.reference_points[0].copy()
        result = self.ra.discrete_orientation(points)
        n_refs = len(self.ra.reference_points)
        n_proper = len(self.ra.proper_rotations)
        self.assertEqual(len(result['all_rotational_distances']), n_refs * n_proper)

    def test_orientation_index_indexes_into_all_rotational_distances(self):
        """orientation_index should correspond to the minimum in all_rotational_distances."""
        points = self.ra.reference_points[0].copy()
        result = self.ra.discrete_orientation(points)
        idx = result['orientation_index']
        self.assertEqual(
            result['rotational_distance'],
            result['all_rotational_distances'][idx])

    def test_discrete_orientation_does_not_mutate_input(self):
        """The input points array should not be modified."""
        points = self.ra.reference_points[0] + 0.5  # offset from origin
        original = points.copy()
        self.ra.discrete_orientation(points)
        np.testing.assert_array_equal(points, original)

    def test_proper_rotation_of_reference_gives_zero_distance(self):
        """Applying a proper rotation to the reference should yield zero rotational distance."""
        R = self.ra.proper_rotations[1]  # non-identity proper rotation
        rotated = (R @ self.ra.reference_points[0].T).T
        result = self.ra.discrete_orientation(rotated)
        self.assertAlmostEqual(result['symmetry_measure'], 0.0, places=5)
        self.assertAlmostEqual(result['rotational_distance'], 0.0, places=5)


if __name__ == '__main__':
    unittest.main()
