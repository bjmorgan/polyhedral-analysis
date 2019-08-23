import unittest
import numpy as np
from polyhedral_analysis.rotation_analyser import RotationAnalyser
from unittest.mock import Mock
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry

class TestRotationAnalyserInit( unittest.TestCase ):

    def test_rotation_analyser_is_initialised(self):
        reference_points = np.array( [[ 1.0, -1.0,  1.0],
                                      [-1.0, -1.0, -1.0],
                                      [ 1.0,  1.0, -1.0],
                                      [-1.0,  1.0,  1.0]] )
        ra = RotationAnalyser(reference_points=reference_points)
        self.assertEqual(ra.reference_points.shape, (1,4,3))
        np.testing.assert_array_equal(ra.reference_points, np.array([reference_points]))

    def test_rotation_analyser_is_initialised_with_multiple_reference_geometries(self):
        reference_points = np.array( [[[ 1.0, -1.0,  1.0],
                                       [-1.0, -1.0, -1.0]],
                                      [[ 1.0,  1.0, -1.0],
                                       [-1.0,  1.0,  1.0]]] )
        ra = RotationAnalyser(reference_points=reference_points)
        self.assertEqual(ra.reference_points.shape, (2,2,3))
        np.testing.assert_array_equal(ra.reference_points, reference_points)

    def test_rotation_analyser_init_raises_error_for_misformed_reference_points(self):
        reference_points = np.array( [[[ 1.0, -1.0,  1.0]],
                                      [[-1.0, -1.0, -1.0],
                                       [ 1.0,  1.0, -1.0],
                                       [-1.0,  1.0,  1.0]]] )
        with self.assertRaises(ValueError):
            ra = RotationAnalyser(reference_points=reference_points)

class TestRotationAnalyser( unittest.TestCase ):

    def setUp(self):
        reference_points = np.array( [[ 1.0, -1.0,  1.0],
                                      [-1.0, -1.0, -1.0],
                                      [ 1.0,  1.0, -1.0],
                                      [-1.0,  1.0,  1.0]] )
        all_points = np.array([ reference_points, reference_points*(-1.0) ])
        self.ra = RotationAnalyser(reference_points=all_points)
       
    def test_polyhedron_orientation(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.abstract_geometry = Mock(spec=AbstractGeometry)
        mock_polyhedron.abstract_geometry.points_wocs_csc = Mock(return_value='foo')
        self.ra.discrete_orientation = Mock(return_value='bar')
        orientation = self.ra.polyhedron_orientation( mock_polyhedron )
        self.ra.discrete_orientation.assert_called_once_with('foo')
        self.assertEqual( orientation, 'bar' )

if __name__ == '__main__':
    unittest.main()

