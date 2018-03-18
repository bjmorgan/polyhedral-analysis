import unittest
import numpy as np
from polyhedral_analysis.symmetry_measure import SymmetryMeasure
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from itertools import permutations

class TestSymmetryMeasureInit( unittest.TestCase ):

    def test_symmetry_measure_is_initialised( self ):
        reference_points = np.array( [ [ 1.0, 2.0, 3.0 ], [ 2.0, 3.0, 4.0 ] ] )
        string = 'polyhedron'
        symmetry_measure = SymmetryMeasure( reference_points=reference_points,
                                            string=string )
        np.testing.assert_array_equal( symmetry_measure.reference_points, reference_points )
        self.assertEqual( symmetry_measure.string, string ) 

    def test_symmetry_measure_from_name( self ):
        octahedral_sm = SymmetryMeasure.from_name( 'Octahedron' )
        expected_reference_points = [[0.0, 0.0, 1.0], 
                                     [0.0, 0.0, -1.0], 
                                     [1.0, 0.0, 0.0], 
                                     [-1.0, 0.0, 0.0], 
                                     [0.0, 1.0, 0.0], 
                                     [0.0, -1.0, 0.0]]
        np.testing.assert_array_equal( octahedral_sm.reference_points, expected_reference_points )
        self.assertEqual( octahedral_sm.string, 'Octahedron' )

class TestSymmetryMeasure( unittest.TestCase ):

    def setUp( self ):
        self.sm = SymmetryMeasure.from_name( 'Octahedron' )
  
    def test_minimum_symmetry_measure( self ):
        # The symmetry measure for a perfect octahedron is 0.0.
        vertex_coordinates = np.array([[0.0, 0.0, 2.0], 
                                       [0.0, 0.0, -2.0], 
                                       [2.0, 0.0, 0.0], 
                                       [0.0, 2.0, 0.0], 
                                       [-2.0, 0.0, 0.0], 
                                       [0.0, -2.0, 0.0]])
        ag = AbstractGeometry( central_site=np.array([0.0, 0.0, 0.0]), 
                               bare_coords=vertex_coordinates)
        self.assertEqual( self.sm.minimum_symmetry_measure( ag=ag ), 0.0 )

    def test_symmetry_measures_from_coordination( self ):
        from polyhedral_analysis.symmetry_measure import symmetry_measures_from_coordination
        expected_reference_points = { 
            'Octahedron': [ [ 0.0,  0.0,  1.0 ],
                            [ 0.0,  0.0, -1.0 ],
                            [ 1.0,  0.0,  0.0 ],
                            [-1.0,  0.0,  0.0 ],
                            [ 0.0,  1.0,  0.0 ],
                            [ 0.0, -1.0,  0.0 ] ],
            'Trigonal prism': [ [-0.654654, -0.377964,  0.654654 ],
                                [ 0.654654, -0.377964,  0.654654 ],
                                [ 0.      ,  0.755929,  0.654654 ],
                                [-0.654653, -0.377964, -0.654653 ], 
                                [ 0.654653, -0.377964, -0.654653 ], 
                                [ 0.0,       0.755928, -0.654653 ] ] }
        for key, value in expected_reference_points.items():
            np.testing.assert_array_almost_equal( symmetry_measures_from_coordination[ len( value ) ][ key ].reference_points, value )

    

if __name__ == '__main__':
    unittest.main()
