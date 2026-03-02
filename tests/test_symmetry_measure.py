import math
import unittest
import numpy as np
from itertools import permutations
from polyhedral_analysis.symmetry_measure import SymmetryMeasure
from polyhedral_analysis.csm import continuous_symmetry_measure

class TestSymmetryMeasureInit( unittest.TestCase ):

    def test_symmetry_measure_is_initialised( self ):
        reference_points = np.array( [ [ 1.0, 2.0, 3.0 ], [ 2.0, 3.0, 4.0 ] ] )
        name = 'polyhedron'
        symmetry_measure = SymmetryMeasure( reference_points=reference_points,
                                            name=name )
        np.testing.assert_array_equal( symmetry_measure.reference_points, reference_points )
        self.assertEqual( symmetry_measure.name, name ) 

    def test_symmetry_measure_from_name( self ):
        octahedral_sm = SymmetryMeasure.from_name( 'Octahedron' )
        expected_reference_points = [[0.0, 0.0, 1.0], 
                                     [0.0, 0.0, -1.0], 
                                     [1.0, 0.0, 0.0], 
                                     [-1.0, 0.0, 0.0], 
                                     [0.0, 1.0, 0.0], 
                                     [0.0, -1.0, 0.0]]
        np.testing.assert_array_equal( octahedral_sm.reference_points, expected_reference_points )
        self.assertEqual( octahedral_sm.name, 'Octahedron' )

class TestSymmetryMeasure( unittest.TestCase ):

    def setUp( self ):
        self.sm = SymmetryMeasure.from_name( 'Octahedron' )
  
    def test_minimum_symmetry_measure( self ):
        # The symmetry measure for a perfect octahedron is 0.0.
        vertex_vectors = np.array([[0.0, 0.0, 2.0],
                                   [0.0, 0.0, -2.0],
                                   [2.0, 0.0, 0.0],
                                   [0.0, 2.0, 0.0],
                                   [-2.0, 0.0, 0.0],
                                   [0.0, -2.0, 0.0]])
        self.assertEqual( self.sm.minimum_symmetry_measure( vertex_vectors ), 0.0 )

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



class TestReducedPermutations(unittest.TestCase):

    def test_tetrahedron_has_one_representative(self):
        sm = SymmetryMeasure.from_name('Tetrahedron')
        self.assertEqual(len(sm.permutations), 1)

    def test_octahedron_has_15_representatives(self):
        sm = SymmetryMeasure.from_name('Octahedron')
        self.assertEqual(len(sm.permutations), 15)

    def test_reduced_permutations_cover_all_orbits(self):
        """For the trigonal bipyramid, verify that the number of
        representatives times the symmetry group order equals N!."""
        sm = SymmetryMeasure.from_name('Trigonal bipyramid')
        # Trigonal bipyramid has |G|=12, so 120/12 = 10 representatives.
        self.assertEqual(len(sm.permutations), 10)
        self.assertEqual(len(sm.permutations) * 12, math.factorial(5))

    def test_reduced_gives_same_csm_as_brute_force(self):
        """Verify that the reduced-permutation CSM matches brute-force
        for a distorted octahedron."""
        sm = SymmetryMeasure.from_name('Octahedron')
        # A distorted octahedron (vertices shifted from ideal positions).
        distorted = np.array([
            [0.0, 0.0, 1.1],
            [0.0, 0.1, -1.0],
            [1.0, 0.0, 0.1],
            [-1.0, 0.1, 0.0],
            [0.0, 1.0, -0.1],
            [0.1, -1.0, 0.0],
        ])
        reduced_csm = sm.minimum_symmetry_measure(distorted)
        # Brute-force over all 6! permutations.
        brute_force_csm = min(
            continuous_symmetry_measure(np.array(p), sm.reference_points).symmetry_measure
            for p in permutations(distorted)
        )
        self.assertAlmostEqual(reduced_csm, brute_force_csm, places=10)


if __name__ == '__main__':
    unittest.main()
