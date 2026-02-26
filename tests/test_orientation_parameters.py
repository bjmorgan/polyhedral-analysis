import unittest

import numpy as np
import math

from polyhedral_analysis.orientation_parameters import cos_theta, oct_rotational_order_parameter

class OrientationParametersTestCase( unittest.TestCase ):

    def test_cos_theta_one( self ):
        a = np.array( [ 0.0, 0.0, 1.0 ] )
        b = np.array( [ 1.0, 0.0, 0.0 ] )
        self.assertEqual( cos_theta( a, b ), 0.0 )

    def test_cos_theta_two( self ):
        a = np.array( [ 0.0, 0.0, 1.0 ] )
        b = np.array( [ 0.0, 0.0, 1.0 ] )
        self.assertEqual( cos_theta( a, b ), 1.0 )

    def test_cos_theta_three( self ):
        a = np.array( [ 0.0, 0.0, 1.0 ] )
        b = np.array( [ 0.0, 1.0, 1.0 ] )
        self.assertTrue( cos_theta( a, b ) - math.sqrt(2)/2.0 < 1e-10 )

class TestOctRotationalOrderParameter(unittest.TestCase):

    def test_perfect_aligned_octahedron_returns_one(self):
        # Vertices along <100> directions, centred on origin.
        points = np.array([
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ])
        self.assertAlmostEqual(oct_rotational_order_parameter(points), 1.0)


if __name__ == '__main__':
    unittest.main()
