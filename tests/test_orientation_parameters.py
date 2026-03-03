import unittest

import numpy as np
import math

from polyhedral_analysis.orientation_parameters import cos_theta, projection_xyz, oct_rotational_order_parameter

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

class TestProjectionXyz(unittest.TestCase):

    def test_vector_along_axis_returns_one(self):
        self.assertAlmostEqual(projection_xyz(np.array([1.0, 0.0, 0.0])), 1.0)
        self.assertAlmostEqual(projection_xyz(np.array([0.0, 1.0, 0.0])), 1.0)
        self.assertAlmostEqual(projection_xyz(np.array([0.0, 0.0, 1.0])), 1.0)

    def test_vector_along_diagonal_returns_minus_one_third(self):
        self.assertAlmostEqual(
            projection_xyz(np.array([1.0, 1.0, 1.0])), -1.0 / 3.0)


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

    def test_45_degree_rotated_octahedron_returns_one_third(self):
        # Rotate around z by 45 degrees: 4 in-plane vertices contribute 0,
        # 2 axial vertices contribute 1; sum = 2/6 = 1/3.
        s = math.sqrt(2) / 2
        points = np.array([
            [ s,  s, 0.0],
            [-s, -s, 0.0],
            [-s,  s, 0.0],
            [ s, -s, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ])
        self.assertAlmostEqual(
            oct_rotational_order_parameter(points), 1.0 / 3.0)


if __name__ == '__main__':
    unittest.main()
