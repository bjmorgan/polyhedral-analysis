import unittest
from polyhedral_analysis.orientation_parameters import cos_theta
import numpy as np
import math

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

if __name__ == '__main__':
    unittest.main()
