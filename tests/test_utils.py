import unittest
import polyhedral_analysis.utils as utils

class TestUtilsFunctions( unittest.TestCase ):

    def test_flatten( self ):
        nested = [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]
        flat = [ 1, 2, 3, 4, 5, 6 ]
        self.assertEqual( flat, utils.flatten( nested ) )

if __name__ == '__main__':
    unittest.main()
