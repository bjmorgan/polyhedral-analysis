import unittest
from unittest.mock import Mock
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.polyhedron_trajectory import PolyhedronTrajectory

class TestPolyhedronTrajectoryInit( unittest.TestCase ):

    def test_polyhedron_trajectory_is_initialised_without_config_numbers( self ):
        mock_polyhedra = [ Mock( spec=CoordinationPolyhedron ),
                           Mock( spec=CoordinationPolyhedron ) ]
        expected_config_numbers = [ 1, 2 ]
        mock_polyhedron_trajectory = PolyhedronTrajectory( polyhedra=mock_polyhedra )
        self.assertEqual( mock_polyhedron_trajectory.polyhedra, mock_polyhedra )
        self.assertEqual( mock_polyhedron_trajectory.config_numbers, expected_config_numbers )

    def test_polyhedron_trajectory_is_initialised_with_config_numbers( self ):
        mock_polyhedra = [ Mock( spec=CoordinationPolyhedron ), 
                           Mock( spec=CoordinationPolyhedron ) ]
        config_numbers = [ 2, 4 ]
        mock_polyhedron_trajectory = PolyhedronTrajectory( polyhedra=mock_polyhedra,
                                                           config_numbers=config_numbers )
        self.assertEqual( mock_polyhedron_trajectory.polyhedra, mock_polyhedra )
        self.assertEqual( mock_polyhedron_trajectory.config_numbers, config_numbers )

    def test_polyhedron_trajectory_raises_ValueError_with_incorrect_config_numbers( self ):
        mock_polyhedra = [ Mock( spec=CoordinationPolyhedron ),
                           Mock( spec=CoordinationPolyhedron ) ]
        config_numbers = [ 2, 4, 6 ]
        with self.assertRaises( ValueError ):
            mock_polyhedron_trajectory = PolyhedronTrajectory( polyhedra=mock_polyhedra,
                                                               config_numbers=config_numbers )

if __name__ == '__main__':
    unittest.main() 
