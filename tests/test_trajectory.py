import unittest
import copy
from unittest.mock import Mock, mock_open, patch
from polyhedral_analysis.trajectory import Trajectory, read_config_numbers_from_xdatcar
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.configuration import Configuration
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen import Structure

class TestTrajectoryInit( unittest.TestCase ):

    def test_trajectory_is_initialised( self ):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure ) ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        trajectory = Trajectory( structures=mock_structures,
                                 configurations=mock_configurations )
        self.assertEqual( trajectory._structures, mock_structures )
        self.assertEqual( trajectory._configurations, mock_configurations )
        self.assertEqual( trajectory.config_numbers, [ 1, 2 ] )

    def test_trajectory_is_initialised_with_config_numbers( self ):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure ) ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        config_numbers = [ 10, 20 ]
        trajectory = Trajectory( structures=mock_structures,
                                 configurations=mock_configurations,
                                 config_numbers=config_numbers )
        self.assertEqual( trajectory._structures, mock_structures )
        self.assertEqual( trajectory._configurations, mock_configurations )
        self.assertEqual( trajectory.config_numbers, config_numbers )

    def test_trajectory_init_raises_TypeError_if_not_given_Structures( self ):
        mock_structures = [ Mock( spec=Structure ), 'foo' ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        with self.assertRaises( TypeError ):
            trajectory = Trajectory( structures=mock_structures,
                                     configurations=mock_configurations )

    def test_trajectory_init_raises_TypeError_if_not_given_Structures( self ):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure ) ]
        mock_configurations = [ 'foo', Mock( spec=Configuration ) ]
        with self.assertRaises( TypeError ):
            trajectory = Trajectory( structures=mock_structures,
                                     configurations=mock_configurations )

    def test_trajectory_init_raises_ValueError_if_config_numbers_are_incorrect_length( self ):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure ) ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        config_numbers = [ 1, 2, 3 ]
        with self.assertRaises( ValueError ):
            trajectory = Trajectory( structures=mock_structures,
                                     configurations=mock_configurations,
                                     config_numbers=config_numbers )

class TestTrajectory( unittest.TestCase ):

    def setUp( self ):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure ) ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        self.trajectory = Trajectory( structures=mock_structures,
                                      configurations=mock_configurations )
        
    def test_extend( self ):
        trajectory1 = self.trajectory
        trajectory2 = copy.deepcopy( trajectory1 )
        trajectory1.extend( trajectory2 )
        self.assertEqual( trajectory1.config_numbers, [ 1, 2, 3, 4 ] )
        self.assertEqual( len( trajectory1.configurations ), 4 )
        self.assertEqual( [ c.config_number for c in trajectory1.configurations ], [ 1, 2, 3, 4 ] )

    def test_extend_with_offset( self ):
        trajectory1 = self.trajectory
        trajectory2 = copy.deepcopy( trajectory1 )
        trajectory1.extend( trajectory2, offset=5 )
        self.assertEqual( trajectory1.config_numbers, [ 1, 2, 8, 9 ] )
        self.assertEqual( len( trajectory1.configurations ), 4 )
        self.assertEqual( [ c.config_number for c in trajectory1.configurations ], [ 1, 2, 8, 9 ] )

    def test_config_numbers_getter( self ):
        trajectory = self.trajectory
        trajectory.configurations[0].config_number = 5
        trajectory.configurations[1].config_number = 10
        self.assertEqual( trajectory.config_numbers, [ 5, 10 ] )

class TestTrajectoryHelperFunctions( unittest.TestCase ):

    def test_read_config_numbers_from_xdatcar( self ):
        example_file = "Direct configuration=  10\nDirect configuration=   20\n"
        with patch( 'polyhedral_analysis.trajectory.zopen', mock_open( read_data=example_file ), create=True ) as m:
            self.assertEqual( read_config_numbers_from_xdatcar( 'foo' ), [ 10, 20 ] )

if __name__ == '__main__':
    unittest.main()
