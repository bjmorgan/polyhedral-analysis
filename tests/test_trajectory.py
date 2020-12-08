import unittest
import copy
from unittest.mock import Mock, mock_open, patch
from polyhedral_analysis.trajectory import Trajectory
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.configuration import Configuration
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.core.structure import Structure

class TestTrajectoryInit( unittest.TestCase ):

    def test_trajectory_is_initialised(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        mock_configurations = [Mock(spec=Configuration), Mock(spec=Configuration)]
        trajectory = Trajectory(structures=mock_structures,
                                configurations=mock_configurations)
        self.assertEqual(trajectory._structures, mock_structures)
        self.assertEqual(trajectory._configurations, mock_configurations)

    def test_trajectory_is_initialised_with_config_numbers(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        mock_configurations = [Mock(spec=Configuration), Mock(spec=Configuration)]
        config_numbers = [10, 20]
        trajectory = Trajectory(structures=mock_structures,
                                configurations=mock_configurations)
        self.assertEqual(trajectory._structures, mock_structures)
        self.assertEqual(trajectory._configurations, mock_configurations)

    def test_trajectory_init_raises_TypeError_if_not_given_Structures(self):
        mock_structures = [ Mock( spec=Structure ), 'foo' ]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        with self.assertRaises( TypeError ):
            trajectory = Trajectory( structures=mock_structures,
                                     configurations=mock_configurations )

    def test_trajectory_init_raises_TypeError_if_not_given_Configurations(self):
        mock_structures = [ Mock( spec=Structure ), Mock( spec=Structure)]
        mock_configurations = [ 'foo', Mock( spec=Configuration ) ]
        with self.assertRaises( TypeError ):
            trajectory = Trajectory( structures=mock_structures,
                                     configurations=mock_configurations )

class TestTrajectory(unittest.TestCase):

    def setUp(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        mock_configurations = [ Mock( spec=Configuration ), Mock( spec=Configuration ) ]
        self.trajectory = Trajectory( structures=mock_structures,
                                      configurations=mock_configurations )
        
    def test_extend(self):
        trajectory1 = self.trajectory
        trajectory2 = copy.deepcopy(trajectory1)
        trajectory1.extend(trajectory2)
        self.assertEqual(len(trajectory1.configurations), 4)

    def test___len__(self):
        trajectory = self.trajectory
        self.assertEqual(len(trajectory), len(trajectory.structures))

if __name__ == '__main__':
    unittest.main()
