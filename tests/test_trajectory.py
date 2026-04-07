import unittest
import copy
from unittest.mock import Mock, mock_open, patch
from polyhedral_analysis.trajectory import Trajectory
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.configuration import Configuration
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.core.structure import Structure
from pymatgen.core import Lattice

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
        self.assertEqual(len(trajectory1.structures), 4)

    def test___len__(self):
        trajectory = self.trajectory
        self.assertEqual(len(trajectory), len(trajectory.structures))

    def test_add_returns_new_trajectory(self):
        trajectory2 = copy.deepcopy(self.trajectory)
        combined = self.trajectory + trajectory2
        self.assertEqual(len(combined), 4)

    def test_add_does_not_modify_originals(self):
        trajectory2 = copy.deepcopy(self.trajectory)
        _ = self.trajectory + trajectory2
        self.assertEqual(len(self.trajectory), 2)
        self.assertEqual(len(trajectory2), 2)

    def test_add_returns_not_implemented_for_non_trajectory(self):
        result = self.trajectory.__add__('not a trajectory')
        self.assertIs(result, NotImplemented)


class TestTrajectoryFromStructures(unittest.TestCase):

    def test_returns_correct_number_of_configurations(self):
        lattice = Lattice.cubic(10.0)
        species = ['Ti'] + ['O'] * 6
        coords = [
            [0.5, 0.5, 0.5],
            [0.6, 0.5, 0.5], [0.4, 0.5, 0.5],
            [0.5, 0.6, 0.5], [0.5, 0.4, 0.5],
            [0.5, 0.5, 0.6], [0.5, 0.5, 0.4],
        ]
        structure = Structure(lattice, species, coords)
        recipe = PolyhedraRecipe(
            method='nearest neighbours',
            central_atoms='Ti',
            vertex_atoms='O',
            n_neighbours=6)
        trajectory = Trajectory.from_structures(
            [structure, structure], [recipe])
        self.assertEqual(len(trajectory), 2)

    def test_polyhedra_constructed_per_recipes(self):
        lattice = Lattice.cubic(10.0)
        species = ['Ti'] + ['O'] * 6
        coords = [
            [0.5, 0.5, 0.5],
            [0.6, 0.5, 0.5], [0.4, 0.5, 0.5],
            [0.5, 0.6, 0.5], [0.5, 0.4, 0.5],
            [0.5, 0.5, 0.6], [0.5, 0.5, 0.4],
        ]
        structure = Structure(lattice, species, coords)
        recipe = PolyhedraRecipe(
            method='nearest neighbours',
            central_atoms='Ti',
            vertex_atoms='O',
            n_neighbours=6)
        trajectory = Trajectory.from_structures(
            [structure], [recipe])
        self.assertEqual(len(trajectory.configurations[0].polyhedra), 1)
        self.assertEqual(
            trajectory.configurations[0].polyhedra[0].coordination_number, 6)

    def test_raises_type_error_for_non_bool_progress(self):
        with self.assertRaises(TypeError):
            Trajectory.from_structures([], [], progress='notebook')


if __name__ == '__main__':
    unittest.main()
