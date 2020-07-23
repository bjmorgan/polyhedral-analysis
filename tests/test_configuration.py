import unittest
from unittest.mock import Mock, patch
from pymatgen import Structure, Site
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom

class TestConfigurationInit(unittest.TestCase):

    def test_configuration_is_initialised(self):
        mock_sites = [Mock(spec=Site) for i in range(4)]
        for m, s in zip(mock_sites, ['A', 'B', 'C', 'D']):
            m.species_string = s
        mock_structure = Mock(spec=Structure)
        mock_structure.sites = mock_sites
        mock_polyhedra = [Mock(spec=CoordinationPolyhedron) for i in range(2)]
        mock_atoms = [Mock(spec=Atom) for i in range(4)]
        for i, m in enumerate(mock_atoms):
            m.index = i
        mock_polyhedra[0].central_atom = mock_atoms[0]
        mock_polyhedra[0].vertices = [mock_atoms[1]]
        mock_polyhedra[1].central_atom = mock_atoms[2]
        mock_polyhedra[1].vertices = [mock_atoms[3]]
        mock_recipe = Mock(spec=PolyhedraRecipe)
        mock_recipe.find_polyhedra = Mock(return_value=mock_polyhedra)
        with patch('polyhedral_analysis.configuration.Atom') as mock_Atom:
            mock_Atom.side_effect = mock_atoms
            configuration = Configuration(
                structure=mock_structure, recipes=[mock_recipe])
        self.assertEqual(configuration.polyhedra, mock_polyhedra)
        self.assertEqual(configuration.central_atoms, [
                         mock_atoms[0], mock_atoms[2]])
        self.assertEqual(configuration.coordination_atoms,
                         [mock_atoms[1], mock_atoms[3]])
        self.assertEqual(configuration.atoms, mock_atoms)

class TestConfiguration(unittest.TestCase):

    def setUp(self):
        mock_sites = [Mock(spec=Site) for i in range(4)]
        for m, s in zip(mock_sites, ['A', 'B', 'C', 'D']):
            m.species_string = s
        mock_structure = Mock(spec=Structure)
        mock_structure.sites = mock_sites
        mock_polyhedra = [Mock(spec=CoordinationPolyhedron) for i in range(2)]
        mock_atoms = [Mock(spec=Atom) for i in range(4)]
        for i, m in enumerate(mock_atoms):
            m.index = i
        mock_polyhedra[0].central_atom = mock_atoms[0]
        mock_polyhedra[0].vertices = [mock_atoms[1]]
        mock_polyhedra[0].label = 'foo'
        mock_polyhedra[1].central_atom = mock_atoms[2]
        mock_polyhedra[1].vertices = [mock_atoms[3]]
        mock_polyhedra[1].label = 'bar'
        mock_recipe = Mock(spec=PolyhedraRecipe)
        mock_recipe.find_polyhedra = Mock(return_value=mock_polyhedra)
        with patch('polyhedral_analysis.configuration.Atom') as mock_Atom:
            mock_Atom.side_effect = mock_atoms
            self.configuration = Configuration(
                structure=mock_structure, recipes=[mock_recipe])
        self.mock_atoms = mock_atoms

    def test_coordination_atom_by_index(self):
        self.assertEqual(
            self.configuration.coordination_atom_by_index(3), self.mock_atoms[3])

    def test_coordination_atom_by_index_returns_none_if_missing(self):
        self.assertEqual(
            self.configuration.coordination_atom_by_index(2), None)

    def test_polyhedra_labels(self):
        self.assertEqual(self.configuration.polyhedra_labels, ['foo', 'bar'])

    def test_polyhedra_by_label(self):
        p = self.configuration.polyhedra_by_label('foo')
        self.assertEqual(len(p), 1)
        self.assertEqual(p[0].label, 'foo')
        self.assertEqual(p[0].central_atom.index, 0)


if __name__ == '__main__':
    unittest.main()
