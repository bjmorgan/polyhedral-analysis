import unittest
import numpy as np
from polyhedral_analysis.polyhedra_recipe import (PolyhedraRecipe,
                                                  matching_sites,
                                                  polyhedra_from_distance_cutoff)
from pymatgen.core import Structure, Lattice
from polyhedral_analysis.atom import Atom
from unittest.mock import Mock, patch

class TestPolyhedraRecipeInit(unittest.TestCase):

    def test_polyhedra_recipe_init_raises_ValueError_for_invalid_method(self):
        invalid_method = 'foo'
        with self.assertRaises(ValueError):
            PolyhedraRecipe(method=invalid_method, central_atoms='foo', vertex_atoms='bar')

    def test_polyhdra_recipe_init(self):
        recipe_args = {'method': 'distance cutoff',
                       'central_atoms': 'S',
                       'vertex_atoms': 'Li',
                       'coordination_cutoff': 5.0}
        with patch('polyhedral_analysis.polyhedra_recipe.generator_from_atom_argument') as mock_generator_from_atom_argument:
            mock_generators = [Mock(), Mock()]
            mock_generator_from_atom_argument.side_effect = mock_generators
            recipe = PolyhedraRecipe(**recipe_args)
            self.assertEqual(recipe.method, recipe_args['method'])
            self.assertEqual(recipe._central_atom_list_generator, mock_generators[0])
            self.assertEqual(recipe._vertex_atom_list_generator, mock_generators[1])
            self.assertEqual(recipe._central_atom_list, None)
            self.assertEqual(recipe._vertex_atom_list, None)
            self.assertEqual(recipe.coordination_cutoff, recipe_args['coordination_cutoff'])
            self.assertEqual(recipe.vertex_graph_cutoff, None)
            self.assertEqual(recipe.n_neighbours, None)
            self.assertEqual(recipe.label, None)
            self.assertEqual(recipe.recalculate, True)

class TestPolyhedraRecipeFunctions(unittest.TestCase):

    def test_matching_sites(self):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.0],
                           [0.0, 0.5, 0.5],
                           [0.5, 0.0, 0.5]])
        atom_list = ['Li'] * len(coords)
        lattice = Lattice.from_parameters(a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)    
        ref_coords = np.array([[0.1, 0.1, 0.1],
                               [0.0, 0.4, 0.5]])
        ref_atom_list = ['Na'] * len(ref_coords)
        ref_structure = Structure(lattice, ref_atom_list, ref_coords)
        matched_sites = matching_sites(structure, ref_structure)
        self.assertEqual(len(matched_sites), 2)
        self.assertEqual(matched_sites[0], (structure[0], 0))
        self.assertEqual(matched_sites[1], (structure[2], 2))

    def test_matching_sites_with_species(self):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.0],
                           [0.0, 0.5, 0.5],
                           [0.5, 0.0, 0.5]])
        atom_list = ['Li', 'Mg', 'Mg', 'Mg']
        lattice = Lattice.from_parameters(a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        ref_coords = np.array([[0.1, 0.1, 0.1],
                               [0.0, 0.4, 0.5]])
        ref_atom_list = ['Na'] * len(ref_coords)
        ref_structure = Structure(lattice, ref_atom_list, ref_coords)
        matched_sites = matching_sites(structure, ref_structure, species=['Li'])
        self.assertEqual(len(matched_sites), 1)
        self.assertEqual(matched_sites[0], (structure[0], 0))

    def test_polyhedra_from_distance_cutoff_with_no_central_atoms_returns_empty_list(self):
        # If an empty list of central atoms is passed in, and empty list of polyhedra
        # should be returned.
        central_atoms = []
        vertex_atoms = [ Mock(spec=Atom), Mock(spec=Atom) ]
        cutoff = 1.0
        polyhedra = polyhedra_from_distance_cutoff( central_atoms=central_atoms,
                                                    vertex_atoms=vertex_atoms,
                                                    cutoff=cutoff )
        self.assertEqual( polyhedra, [] )

if __name__ == '__main__':
    unittest.main()

