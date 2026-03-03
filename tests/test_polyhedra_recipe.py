import unittest
import numpy as np
from polyhedral_analysis.polyhedra_recipe import (PolyhedraRecipe,
                                                  matching_sites,
                                                  polyhedra_from_distance_cutoff,
                                                  generator_from_atom_argument,
                                                  polyhedra_from_nearest_neighbours,
                                                  polyhedra_from_closest_centre,
                                                  polyhedra_from_atom_indices)
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from pymatgen.core import Structure, Lattice
from pymatgen.core.sites import PeriodicSite
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
        polyhedra = polyhedra_from_distance_cutoff(central_atoms=central_atoms,
                                                   vertex_atoms=vertex_atoms,
                                                   cutoff=cutoff)
        self.assertEqual(polyhedra, [])
        
class TestPolyhedraRecipe(unittest.TestCase):
    
    def setUp(self):
        self.mock_structure = Mock(spec=Structure)
        self.recipe = PolyhedraRecipe(method='distance cutoff',
                                      central_atoms='Li',
                                      vertex_atoms=['O', 'F'],
                                      coordination_cutoff=2.5)
                                      
    def test_central_atom_list(self):
        mock_generator = Mock(return_value=[0, 1, 2])
        self.recipe._central_atom_list_generator = mock_generator        
        result = self.recipe.central_atom_list(self.mock_structure)        
        self.assertEqual(result, [0, 1, 2])
        mock_generator.assert_called_once_with(self.mock_structure)

    def test_vertex_atom_list(self):
        mock_generator = Mock(return_value=[3, 4, 5])
        self.recipe._vertex_atom_list_generator = mock_generator        
        result = self.recipe.vertex_atom_list(self.mock_structure)        
        self.assertEqual(result, [3, 4, 5])
        mock_generator.assert_called_once_with(self.mock_structure)

    def test_central_atom_list_no_recalculate(self):
        self.recipe.recalculate = False
        self.recipe._central_atom_list = [0, 1, 2]        
        result = self.recipe.central_atom_list()        
        self.assertEqual(result, [0, 1, 2])
    
    def test_find_polyhedra(self):
        mock_atoms = [Mock(spec=Atom, index=i) for i in range(6)]
        self.recipe.central_atom_list = Mock(return_value=[0, 1])
        self.recipe.vertex_atom_list = Mock(return_value=[2, 3, 4, 5])
        
        with patch('polyhedral_analysis.polyhedra_recipe.polyhedra_from_distance_cutoff') as mock_polyhedra_func:
            mock_polyhedra_func.return_value = [Mock(spec=CoordinationPolyhedron)]
            
            result = self.recipe.find_polyhedra(mock_atoms, self.mock_structure)
            
            mock_polyhedra_func.assert_called_once_with(
                central_atoms=[mock_atoms[0], mock_atoms[1]],
                vertex_atoms=[mock_atoms[2], mock_atoms[3], mock_atoms[4], mock_atoms[5]],
                cutoff=2.5,
                label=None
            )
            self.assertEqual(result, mock_polyhedra_func.return_value)

class TestPolyhedraGenerationFunctions(unittest.TestCase):
    def setUp(self):
        self.mock_central_atoms = [Mock(spec=Atom) for _ in range(2)]
        self.mock_vertex_atoms = [Mock(spec=Atom) for _ in range(4)]

        # Assign coordinates to atoms
        for i, atom in enumerate(self.mock_central_atoms + self.mock_vertex_atoms):
            atom.coords = np.array([i, i, i], dtype=float)

        # Implement distance method for mock atoms
        def distance(self, other):
            return np.linalg.norm(self.coords - other.coords)

        for atom in self.mock_central_atoms + self.mock_vertex_atoms:
            atom.distance = distance.__get__(atom)


class TestPolyhedraFromNearestNeighbours(unittest.TestCase):

    def test_selects_n_nearest_vertex_atoms(self):
        lattice = Lattice.cubic(10.0)
        central = Atom(0, PeriodicSite('Ti', [0.5, 0.5, 0.5], lattice))
        close_coords = [
            [0.6, 0.5, 0.5], [0.4, 0.5, 0.5],
            [0.5, 0.6, 0.5], [0.5, 0.4, 0.5],
            [0.5, 0.5, 0.6], [0.5, 0.5, 0.4],
        ]
        far_coords = [[0.8, 0.5, 0.5], [0.2, 0.5, 0.5]]
        vertex_atoms = [Atom(i + 1, PeriodicSite('O', c, lattice))
                        for i, c in enumerate(close_coords + far_coords)]
        polyhedra = polyhedra_from_nearest_neighbours(
            central_atoms=[central], vertex_atoms=vertex_atoms, nn=6)
        self.assertEqual(len(polyhedra), 1)
        self.assertEqual(polyhedra[0].coordination_number, 6)
        self.assertEqual(sorted(polyhedra[0].vertex_indices), [1, 2, 3, 4, 5, 6])

    def test_each_central_atom_gets_n_neighbours(self):
        lattice = Lattice.cubic(10.0)
        c1 = Atom(0, PeriodicSite('Ti', [0.2, 0.5, 0.5], lattice))
        c2 = Atom(1, PeriodicSite('Ti', [0.8, 0.5, 0.5], lattice))
        vertex_coords = [
            [0.15, 0.5, 0.5], [0.25, 0.5, 0.5],
            [0.2, 0.55, 0.5], [0.2, 0.45, 0.5],
            [0.75, 0.5, 0.5], [0.85, 0.5, 0.5],
            [0.8, 0.55, 0.5], [0.8, 0.45, 0.5],
        ]
        vertex_atoms = [Atom(i + 2, PeriodicSite('O', c, lattice))
                        for i, c in enumerate(vertex_coords)]
        polyhedra = polyhedra_from_nearest_neighbours(
            central_atoms=[c1, c2], vertex_atoms=vertex_atoms, nn=4)
        self.assertEqual(len(polyhedra), 2)
        for p in polyhedra:
            self.assertEqual(p.coordination_number, 4)


class TestPolyhedraFromClosestCentre(unittest.TestCase):

    def test_each_vertex_assigned_to_nearest_centre(self):
        lattice = Lattice.cubic(10.0)
        c1 = Atom(0, PeriodicSite('Ti', [0.2, 0.5, 0.5], lattice))
        c2 = Atom(1, PeriodicSite('Ti', [0.8, 0.5, 0.5], lattice))
        v1 = Atom(2, PeriodicSite('O', [0.25, 0.5, 0.5], lattice))
        v2 = Atom(3, PeriodicSite('O', [0.15, 0.5, 0.5], lattice))
        v3 = Atom(4, PeriodicSite('O', [0.75, 0.5, 0.5], lattice))
        v4 = Atom(5, PeriodicSite('O', [0.85, 0.5, 0.5], lattice))
        polyhedra = polyhedra_from_closest_centre(
            central_atoms=[c1, c2], vertex_atoms=[v1, v2, v3, v4])
        self.assertEqual(len(polyhedra), 2)
        self.assertEqual(sorted(polyhedra[0].vertex_indices), [2, 3])
        self.assertEqual(sorted(polyhedra[1].vertex_indices), [4, 5])

    def test_no_vertex_assigned_to_multiple_polyhedra(self):
        lattice = Lattice.cubic(10.0)
        c1 = Atom(0, PeriodicSite('Ti', [0.2, 0.5, 0.5], lattice))
        c2 = Atom(1, PeriodicSite('Ti', [0.8, 0.5, 0.5], lattice))
        vertex_atoms = [
            Atom(2, PeriodicSite('O', [0.25, 0.5, 0.5], lattice)),
            Atom(3, PeriodicSite('O', [0.15, 0.5, 0.5], lattice)),
            Atom(4, PeriodicSite('O', [0.75, 0.5, 0.5], lattice)),
            Atom(5, PeriodicSite('O', [0.85, 0.5, 0.5], lattice)),
        ]
        polyhedra = polyhedra_from_closest_centre(
            central_atoms=[c1, c2], vertex_atoms=vertex_atoms)
        all_vertex_indices = []
        for p in polyhedra:
            all_vertex_indices.extend(p.vertex_indices)
        self.assertEqual(len(all_vertex_indices), len(set(all_vertex_indices)))


class TestPolyhedraFromAtomIndices(unittest.TestCase):

    def test_correct_assignment_from_explicit_indices(self):
        lattice = Lattice.cubic(10.0)
        c1 = Atom(0, PeriodicSite('Ti', [0.2, 0.5, 0.5], lattice))
        c2 = Atom(1, PeriodicSite('Ti', [0.8, 0.5, 0.5], lattice))
        v1 = Atom(2, PeriodicSite('O', [0.25, 0.5, 0.5], lattice))
        v2 = Atom(3, PeriodicSite('O', [0.15, 0.5, 0.5], lattice))
        v3 = Atom(4, PeriodicSite('O', [0.75, 0.5, 0.5], lattice))
        v4 = Atom(5, PeriodicSite('O', [0.85, 0.5, 0.5], lattice))
        polyhedra = polyhedra_from_atom_indices(
            central_atoms=[c1, c2],
            vertex_atoms=[v1, v2, v3, v4],
            central_indices=[0, 1],
            vertex_indices=[[2, 3], [4, 5]])
        self.assertEqual(len(polyhedra), 2)
        self.assertEqual(sorted(polyhedra[0].vertex_indices), [2, 3])
        self.assertEqual(sorted(polyhedra[1].vertex_indices), [4, 5])

    def test_raises_for_mismatched_index_lengths(self):
        lattice = Lattice.cubic(10.0)
        c1 = Atom(0, PeriodicSite('Ti', [0.5, 0.5, 0.5], lattice))
        v1 = Atom(1, PeriodicSite('O', [0.6, 0.5, 0.5], lattice))
        with self.assertRaises(ValueError):
            polyhedra_from_atom_indices(
                central_atoms=[c1],
                vertex_atoms=[v1],
                central_indices=[0, 1],
                vertex_indices=[[1]])


class TestFindPolyhedraNearestNeighbours(unittest.TestCase):

    def test_end_to_end_nearest_neighbours(self):
        lattice = Lattice.cubic(10.0)
        species = ['Ti'] + ['O'] * 8
        coords = [
            [0.5, 0.5, 0.5],
            [0.6, 0.5, 0.5], [0.4, 0.5, 0.5],
            [0.5, 0.6, 0.5], [0.5, 0.4, 0.5],
            [0.5, 0.5, 0.6], [0.5, 0.5, 0.4],
            [0.8, 0.5, 0.5], [0.2, 0.5, 0.5],
        ]
        structure = Structure(lattice, species, coords)
        atoms = [Atom(i, site) for i, site in enumerate(structure.sites)]
        recipe = PolyhedraRecipe(
            method='nearest neighbours',
            central_atoms='Ti',
            vertex_atoms='O',
            n_neighbours=6)
        polyhedra = recipe.find_polyhedra(atoms, structure)
        self.assertEqual(len(polyhedra), 1)
        self.assertEqual(polyhedra[0].coordination_number, 6)
        self.assertEqual(sorted(polyhedra[0].vertex_indices), [1, 2, 3, 4, 5, 6])


class TestFindPolyhedraClosestCentre(unittest.TestCase):

    def test_end_to_end_closest_centre(self):
        lattice = Lattice.cubic(10.0)
        species = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
        coords = [
            [0.2, 0.5, 0.5], [0.8, 0.5, 0.5],
            [0.25, 0.5, 0.5], [0.15, 0.5, 0.5],
            [0.75, 0.5, 0.5], [0.85, 0.5, 0.5],
        ]
        structure = Structure(lattice, species, coords)
        atoms = [Atom(i, site) for i, site in enumerate(structure.sites)]
        recipe = PolyhedraRecipe(
            method='closest centre',
            central_atoms='Ti',
            vertex_atoms='O')
        polyhedra = recipe.find_polyhedra(atoms, structure)
        self.assertEqual(len(polyhedra), 2)
        self.assertEqual(sorted(polyhedra[0].vertex_indices), [2, 3])
        self.assertEqual(sorted(polyhedra[1].vertex_indices), [4, 5])


if __name__ == '__main__':
    unittest.main()

