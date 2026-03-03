import unittest
from unittest.mock import Mock

import numpy as np

import polyhedral_analysis.utils as utils


class TestFlatten(unittest.TestCase):

    def test_flatten(self):
        nested = [[1, 2, 3], [4, 5, 6]]
        self.assertEqual(utils.flatten(nested), [1, 2, 3, 4, 5, 6])


class TestPruneNeighbourList(unittest.TestCase):

    def test_prunes_keys_not_in_indices(self):
        neighbours = {1: (2, 3), 2: (1, 3), 4: (1,)}
        result = utils.prune_neighbour_list(neighbours, [1, 2])
        self.assertNotIn(4, result)

    def test_prunes_values_not_in_indices(self):
        neighbours = {1: (2, 3, 5), 2: (1, 3)}
        result = utils.prune_neighbour_list(neighbours, [1, 2])
        self.assertNotIn(3, result[1])
        self.assertNotIn(5, result[1])

    def test_retains_valid_neighbours(self):
        neighbours = {1: (2, 3), 2: (1, 3), 3: (1, 2)}
        result = utils.prune_neighbour_list(neighbours, [1, 2, 3])
        self.assertEqual(set(result[1]), {2, 3})
        self.assertEqual(set(result[2]), {1, 3})
        self.assertEqual(set(result[3]), {1, 2})

    def test_empty_indices_returns_empty_dict(self):
        neighbours = {1: (2, 3), 2: (1, 3)}
        result = utils.prune_neighbour_list(neighbours, [])
        self.assertEqual(result, {})


class TestLatticeMcString(unittest.TestCase):

    def test_format(self):
        poly = Mock()
        poly.index = 5
        poly.central_atom.coords = np.array([1.5, 2.5, 3.5])
        poly.label = 'oct'
        neighbour_list = {5: (7, 12)}
        result = utils.lattice_mc_string(poly, neighbour_list)
        lines = result.strip().split('\n')
        self.assertEqual(lines[0], 'site: 5')
        self.assertTrue(lines[1].startswith('centre: '))
        self.assertIn('1.50000000', lines[1])
        self.assertEqual(lines[2], 'neighbours: 7 12')
        self.assertEqual(lines[3], 'label: oct')


if __name__ == '__main__':
    unittest.main()
