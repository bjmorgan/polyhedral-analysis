import unittest
from unittest.mock import Mock, patch
from polyhedral_analysis.octahedral_analysis import check_octahedra, opposite_vertex_pairs
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom

class TestOctahedralAnalysis(unittest.TestCase):

    def test_check_octahedra_passes_if_an_octahedron(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.best_fit_geometry = {'geometry': 'Octahedron'}
        check_octahedra(mock_polyhedron)

    def test_check_octahedra_raises_ValueError_if_not_an_octahedron(self):
        mock_polyhedron = Mock( spec=CoordinationPolyhedron )
        mock_polyhedron.best_fit_geometry = {'geometry': 'Tetrahedron'}
        with self.assertRaises(ValueError):
            check_octahedra(mock_polyhedron)

    def test_opposite_vertex_pairs(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.index = 0
        mock_vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                         Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        neighbour_lists = [{0: [1, 2, 3, 4]},
                           {0: [0, 2, 4, 5]},
                           {0: [0, 1, 3, 5]},
                           {0: [0, 3, 1, 5]},
                           {0: [0, 4, 2, 5]},
                           {0: [1, 2, 3, 4]}]
        mock_vertices_returns = []
        for n_list in neighbour_lists:
            mock_vertices_returns.append([mock_vertices[i] for i in n_list[0]])
        mock_polyhedron.vertices_by_indices = Mock(side_effect=mock_vertices_returns)
        for i, v in enumerate(mock_vertices):
            v.index = i
            v.neighbours = neighbour_lists[i]
        mock_polyhedron.vertices = mock_vertices
        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra:
            mock_check_octahedra.return_value = True
            output = opposite_vertex_pairs(mock_polyhedron)
            expected_output = ((mock_vertices[0], mock_vertices[5]),
                               (mock_vertices[1], mock_vertices[3]),
                               (mock_vertices[2], mock_vertices[4]))
            self.assertEqual(output, expected_output)
          
    
if __name__ == '__main__':
    unittest.main()


