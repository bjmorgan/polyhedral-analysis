import unittest
from unittest.mock import Mock, patch
from polyhedral_analysis.octahedral_analysis import check_octahedra
from polyhedral_analysis.octahedral_analysis import opposite_vertex_pairs
from polyhedral_analysis.octahedral_analysis import opposite_vertex_distances
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom

class TestOctahedralAnalysis(unittest.TestCase):

    def test_check_octahedra_passes_if_an_octahedron(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.best_fit_geometry = {'geometry': 'Octahedron'}
        check_octahedra(mock_polyhedron)

    def test_check_octahedra_raises_ValueError_if_not_an_octahedron(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
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
        mock_polyhedron.edge_graph = {i: n_list[0] for i, n_list in enumerate(neighbour_lists)}
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
   
    def test_opposite_vertex_pairs_raises_AssertionError_if_not_octahedral(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        neighbour_lists = [{0: [1, 2, 3, 4]},
                           {0: [0, 2, 4, 5]},
                           {0: [0, 1, 3, 5]},
                           {0: [0, 3, 1]},
                           {0: [0, 4, 2, 5]},
                           {0: [1, 2, 4]}]
        mock_polyhedron.edge_graph = {i: n_list[0] for i, n_list in enumerate(neighbour_lists)}
        with self.assertRaises(AssertionError):
            opposite_vertex_pairs(mock_polyhedron, check=False)

    def test_opposite_vertex_distances(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                    Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        vertices[0].distance = Mock(return_value = 1.0)
        vertices[2].distance = Mock(return_value = 2.0)
        vertices[4].distance = Mock(return_value = 3.0)
        vertex_pairs = ((vertices[0], vertices[1]),
                        (vertices[2], vertices[3]),
                        (vertices[4], vertices[5]))
        with patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
            mock_opposite_vertex_pairs.return_value = vertex_pairs
            self.assertEqual(opposite_vertex_distances(mock_polyhedron), (1.0, 2.0, 3.0)) 

if __name__ == '__main__':
    unittest.main()


