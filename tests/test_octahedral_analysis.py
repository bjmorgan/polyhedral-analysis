import unittest
from unittest.mock import Mock, patch
from polyhedral_analysis.octahedral_analysis import check_octahedra
from polyhedral_analysis.octahedral_analysis import opposite_vertex_pairs
from polyhedral_analysis.octahedral_analysis import opposite_vertex_distances
from polyhedral_analysis.octahedral_analysis import trans_vertex_vectors
from polyhedral_analysis.octahedral_analysis import trans_vector_orthogonality
from polyhedral_analysis.octahedral_analysis import isomer_is_trans, isomer_is_cis, isomer_is_fac, isomer_is_mer
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from collections import Counter
from itertools import permutations
import numpy as np

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
        for i, mock_vertex in zip([10, 11, 12, 13, 14, 15], mock_vertices):
            mock_vertex.index = i
        neighbour_lists = [{0: [11, 12, 13, 14]},
                           {0: [10, 12, 14, 15]},
                           {0: [10, 11, 13, 15]},
                           {0: [10, 13, 11, 15]},
                           {0: [10, 14, 12, 15]},
                           {0: [11, 12, 13, 14]}]
        mock_vertices_returns = []
        for n_list in neighbour_lists:
            mock_vertices_returns.append([next(v for v in mock_vertices if v.index == i) for i in n_list[0]])
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

    def test_isomer_is_trans_raises_ValueError_if_not_2_plus_4_coordinate(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.vertex_count = Counter(['A', 'A', 'A', 'B', 'B', 'B'])
        with self.assertRaises(ValueError):
            isomer_is_trans(mock_polyhedron)

    def test_isomer_is_trans_return_True_for_trans_coordination(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        vertex_labels = ['A', 'A', 'A', 'A', 'B', 'B']
        mock_polyhedron.vertex_count = Counter(vertex_labels)
        vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                    Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        for v, l in zip(vertices, vertex_labels):
            v.label = l
        vertex_pairs = ((vertices[0], vertices[1]),
                        (vertices[2], vertices[3]),
                        (vertices[4], vertices[5]))
        with patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
            mock_opposite_vertex_pairs.return_value = vertex_pairs
            self.assertEqual(isomer_is_trans(mock_polyhedron), True)

    def test_isomer_is_trans_return_False_for_cisomer_is_coordination(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        vertex_labels = ['A', 'A', 'A', 'B', 'A', 'B']
        mock_polyhedron.vertex_count = Counter(vertex_labels)
        vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                    Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        for v, l in zip(vertices, vertex_labels):
            v.label = l
        vertex_pairs = ((vertices[0], vertices[1]),
                        (vertices[2], vertices[3]),
                        (vertices[4], vertices[5]))
        with patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
            mock_opposite_vertex_pairs.return_value = vertex_pairs
            self.assertEqual(isomer_is_trans(mock_polyhedron), False)

    def test_isomer_is_cisomer_is_returns_not_isomer_is_trans(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        with patch('polyhedral_analysis.octahedral_analysis.isomer_is_trans') as mock_isomer_is_trans:
            mock_isomer_is_trans.return_value = True
            self.assertEqual(isomer_is_cis(mock_polyhedron), False)
            mock_isomer_is_trans.return_value = False
            self.assertEqual(isomer_is_cis(mock_polyhedron), True)

    def test_isomer_is_fac_raises_ValueError_if_not_3_plus_3_coordinate(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.vertex_count = Counter(['A', 'A', 'A', 'A', 'B', 'B'])
        with self.assertRaises(ValueError):
            isomer_is_fac(mock_polyhedron)

    def test_isomer_is_fac_return_True_for_fac_coordination(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        vertex_labels = ['A', 'B', 'A', 'B', 'A', 'B']
        mock_polyhedron.vertex_count = Counter(vertex_labels)
        vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                    Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        for v, l in zip(vertices, vertex_labels):
            v.label = l
        vertex_pairs = ((vertices[0], vertices[1]),
                        (vertices[2], vertices[3]),
                        (vertices[4], vertices[5]))
        with patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
            mock_opposite_vertex_pairs.return_value = vertex_pairs
            self.assertEqual(isomer_is_fac(mock_polyhedron), True)

    def test_isomer_is_fac_return_False_for_mer_coordination(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        vertex_labels = ['A', 'A', 'B', 'A', 'B', 'B']
        mock_polyhedron.vertex_count = Counter(vertex_labels)
        vertices = [Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom),
                    Mock(spec=Atom), Mock(spec=Atom), Mock(spec=Atom)]
        for v, l in zip(vertices, vertex_labels):
            v.label = l
        vertex_pairs = ((vertices[0], vertices[1]),
                        (vertices[2], vertices[3]),
                        (vertices[4], vertices[5]))
        with patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
            mock_opposite_vertex_pairs.return_value = vertex_pairs
            self.assertEqual(isomer_is_fac(mock_polyhedron), False)

    def test_isomer_is_mer_returns_not_isomer_is_fac(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        with patch('polyhedral_analysis.octahedral_analysis.isomer_is_fac') as mock_isomer_is_fac:
            mock_isomer_is_fac.return_value = True
            self.assertEqual(isomer_is_mer(mock_polyhedron), False)
            mock_isomer_is_fac.return_value = False
            self.assertEqual(isomer_is_mer(mock_polyhedron), True)

    def test_trans_vertex_vectors(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        mock_polyhedron.vertex_vectors = Mock(return_value=np.array([
            [1, 0, 0],   # 0
            [0, 1, 0],   # 1
            [0, 0, 1],   # 2
            [-1, 0, 0],  # 3
            [0, -1, 0],  # 4
            [0, 0, -1]   # 5
        ]))
        
        # Test cases with different orderings of trans pairs
        test_cases = [
            # Case 1: Sequential order
            {
                'vertex_pairs': [(Mock(index=0), Mock(index=3)),
                                 (Mock(index=1), Mock(index=4)),
                                 (Mock(index=2), Mock(index=5))],
                'expected_vectors': [np.array([-2, 0, 0]),
                                     np.array([0, -2, 0]),
                                     np.array([0, 0, -2])]
            },
            # Case 2: Non-sequential order
            {
                'vertex_pairs': [(Mock(index=2), Mock(index=5)),
                                 (Mock(index=0), Mock(index=3)),
                                 (Mock(index=1), Mock(index=4))],
                'expected_vectors': [np.array([0, 0, -2]),
                                     np.array([-2, 0, 0]),
                                     np.array([0, -2, 0])]
            },
            # Case 3: Reversed pairs
            {
                'vertex_pairs': [(Mock(index=3), Mock(index=0)),
                                 (Mock(index=4), Mock(index=1)),
                                 (Mock(index=5), Mock(index=2))],
                'expected_vectors': [np.array([2, 0, 0]),
                                     np.array([0, 2, 0]),
                                     np.array([0, 0, 2])]
            },
            # Case 4: Mixed order and direction
            {
                'vertex_pairs': [(Mock(index=1), Mock(index=4)),
                                 (Mock(index=5), Mock(index=2)),
                                 (Mock(index=3), Mock(index=0))],
                'expected_vectors': [np.array([0, -2, 0]),
                                     np.array([0, 0, 2]),
                                     np.array([2, 0, 0])]
            }
        ]
        
        for case in test_cases:
            with self.subTest(vertex_pairs=case['vertex_pairs']):
                with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra, \
                     patch('polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs') as mock_opposite_vertex_pairs:
                    
                    mock_check_octahedra.return_value = None
                    mock_opposite_vertex_pairs.return_value = case['vertex_pairs']
                    
                    mock_polyhedron.vertex_internal_index_from_global_index = lambda x: x
                    
                    vectors = trans_vertex_vectors(mock_polyhedron)
                    
                    for v, e in zip(vectors, case['expected_vectors']):
                        np.testing.assert_array_almost_equal(v, e)
                    
                    mock_check_octahedra.assert_called_once()
                    mock_opposite_vertex_pairs.assert_called_once()

    def test_trans_vertex_vectors_not_octahedron(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        
        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra:
            mock_check_octahedra.side_effect = ValueError("Not an octahedron")
            
            with self.assertRaises(ValueError):
                trans_vertex_vectors(mock_polyhedron)

    def test_trans_vector_orthogonality(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        
        # Mock trans vectors for a perfect octahedron
        mock_trans_vectors = [
            np.array([2, 0, 0]),
            np.array([0, 2, 0]),
            np.array([0, 0, 2])
        ]
        
        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra, \
             patch('polyhedral_analysis.octahedral_analysis.trans_vertex_vectors') as mock_trans_vertex_vectors:
            
            mock_check_octahedra.return_value = None
            mock_trans_vertex_vectors.return_value = mock_trans_vectors
            
            angles = trans_vector_orthogonality(mock_polyhedron)
            
            self.assertEqual(len(angles), 3)
            for angle in angles:
                self.assertAlmostEqual(angle, 0.0, places=6)
            
            mock_check_octahedra.assert_called_once()
            mock_trans_vertex_vectors.assert_called_once()

    def test_trans_vector_orthogonality_distorted(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        
        # Mock trans vectors for a slightly distorted octahedron
        mock_trans_vectors = [
            np.array([2, 0, 0]),
            np.array([0, 2, 0.1]),
            np.array([0.1, 0.1, 2])  # Slightly off from perfect
        ]
        
        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra, \
             patch('polyhedral_analysis.octahedral_analysis.trans_vertex_vectors') as mock_trans_vertex_vectors:
            
            mock_check_octahedra.return_value = None
            mock_trans_vertex_vectors.return_value = mock_trans_vectors
            
            angles = trans_vector_orthogonality(mock_polyhedron)
            expected_angles = (2.8731459157778483, 5.724810452223508, 6.396811661684919)

            np.testing.assert_allclose(sorted(angles), sorted(expected_angles), rtol=1e-6)

            mock_check_octahedra.assert_called_once()
            mock_trans_vertex_vectors.assert_called_once()

    def test_trans_vector_orthogonality_not_octahedron(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        
        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra:
            mock_check_octahedra.side_effect = ValueError("Not an octahedron")
            
            with self.assertRaises(ValueError):
                trans_vector_orthogonality(mock_polyhedron)

    def test_trans_vector_orthogonality_permutation_invariance(self):
        mock_polyhedron = Mock(spec=CoordinationPolyhedron)
        
        # Define mock trans vectors
        v1 = np.array([2, 0, 0])
        v2 = np.array([0, 2, 0.1])
        v3 = np.array([0.1, 0.1, 2])  # Slightly distorted
        
        base_vectors = [v1, v2, v3]
        expected_angles = (2.8731459157778483, 5.724810452223508, 6.396811661684919)

        with patch('polyhedral_analysis.octahedral_analysis.check_octahedra') as mock_check_octahedra, \
             patch('polyhedral_analysis.octahedral_analysis.trans_vertex_vectors') as mock_trans_vertex_vectors:

            mock_check_octahedra.return_value = None

            # Test all permutations of the vectors
            all_angles = []
            for perm in permutations(base_vectors):
                mock_trans_vertex_vectors.return_value = list(perm)
                angles = trans_vector_orthogonality(mock_polyhedron)
                all_angles.append(tuple(sorted(angles)))  # Store as a sorted tuple

            # Check that all permutations give the same angles
            self.assertEqual(len(set(all_angles)), 1)

            # Check that the angles match the expected values
            for angles in all_angles:
                self.assertEqual(len(angles), 3)  # Should always have 3 angles
                np.testing.assert_allclose(sorted(angles), expected_angles, rtol=1e-6)

            # Check that mock methods were called the correct number of times
            self.assertEqual(mock_check_octahedra.call_count, 6)  # Once for each permutation
            self.assertEqual(mock_trans_vertex_vectors.call_count, 6)  # Once for each permutation

if __name__ == '__main__':
    unittest.main()


