import unittest
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from pymatgen.core import Site
from unittest.mock import Mock, patch, PropertyMock, call
import copy
import numpy as np
from scipy.spatial import ConvexHull

module_str = 'polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron'


def mock_atom_lt(self, other):
    return self.index < other.index


def mock_atom_eq(self, other):
    return self.index == other.index

class TestCoordinationPolyhedronInit(unittest.TestCase):

    def test_coordination_polyhedron_is_initialised(self):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 10
        mock_central_atom.label = 'Li'
        mock_vertices = [Mock(spec=Atom) for i in range(6)]
        for i, v in enumerate(mock_vertices, 1):
            v.neighbours = None
            v.__lt__ = mock_atom_lt
            v.index = i
            v.in_polyhedra = []
        CoordinationPolyhedron(central_atom=mock_central_atom,
                               vertices=mock_vertices)

class TestCoordinationPolyhedron(unittest.TestCase):

    def setUp(self):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 0
        mock_central_atom.label = 'A'
        mock_central_atom.__eq__ = mock_atom_eq
        mock_vertices = [Mock(spec=Atom) for i in range(6)]
        for i, v in enumerate(mock_vertices, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.__eq__ = mock_atom_eq
            v.in_polyhedra = []
        self.coordination_polyhedron = CoordinationPolyhedron(
            central_atom=mock_central_atom,
            vertices=mock_vertices)

    def test_equal_members(self):
        other_coordination_polyhedron = copy.deepcopy(
            self.coordination_polyhedron)
        other_coordination_polyhedron.vertices[0].neighbours = {0: [1, 2, 3]}
        other_coordination_polyhedron.vertices[4].neighbours = {4: [1, 3, 5]}
        self.assertTrue(self.coordination_polyhedron.equal_members(
            other_coordination_polyhedron))

    def test_vertex_vectors(self):
        mock_min_image_coords = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 2.0]
        ])
        self.coordination_polyhedron.minimum_image_vertex_coordinates = Mock(return_value=mock_min_image_coords)

        mock_centroid = np.array([-0.1, 0.1, 0.2])
        self.coordination_polyhedron.centroid = Mock(return_value=mock_centroid)

        mock_central_atom_coords = np.array([0.1, -0.1, -0.2])
        self.coordination_polyhedron.central_atom.coords = mock_central_atom_coords

        expected_vectors_centroid = mock_min_image_coords - mock_centroid
        actual_vectors_centroid = self.coordination_polyhedron.vertex_vectors(reference='centroid')
        np.testing.assert_array_almost_equal(actual_vectors_centroid, expected_vectors_centroid)

        expected_vectors_central_atom = mock_min_image_coords - mock_central_atom_coords
        actual_vectors_central_atom = self.coordination_polyhedron.vertex_vectors(reference='central_atom')
        np.testing.assert_array_almost_equal(actual_vectors_central_atom, expected_vectors_central_atom)

        self.coordination_polyhedron.centroid.assert_called_once()
        self.assertEqual(self.coordination_polyhedron.minimum_image_vertex_coordinates.call_count, 2)

    def test_vertex_vectors_invalid_reference(self):
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.vertex_vectors(reference="invalid")


    def test_angles(self):
        vertex_vectors = np.array([[1.0, 0.0, 0.0],
                                   [0.0, 1.0, 0.0],
                                   [0.0, -1.0, 0.0]])
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors') as mock_vertex_vectors:
            mock_vertex_vectors.return_value = vertex_vectors
            angles = self.coordination_polyhedron.angles()
            np.testing.assert_equal(angles, [90.0, 90.0, 180.0])

    def test_vertex_distances(self):
        mock_vectors = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3], [3, 4, 0]])
        expected_distances = np.array([1.0, 2.0, 3.0, 5.0])
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vectors):
            # Test with default reference (central_atom)
            distances = self.coordination_polyhedron.vertex_distances()
            print(distances)
            np.testing.assert_equal(distances, expected_distances)
            # Test with centroid reference
            distances_centroid = self.coordination_polyhedron.vertex_distances(reference='centroid')
            np.testing.assert_equal(distances_centroid, expected_distances)
        # Test with invalid reference
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.vertex_distances(reference='invalid')

    def test_vertex_distances_and_labels(self):
        # Mock the vertex_distances method and vertex_labels property
        mock_distances = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        mock_labels = ['A', 'B', 'C', 'D', 'E', 'F']
        
        self.coordination_polyhedron.vertex_distances = Mock(return_value=mock_distances)
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_labels', new_callable=PropertyMock) as mock_vertex_labels:
            mock_vertex_labels.return_value = mock_labels

            # Test with default reference (central_atom)
            result = self.coordination_polyhedron.vertex_distances_and_labels()
            expected = ((1.0, 'A'), (2.0, 'B'), (3.0, 'C'), (4.0, 'D'), (5.0, 'E'), (6.0, 'F'))
            self.assertEqual(result, expected)
            self.coordination_polyhedron.vertex_distances.assert_called_with(reference='central_atom')

            # Test with centroid reference
            result_centroid = self.coordination_polyhedron.vertex_distances_and_labels(reference='centroid')
            self.assertEqual(result_centroid, expected)
            self.coordination_polyhedron.vertex_distances.assert_called_with(reference='centroid')

            # Test with invalid reference
            with self.assertRaises(ValueError):
                self.coordination_polyhedron.vertex_distances_and_labels(reference='invalid')

    def test_update_vertex_neighbours(self):
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.edge_graph', new_callable=PropertyMock) as mock_edge_graph:
            edge_graph = {1: [2, 3, 4, 5],
                          2: [1, 3, 5, 6],
                          3: [1, 2, 4, 6],
                          4: [1, 3, 5, 6],
                          5: [1, 2, 4, 6],
                          6: [2, 3, 4, 5]}
            mock_edge_graph.return_value = edge_graph
            self.coordination_polyhedron.update_vertex_neighbours()
            for v, n in zip(self.coordination_polyhedron.vertices, edge_graph.values()):
                self.assertEqual(
                    v._neighbours[self.coordination_polyhedron.index], n)

    def test_vertices_by_indices(self):
        vertices = self.coordination_polyhedron.vertices_by_indices([2, 4])
        self.assertEqual(sorted([v.index for v in vertices]), [2, 4])

    def test_vertex_indices(self):
        self.assertEqual(self.coordination_polyhedron.vertex_indices, [
                         1, 2, 3, 4, 5, 6])

    def test_vertex_coords(self):
        coords = np.array([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0],
                           [-1.0, 0.0, 0.0],
                           [0.0, -1.0, 0.0],
                           [0.0, 0.0, -1.0]])
        for v, c in zip(self.coordination_polyhedron.vertices, coords):
            v.coords = c
        np.testing.assert_array_equal(
            self.coordination_polyhedron.vertex_coords, coords)

    def test_vertex_labels(self):
        labels = ['A', 'B', 'C', 'D', 'E', 'F']
        for v, l in zip(self.coordination_polyhedron.vertices, labels):
            v.label = l
        self.assertEqual(self.coordination_polyhedron.vertex_labels, labels)

    def test_coordination_number(self):
        self.assertEqual(self.coordination_polyhedron.coordination_number,
                         len(self.coordination_polyhedron.vertices))

    def test_index(self):
        self.assertEqual(self.coordination_polyhedron.index,
                         self.coordination_polyhedron.central_atom.index)

    def test_edge_vertex_indices(self):
        mock_edge_graph = {
            1: [2, 3, 4],
            2: [1, 3, 4],
            3: [1, 2, 4],
            4: [1, 2, 3]
        }
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.edge_graph', new_callable=PropertyMock) as mock_edge_graph_property:
            mock_edge_graph_property.return_value = mock_edge_graph
            expected_edges = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
            self.assertEqual(self.coordination_polyhedron.edge_vertex_indices(), expected_edges)

    def test_edge_graph_if_cached(self):
        edge_graph = {1: [2, 3, 4, 5],
                      2: [1, 3, 5, 6],
                      3: [1, 2, 4, 6],
                      4: [1, 3, 5, 6],
                      5: [1, 2, 4, 6],
                      6: [2, 3, 4, 5]}
        self.coordination_polyhedron._edge_graph = edge_graph
        self.assertEqual(self.coordination_polyhedron.edge_graph, edge_graph)

    def test_edge_graph_if_not_cached(self):
        edge_graph = {1: [2, 3, 4, 5],
                      2: [1, 3, 5, 6],
                      3: [1, 2, 4, 6],
                      4: [1, 3, 5, 6],
                      5: [1, 2, 4, 6],
                      6: [2, 3, 4, 5]}
        self.coordination_polyhedron.construct_edge_graph = Mock(
            return_value=edge_graph)
        self.coordination_polyhedron.update_vertex_neighbours = Mock()
        self.assertEqual(self.coordination_polyhedron.edge_graph, edge_graph)
        self.assertEqual(
            self.coordination_polyhedron.update_vertex_neighbours.call_count, 1)
        self.assertEqual(
            self.coordination_polyhedron.construct_edge_graph.call_count, 1)

    def test_symmetry_measure_delegates_to_vertex_vectors(self):
        polyhedron = self.coordination_polyhedron
        mock_vectors = np.array([
            [0.0, 0.0, 2.0],
            [0.0, 0.0, -2.0],
            [2.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, -2.0, 0.0],
        ])
        polyhedron.vertex_vectors = Mock(return_value=mock_vectors)
        sm = polyhedron.symmetry_measure
        polyhedron.vertex_vectors.assert_called_with(reference='central_atom')
        self.assertAlmostEqual(sm['Octahedron'], 0.0)

    def test_symmetry_measure_is_cached(self):
        polyhedron = self.coordination_polyhedron
        mock_vectors = np.array([
            [0.0, 0.0, 2.0],
            [0.0, 0.0, -2.0],
            [2.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, -2.0, 0.0],
        ])
        polyhedron.vertex_vectors = Mock(return_value=mock_vectors)
        _ = polyhedron.symmetry_measure
        call_count_after_first = polyhedron.vertex_vectors.call_count
        _ = polyhedron.symmetry_measure
        self.assertEqual(polyhedron.vertex_vectors.call_count, call_count_after_first)

    def test_faces(self):
        polyhedron = self.coordination_polyhedron
        polyhedron.convex_hull = Mock(return_value=Mock(spec=ConvexHull))
        returned_vertex_indices = [1, 2, 3, 4, 5, 6]
        with patch(f'{module_str}.vertex_indices', new_callable=PropertyMock) as mock_vertex_indices:
            mock_vertex_indices.return_value = returned_vertex_indices
            simplices = [np.array([0, 3, 1]),
                         np.array([0, 3, 4]),
                         np.array([5, 3, 1]),
                         np.array([5, 3, 4]),
                         np.array([2, 0, 4]),
                         np.array([2, 0, 1]),
                         np.array([2, 5, 4]),
                         np.array([2, 5, 1])]
            with patch('polyhedral_analysis.coordination_polyhedron.merge_coplanar_simplices') as mock_merge_coplanar_simplices:
                mock_merge_coplanar_simplices.return_value = simplices
                faces = polyhedron.faces()
                self.assertEqual(faces, ((1, 2, 4), (1, 4, 5), (2, 4, 6), (4, 5, 6),
                                         (1, 3, 5), (1, 2, 3), (3, 5, 6), (2, 3, 6)))

    def test_neighbours(self):
        mock_central_atom_i = Mock(spec=Atom)
        mock_central_atom_i.in_polyhedra = []
        mock_central_atom_i.index = 3
        mock_central_atom_i.label = 'Li'
        mock_central_atom_j = Mock(spec=Atom)
        mock_central_atom_j.in_polyhedra = []
        mock_central_atom_j.index = 7
        mock_central_atom_j.label = 'Li'
        mock_vertices = [Mock(spec=Atom) for i in range(10)]
        for i, v in enumerate(mock_vertices, 1):
            v.neighbours = None
            v.__lt__ = mock_atom_lt
            v.index = i
            v.in_polyhedra = []
        polyhedron_i = CoordinationPolyhedron(central_atom=mock_central_atom_i,
                                              vertices=mock_vertices[:6])
        polyhedron_j = CoordinationPolyhedron(central_atom=mock_central_atom_j,
                                              vertices=mock_vertices[-6:])
        polyhedron_i._edge_graph = Mock()
        polyhedron_j._edge_graph = Mock()
        self.assertEqual(polyhedron_i.neighbours(), (polyhedron_j,))
        self.assertEqual(polyhedron_j.neighbours(), (polyhedron_i,))

    def test_interection(self):
        mock_central_atom_i = Mock(spec=Atom)
        mock_central_atom_i.in_polyhedra = []
        mock_central_atom_i.index = 3
        mock_central_atom_i.label = 'Li'
        mock_central_atom_j = Mock(spec=Atom)
        mock_central_atom_j.in_polyhedra = []
        mock_central_atom_j.index = 7
        mock_central_atom_j.label = 'Li'
        mock_vertices = [Mock(spec=Atom) for i in range(10)]
        for i, v in enumerate(mock_vertices, 1):
            v.neighbours = None
            v.__lt__ = mock_atom_lt
            v.index = i
            v.in_polyhedra = []
        polyhedron_i = CoordinationPolyhedron(central_atom=mock_central_atom_i,
                                              vertices=mock_vertices[:6])
        polyhedron_j = CoordinationPolyhedron(central_atom=mock_central_atom_j,
                                              vertices=mock_vertices[-6:])
        self.assertEqual(polyhedron_i.intersection(polyhedron_j), (5, 6))
        self.assertEqual(polyhedron_j.intersection(polyhedron_i), (5, 6))

    def test_neighbours_by_index_and_shared_vertices(self):
        polyhedron_i = self.coordination_polyhedron
        polyhedron_j = copy.deepcopy(polyhedron_i)
        polyhedron_i.neighbours = Mock(return_value=(polyhedron_j,))
        polyhedron_j.neighbours = Mock(return_value=(polyhedron_i,))
        polyhedron_i.intersection = Mock(return_value=(5, 6))
        polyhedron_j.intersection = Mock(return_value=(5, 6))
        self.assertEqual(polyhedron_i.neighbours_by_index_and_shared_vertices(), 
                         {polyhedron_j.index: (5, 6)})
        self.assertEqual(polyhedron_j.neighbours_by_index_and_shared_vertices(), 
                         {polyhedron_i.index: (5, 6)})

    def test_vertex_internal_index_from_global_index(self):
        # Set up mock vertex indices
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_indices', new_callable=PropertyMock) as mock_vertex_indices:
            mock_vertex_indices.return_value = [10, 20, 30, 40, 50, 60]
        
            # Test successful cases
            self.assertEqual(self.coordination_polyhedron.vertex_internal_index_from_global_index(10), 0)
            self.assertEqual(self.coordination_polyhedron.vertex_internal_index_from_global_index(30), 2)
            self.assertEqual(self.coordination_polyhedron.vertex_internal_index_from_global_index(60), 5)
        
            # Test case where global index doesn't exist
            with self.assertRaises(ValueError):
                self.coordination_polyhedron.vertex_internal_index_from_global_index(70)
        
            # Verify that vertex_indices was called
            mock_vertex_indices.assert_called()

    def test_centroid(self):
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.minimum_image_vertex_coordinates') as mock_min_image:
            mock_min_image.return_value = np.array([
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0]
            ])
            centroid = self.coordination_polyhedron.centroid()
            # The centroid should be (0, 0, 0)
            np.testing.assert_array_almost_equal(centroid, np.array([0.0, 0.0, 0.0]))

    def test_centroid_with_offset(self):
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.minimum_image_vertex_coordinates') as mock_min_image:
            mock_min_image.return_value = np.array([
                [2.0, 1.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 2.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 2.0],
                [1.0, 1.0, 0.0]
            ])
            centroid = self.coordination_polyhedron.centroid()
            # The centroid should be (1, 1, 1)
            np.testing.assert_array_almost_equal(centroid, np.array([1.0, 1.0, 1.0]))

    def test_centroid_to_central_atom_vector(self):
        self.coordination_polyhedron.central_atom.coords = np.array([0.0, 0.0, 0.0])
        with patch.object(self.coordination_polyhedron, 'centroid') as mock_centroid:
            mock_centroid.return_value = np.array([0.0, 0.0, 0.0])
            displacement = self.coordination_polyhedron.centroid_to_central_atom_vector()
            # The centroid is at (0, 0, 0), so the displacement should be (0, 0, 0)
            np.testing.assert_array_almost_equal(displacement, np.array([0.0, 0.0, 0.0]))

    def test_centroid_to_central_atom_vector_with_offset(self):
        self.coordination_polyhedron.central_atom.coords = np.array([0.1, 0.2, 0.3])
        with patch.object(self.coordination_polyhedron, 'centroid') as mock_centroid:
            mock_centroid.return_value = np.array([0.0, 0.0, 0.0])
            displacement = self.coordination_polyhedron.centroid_to_central_atom_vector()
            # The centroid is at (0, 0, 0), so the displacement should be (0.1, 0.2, 0.3)
            np.testing.assert_array_almost_equal(displacement, np.array([0.1, 0.2, 0.3]))

    def test_radial_distortion_parameter(self):
        # Mock the vertex_vectors method
        mock_distances = np.array([1.0, 2.0, 3.0, 4.0])
        self.coordination_polyhedron.vertex_distances = Mock(return_value=mock_distances)

        # Test MSD (default)
        msd_normalized = self.coordination_polyhedron.radial_distortion_parameter()
        self.assertAlmostEqual(msd_normalized, 0.2, places=6)

        # Test MAD
        mad_normalized = self.coordination_polyhedron.radial_distortion_parameter(method='MAD')
        self.assertAlmostEqual(mad_normalized, 0.4, places=6)

        # Test non-normalized MSD
        msd_non_normalized = self.coordination_polyhedron.radial_distortion_parameter(normalize=False)
        self.assertAlmostEqual(msd_non_normalized, 1.25, places=6)

        # Test non-normalized MAD
        mad_non_normalized = self.coordination_polyhedron.radial_distortion_parameter(normalize=False, method='MAD')
        self.assertAlmostEqual(mad_non_normalized, 1.0, places=6)

        # Verify that vertex_distances was called with the correct reference
        self.coordination_polyhedron.vertex_distances.assert_called_with(reference='centroid')

    def test_radial_distortion_parameter_perfect_polyhedron(self):
        # Mock a perfect polyhedron where all vertices are equidistant from the reference
        mock_vectors = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0]
        ])
        self.coordination_polyhedron.vertex_vectors = Mock(return_value=mock_vectors)

        self.assertAlmostEqual(self.coordination_polyhedron.radial_distortion_parameter(method='MAD'), 0.0)
        self.assertAlmostEqual(self.coordination_polyhedron.radial_distortion_parameter(method='MSD'), 0.0)
        self.assertAlmostEqual(self.coordination_polyhedron.radial_distortion_parameter(normalize=False, method='MAD'), 0.0)
        self.assertAlmostEqual(self.coordination_polyhedron.radial_distortion_parameter(normalize=False, method='MSD'), 0.0)

    def test_radial_distortion_parameter_invalid_inputs(self):
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.radial_distortion_parameter(reference="invalid")
        
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.radial_distortion_parameter(method="invalid")

    def test_sharing_neighbour_lists(self):
        mock_neighbours = {
            1: (1,),
            2: (1, 2),
            3: (1, 2, 3),
            4: (1, 2, 3, 4)
        }
        self.coordination_polyhedron.neighbours_by_index_and_shared_vertices = Mock(return_value=mock_neighbours)
    
        self.assertEqual(self.coordination_polyhedron.corner_sharing_neighbour_list(), (1,))
        self.assertEqual(self.coordination_polyhedron.edge_sharing_neighbour_list(), (2,))
        self.assertEqual(self.coordination_polyhedron.face_sharing_neighbour_list(), (3, 4))

    def test_vertex_vector_projections(self):
        mock_vertex_vectors = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-1, 0, 0]
        ])
    
        input_vectors = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 1]
        ])
    
        expected_projections = np.array([
            [1, 0, 1/np.sqrt(3)],
            [0, 1, 1/np.sqrt(3)],
            [0, 0, 1/np.sqrt(3)],
            [-1, 0, -1/np.sqrt(3)]
        ])
    
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vertex_vectors):
            projections = self.coordination_polyhedron.vertex_vector_projections(input_vectors)
            np.testing.assert_array_almost_equal(projections, expected_projections)

        # Test with a single vector
        single_vector = np.array([1, 1, 1])
        expected_single_projection = np.array([[1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3)]]).T
    
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vertex_vectors):
            single_projection = self.coordination_polyhedron.vertex_vector_projections(single_vector)
            np.testing.assert_array_almost_equal(single_projection, expected_single_projection)

        # Test with different reference
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vertex_vectors):
            projections_central = self.coordination_polyhedron.vertex_vector_projections(input_vectors, reference='central_atom')
            np.testing.assert_array_almost_equal(projections_central, expected_projections)

    def test_coordination_distances(self):
        mock_distances = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    
        with patch.object(self.coordination_polyhedron, 'vertex_distances', return_value=mock_distances) as mock_vertex_distances:
            distances = self.coordination_polyhedron.coordination_distances()
        
            self.assertEqual(distances, list(mock_distances))
            self.assertIsInstance(distances, list)
            mock_vertex_distances.assert_called_once_with(reference='central_atom')

    def test_from_sites(self):
        mock_central_site = Mock()
        mock_vertex_sites = [Mock() for _ in range(3)]
    
        with patch('polyhedral_analysis.coordination_polyhedron.Atom') as mock_atom:
            def create_mock_atom(index, site, label=None):
                mock = Mock(spec=Atom)
                mock.index = index
                mock.site = site
                mock.label = label if label else site.species_string
                mock.in_polyhedra = []
                mock._neighbours = {}
                mock.__lt__ = lambda self, other: self.index < other.index
                return mock
        
            mock_atom.side_effect = create_mock_atom
        
            polyhedron = CoordinationPolyhedron.from_sites(mock_central_site, mock_vertex_sites)
        
            self.assertEqual(polyhedron.central_atom.index, -1)
            self.assertEqual(polyhedron.central_atom.site, mock_central_site)
            self.assertEqual(len(polyhedron.vertices), 3)
            for i, vertex in enumerate(polyhedron.vertices):
                self.assertEqual(vertex.index, i)
                self.assertEqual(vertex.site, mock_vertex_sites[i])
                self.assertEqual(vertex.in_polyhedra, [polyhedron])
                self.assertEqual(vertex._neighbours, {})

            # Check that vertices are sorted
            self.assertEqual([v.index for v in polyhedron.vertices], [0, 1, 2])


    def test_vertex_angles_central_atom(self):
        mock_vectors = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, -1]
        ])
        self.coordination_polyhedron.vertex_vectors = Mock(return_value = mock_vectors)
        self.coordination_polyhedron.vertex_internal_index_from_global_index = Mock(side_effect = (0, 1, 1, 2, 2, 3))
        
        vertex_pairs = ((0, 1), (1, 2), (2, 3))
        expected_angles = np.array([90, 90, 180])
        
        angles = self.coordination_polyhedron.vertex_angles(vertex_pairs)
        
        np.testing.assert_array_almost_equal(angles, expected_angles, decimal=5)
        self.coordination_polyhedron.vertex_vectors.assert_called_once_with(reference='central_atom')

    def test_vertex_angles_centroid(self):
        mock_vectors = np.array([
            [1, 1, 0],
            [-1, 1, 0],
            [0, -1, 1],
            [0, -1, -1]
        ])
        self.coordination_polyhedron.vertex_vectors = Mock(return_value = mock_vectors)
        self.coordination_polyhedron.vertex_internal_index_from_global_index = Mock(side_effect = (0, 1, 1, 2, 2, 3))
        
        vertex_pairs = ((0, 1), (1, 2), (2, 3))
        expected_angles = np.array([90, 120, 90])
        
        angles = self.coordination_polyhedron.vertex_angles(vertex_pairs, reference='centroid')
        
        np.testing.assert_array_almost_equal(angles, expected_angles, decimal=5)
        self.coordination_polyhedron.vertex_vectors.assert_called_once_with(reference='centroid')

    def test_vertex_angles_invalid_index(self):
        mock_vectors = np.array([
            [1, 0, 0], 
            [0, 1, 0]]
        )
        self.coordination_polyhedron.vertex_vectors = Mock(return_value = mock_vectors)
        self.coordination_polyhedron.vertex_internal_index_from_global_index = Mock(side_effect=[0, ValueError])
        
        vertex_pairs = ((0, 2),)
        
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.vertex_angles(vertex_pairs)

    def test_off_centre_displacement(self):
        mock_displacement_vector = np.array([0.1, 0.2, 0.3])
        self.coordination_polyhedron.centroid_to_central_atom_vector = Mock(return_value = mock_displacement_vector)
        expected_displacement = np.linalg.norm(mock_displacement_vector)
        self.assertEqual(self.coordination_polyhedron.off_centre_displacement, expected_displacement)

    def test_vertex_vector_orientations(self):
        # Mock vertex vectors
        mock_vectors = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-1, -1, -1]
        ])
        
        # Calculate the actual angle for [-1, -1, -1] vector
        actual_theta = 180 - np.degrees(np.arccos(1/np.sqrt(3)))
        
        # Expected results for degrees
        expected_centroid_deg = [
            (90.0, 0.0),
            (90.0, 90.0),
            (0.0, 0.0),
            (actual_theta, -135.0)
        ]
        
        # Expected results for radians
        expected_centroid_rad = [
            (np.pi/2, 0.0),
            (np.pi/2, np.pi/2),
            (0.0, 0.0),
            (np.pi - np.arccos(1/np.sqrt(3)), -3*np.pi/4)
        ]
        
        # Test with centroid reference (default)
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vectors):
            orientations_deg = self.coordination_polyhedron.vertex_vector_orientations()
            orientations_rad = self.coordination_polyhedron.vertex_vector_orientations(units='radians')
            
            np.testing.assert_array_almost_equal(orientations_deg, expected_centroid_deg)
            np.testing.assert_array_almost_equal(orientations_rad, expected_centroid_rad)
        
        # Test with central_atom reference
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vectors):
            orientations_central = self.coordination_polyhedron.vertex_vector_orientations(reference='central_atom')
            
            np.testing.assert_array_almost_equal(orientations_central, expected_centroid_deg)
        
        # Test with return_distance=True
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', return_value=mock_vectors):
            orientations_with_dist = self.coordination_polyhedron.vertex_vector_orientations(return_distance=True)
            
            expected_with_dist = [
                (90.0, 0.0, 1.0),
                (90.0, 90.0, 1.0),
                (0.0, 0.0, 1.0),
                (actual_theta, -135.0, np.sqrt(3))
            ]
            np.testing.assert_array_almost_equal(orientations_with_dist, expected_with_dist)
        
        # Test with invalid reference
        with self.assertRaises(ValueError):
            self.coordination_polyhedron.vertex_vector_orientations(reference='invalid')

class TestCoordinationPolyhedronLabel(unittest.TestCase):

    def test_default_label_from_central_atom(self):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 0
        mock_central_atom.label = 'Ti'
        mock_vertices = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly = CoordinationPolyhedron(
            central_atom=mock_central_atom, vertices=mock_vertices)
        self.assertEqual(poly.label, 'Ti')

    def test_explicit_label_overrides_central_atom(self):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 0
        mock_central_atom.label = 'Ti'
        mock_vertices = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly = CoordinationPolyhedron(
            central_atom=mock_central_atom, vertices=mock_vertices,
            label='oct')
        self.assertEqual(poly.label, 'oct')

    def test_set_label(self):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 0
        mock_central_atom.label = 'Ti'
        mock_vertices = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly = CoordinationPolyhedron(
            central_atom=mock_central_atom, vertices=mock_vertices)
        poly.set_label('tet')
        self.assertEqual(poly.label, 'tet')


class TestCoordinationPolyhedronEquality(unittest.TestCase):

    def setUp(self):
        self.mock_central_atom_a = Mock(spec=Atom)
        self.mock_central_atom_a.in_polyhedra = []
        self.mock_central_atom_a.index = 0
        self.mock_central_atom_a.label = 'A'
        self.mock_central_atom_a.__eq__ = mock_atom_eq
        self.mock_vertices_a = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(self.mock_vertices_a, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.__eq__ = mock_atom_eq
            v.in_polyhedra = []
        self.poly_a = CoordinationPolyhedron(
            central_atom=self.mock_central_atom_a,
            vertices=self.mock_vertices_a)

    def test_equal_vertices_true_for_same_indices(self):
        mock_central_atom_b = Mock(spec=Atom)
        mock_central_atom_b.in_polyhedra = []
        mock_central_atom_b.index = 10
        mock_central_atom_b.label = 'B'
        mock_vertices_b = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices_b, 1):
            v._neighbours = {}
            v.index = i  # same indices as poly_a
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly_b = CoordinationPolyhedron(
            central_atom=mock_central_atom_b, vertices=mock_vertices_b)
        self.assertTrue(self.poly_a.equal_vertices(poly_b))

    def test_equal_vertices_false_for_different_indices(self):
        mock_central_atom_b = Mock(spec=Atom)
        mock_central_atom_b.in_polyhedra = []
        mock_central_atom_b.index = 10
        mock_central_atom_b.label = 'B'
        mock_vertices_b = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices_b, 11):
            v._neighbours = {}
            v.index = i  # different indices
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly_b = CoordinationPolyhedron(
            central_atom=mock_central_atom_b, vertices=mock_vertices_b)
        self.assertFalse(self.poly_a.equal_vertices(poly_b))

    def test_eq_delegates_to_equal_edge_graph(self):
        poly_b = copy.deepcopy(self.poly_a)
        edge_graph = {1: [2, 3], 2: [1, 3], 3: [1, 2], 4: []}
        self.poly_a._edge_graph = edge_graph
        poly_b._edge_graph = edge_graph
        self.assertEqual(self.poly_a, poly_b)

    def test_eq_false_for_different_edge_graphs(self):
        poly_b = copy.deepcopy(self.poly_a)
        self.poly_a._edge_graph = {1: [2, 3], 2: [1, 3], 3: [1, 2], 4: []}
        poly_b._edge_graph = {1: [2], 2: [1], 3: [4], 4: [3]}
        self.assertNotEqual(self.poly_a, poly_b)

    def test_eq_returns_not_implemented_for_non_polyhedron(self):
        result = self.poly_a.__eq__('not a polyhedron')
        self.assertIs(result, NotImplemented)


class TestCoordinationPolyhedronFromSitesOctahedron(unittest.TestCase):
    """Tests using a real octahedron built from pymatgen PeriodicSite objects."""

    def setUp(self):
        from pymatgen.core import Lattice
        from pymatgen.core.sites import PeriodicSite
        lattice = Lattice.cubic(10.0)
        central = PeriodicSite('Ti', [0.5, 0.5, 0.5], lattice)
        vertices = [
            PeriodicSite('O', [0.6, 0.5, 0.5], lattice),
            PeriodicSite('O', [0.4, 0.5, 0.5], lattice),
            PeriodicSite('O', [0.5, 0.6, 0.5], lattice),
            PeriodicSite('O', [0.5, 0.4, 0.5], lattice),
            PeriodicSite('O', [0.5, 0.5, 0.6], lattice),
            PeriodicSite('O', [0.5, 0.5, 0.4], lattice),
        ]
        self.poly = CoordinationPolyhedron.from_sites(central, vertices)

    def test_best_fit_geometry_is_octahedron(self):
        result = self.poly.best_fit_geometry
        self.assertEqual(result['geometry'], 'Octahedron')
        self.assertAlmostEqual(result['symmetry_measure'], 0.0)

    def test_best_fit_geometry_has_expected_keys(self):
        result = self.poly.best_fit_geometry
        self.assertIn('geometry', result)
        self.assertIn('symmetry_measure', result)

    def test_volume_of_regular_octahedron(self):
        # Vertices at ±1 from centre: volume = 4/3 * a^3 where a = 1
        self.assertAlmostEqual(self.poly.volume, 4.0 / 3.0, places=5)

    def test_edge_graph_each_vertex_has_four_neighbours(self):
        for neighbours in self.poly.edge_graph.values():
            self.assertEqual(len(neighbours), 4)

    def test_edge_graph_opposite_vertices_not_connected(self):
        edge_graph = self.poly.edge_graph
        vi = self.poly.vertex_indices
        # Vertices 0,1 are ±x; 2,3 are ±y; 4,5 are ±z
        opposite_pairs = [(vi[0], vi[1]), (vi[2], vi[3]), (vi[4], vi[5])]
        for v1, v2 in opposite_pairs:
            self.assertNotIn(v2, edge_graph[v1])
            self.assertNotIn(v1, edge_graph[v2])

    def test_faces_returns_eight_triangles(self):
        faces = self.poly.faces()
        self.assertEqual(len(faces), 8)
        for face in faces:
            self.assertEqual(len(face), 3)

    def test_faces_are_sorted_tuples(self):
        for face in self.poly.faces():
            self.assertIsInstance(face, tuple)
            self.assertEqual(list(face), sorted(face))


class TestMinimumImageVertexCoordinates(unittest.TestCase):

    def test_vertex_wraps_across_periodic_boundary(self):
        from pymatgen.core import Lattice
        from pymatgen.core.sites import PeriodicSite
        lattice = Lattice.cubic(10.0)
        central = PeriodicSite('Ti', [0.95, 0.5, 0.5], lattice)
        vertices = [
            PeriodicSite('O', [0.05, 0.5, 0.5], lattice),  # wraps across x
            PeriodicSite('O', [0.95, 0.6, 0.5], lattice),
            PeriodicSite('O', [0.95, 0.4, 0.5], lattice),
            PeriodicSite('O', [0.95, 0.5, 0.6], lattice),
            PeriodicSite('O', [0.95, 0.5, 0.4], lattice),
            PeriodicSite('O', [0.85, 0.5, 0.5], lattice),
        ]
        poly = CoordinationPolyhedron.from_sites(central, vertices)
        min_image_coords = poly.minimum_image_vertex_coordinates()
        # Vertex at frac [0.05] wraps to image at frac [1.05] → 10.5 Angstrom
        # Central atom is at 9.5 Angstrom, so distance is 1.0 Angstrom
        np.testing.assert_array_almost_equal(
            min_image_coords[0], [10.5, 5.0, 5.0])


class TestIntersectionNoSharedVertices(unittest.TestCase):

    def test_no_shared_vertices_returns_empty_tuple(self):
        mock_central_a = Mock(spec=Atom)
        mock_central_a.in_polyhedra = []
        mock_central_a.index = 0
        mock_central_a.label = 'Li'
        mock_central_b = Mock(spec=Atom)
        mock_central_b.in_polyhedra = []
        mock_central_b.index = 1
        mock_central_b.label = 'Li'
        verts_a = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(verts_a, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        verts_b = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(verts_b, 11):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly_a = CoordinationPolyhedron(
            central_atom=mock_central_a, vertices=verts_a)
        poly_b = CoordinationPolyhedron(
            central_atom=mock_central_b, vertices=verts_b)
        self.assertEqual(poly_a.intersection(poly_b), ())


class TestPolyhedronNotOwnNeighbour(unittest.TestCase):

    def test_polyhedron_is_not_its_own_neighbour(self):
        mock_central = Mock(spec=Atom)
        mock_central.in_polyhedra = []
        mock_central.index = 0
        mock_central.label = 'Li'
        mock_vertices = [Mock(spec=Atom) for _ in range(4)]
        for i, v in enumerate(mock_vertices, 1):
            v._neighbours = {}
            v.index = i
            v.__lt__ = mock_atom_lt
            v.in_polyhedra = []
        poly = CoordinationPolyhedron(
            central_atom=mock_central, vertices=mock_vertices)
        poly._edge_graph = Mock()
        self.assertNotIn(poly, poly.neighbours())


if __name__ == '__main__':
    unittest.main()
