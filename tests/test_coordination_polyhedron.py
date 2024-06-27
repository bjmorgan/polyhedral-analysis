import unittest
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
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
        mock_vertex_distances = [2.0, 1.0, 1.0, 1.0, 1.0, 1.5]
        self.coordination_polyhedron.central_atom.distance = \
            Mock(side_effect=mock_vertex_distances)
        vertex_distances = self.coordination_polyhedron.vertex_distances()
        np.testing.assert_equal(vertex_distances, mock_vertex_distances)

    def test_vertex_distances_and_labels(self):
        mock_vertex_distances = (2.0, 1.0, 1.0, 1.0, 1.0, 1.5)
        mock_labels = ['O', 'O', 'F', 'F', 'F', 'F']
        self.coordination_polyhedron.central_atom = Mock(spec=Atom)
        self.coordination_polyhedron.vertex_distances = \
            Mock(return_value=mock_vertex_distances)
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_labels', new_callable=PropertyMock) as mock_vertex_labels:
            mock_vertex_labels.__get__ = Mock(return_value=mock_labels)
            output = self.coordination_polyhedron.vertex_distances_and_labels()
        np.testing.assert_equal(output,
                                list(zip(mock_vertex_distances, mock_labels)))

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

    def test_abstract_geometry_if_cached(self):
        abstract_geometry = Mock(spec=AbstractGeometry)
        self.coordination_polyhedron._abstract_geometry = abstract_geometry
        self.assertEqual(
            self.coordination_polyhedron.abstract_geometry, abstract_geometry)

    def test_abstract_geometry_if_not_cached(self):
        abstract_geometry = Mock(spec=AbstractGeometry)
        self.coordination_polyhedron.construct_abstract_geometry = Mock(
            return_value=abstract_geometry)
        self.assertEqual(
            self.coordination_polyhedron.abstract_geometry, abstract_geometry)
        self.assertEqual(
            self.coordination_polyhedron.construct_abstract_geometry.call_count, 1)

    def test_construct_abstract_geometry(self):
        polyhedron = self.coordination_polyhedron
        polyhedron.central_atom.coords = np.array([1.0, 2.0, 3.0])
        polyhedron.minimum_image_vertex_coordinates = Mock(return_value='foo')
        abstract_geometry = Mock(spec=AbstractGeometry)
        with patch(f'polyhedral_analysis.coordination_polyhedron.AbstractGeometry') as mock_abstract_geometry:
            mock_abstract_geometry.return_value = abstract_geometry
            ag = polyhedron.construct_abstract_geometry()
            self.assertEqual(ag, abstract_geometry)
            mock_abstract_geometry.assert_called_with(central_site=polyhedron.central_atom.coords,
                                                      bare_coords='foo',
                                                      include_central_site_in_centroid=False)

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
        mock_vectors = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
            [-1.0, 0.0, 0.0]
        ])
        self.coordination_polyhedron.vertex_vectors = Mock(return_value=mock_vectors)

        # Calculate expected distortions
        distances = np.array([1.0, 2.0, 3.0, 1.0])
        d_mean = np.mean(distances)
        relative_deviations = (distances - d_mean) / d_mean
        
        expected_mad_normalized = np.mean(np.abs(relative_deviations))
        expected_msd_normalized = np.mean(relative_deviations**2)
        expected_mad_non_normalized = np.mean(np.abs(distances - d_mean))
        expected_msd_non_normalized = np.mean((distances - d_mean)**2)

        # Test MSD (default)
        self.assertAlmostEqual(
            self.coordination_polyhedron.radial_distortion_parameter(),
            expected_msd_normalized
        )

        # Test MAD
        self.assertAlmostEqual(
            self.coordination_polyhedron.radial_distortion_parameter(method='MAD'),
            expected_mad_normalized
        )

        # Test non-normalized MSD
        self.assertAlmostEqual(
            self.coordination_polyhedron.radial_distortion_parameter(normalize=False),
            expected_msd_non_normalized
        )

        # Test non-normalized MAD
        self.assertAlmostEqual(
            self.coordination_polyhedron.radial_distortion_parameter(normalize=False, method='MAD'),
            expected_mad_non_normalized
        )


        # Verify that vertex_vectors was called
        self.coordination_polyhedron.vertex_vectors.assert_called_with(reference='centroid')

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

if __name__ == '__main__':
    unittest.main()
