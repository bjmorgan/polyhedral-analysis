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
        vectors = [np.array([1.0, 0.0, 0.0]),
                   np.array([0.0, 1.0, 2.0])]
        self.coordination_polyhedron._abstract_geometry = Mock(
            spec=AbstractGeometry)
        self.coordination_polyhedron.abstract_geometry.points_wocs_ctwocc.return_value = vectors
        returned_vectors = self.coordination_polyhedron.vertex_vectors
        for v1, v2 in zip(vectors, returned_vectors):
            np.testing.assert_equal(v1, v2)

    def test_angles(self):
        vertex_vectors = np.array([[1.0, 0.0, 0.0],
                                   [0.0, 1.0, 0.0],
                                   [0.0, -1.0, 0.0]])
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', new_callable=PropertyMock) as mock_vertex_vectors:
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

    def test_neighbours_and_shared_vertices(self):
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
        self.assertEqual(polyhedron_i.neighbours_by_index_and_shared_vertices(), 
                         {polyhedron_j.index: (5, 6)})
        self.assertEqual(polyhedron_j.neighbours_by_index_and_shared_vertices(), 
                         {polyhedron_i.index: (5, 6)})

    def test_shares_face_returns_True_when_faces_match(self):
        polyhedron = self.coordination_polyhedron
        other_polyhedron = copy.deepcopy(self.coordination_polyhedron)
        # consider two face-sharing tetrahedra
        polyhedron.faces = Mock(return_value=(
            (1, 2, 3), (2, 3, 4), (1, 3, 4), (1, 2, 4)))
        other_polyhedron.faces = Mock(return_value=(
            (1, 2, 3), (2, 3, 5), (1, 3, 5), (1, 3, 5)))
        self.assertEqual(polyhedron.shares_face(other_polyhedron), True)

    def test_shares_face_returs_False_when_no_faces_match(self):
        polyhedron = self.coordination_polyhedron
        other_polyhedron = copy.deepcopy(self.coordination_polyhedron)
        # consider two edge-sharing tetrahedra
        polyhedron.faces = Mock(return_value=(
            (1, 2, 3), (2, 3, 4), (1, 3, 4), (1, 2, 4)))
        other_polyhedron.faces = Mock(return_value=(
            (1, 2, 6), (2, 6, 5), (1, 6, 5), (1, 6, 5)))
        self.assertEqual(polyhedron.shares_face(other_polyhedron), False)

    def test_shares_face_raises_TypeError_for_incorrect_type(self):
        with self.assertRaises(TypeError):
            self.coordination_polyhedron.shares_face('foo')


if __name__ == '__main__':
    unittest.main()
