import unittest
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from unittest.mock import Mock, MagicMock, patch, PropertyMock
import copy
import numpy as np
from pymatgen.core.sites import Site

def mock_atom_lt( self, other ):
    return self.index < other.index

def mock_atom_eq( self, other ):
    return self.index == other.index

class TestCoordinationPolyhedronInit( unittest.TestCase ):

    def test_coordination_polyhedron_is_initialised( self ):
        mock_central_atom = Mock( spec=Atom )
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 10
        mock_central_atom.label = 'Li'
        mock_vertices = [ Mock( spec=Atom ) for i in range(6) ]
        for i, v in enumerate( mock_vertices, 1 ):
            v.neighbours = None
            v.__lt__ = mock_atom_lt
            v.index = i
            v.in_polyhedra = []
        with patch( 'polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.construct_edge_graph' ) as mock_construct_edge_graph:
            mock_construct_edge_graph.return_value = { 0: [ 1, 2, 3, 4 ],
                                                       1: [ 0, 2, 3, 5 ],
                                                       2: [ 0, 1, 3, 5 ],
                                                       3: [ 0, 2, 4, 5 ],
                                                       4: [ 0, 1, 3, 5 ],
                                                       5: [ 1, 2, 3, 4 ] }
            with patch( 'polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.construct_abstract_geometry' ) as mock_construct_abstract_geometry:
                mock_construct_abstract_geometry.return_value = Mock( spec=AbstractGeometry )
                CoordinationPolyhedron( central_atom=mock_central_atom, 
                                        vertices=mock_vertices )

class TestCoordinationPolyhedron(unittest.TestCase):

    def setUp( self ):
        mock_central_atom = Mock(spec=Atom)
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 0
        mock_central_atom.label = 'A'
        mock_central_atom.__eq__ = mock_atom_eq
        mock_vertices = [Mock(spec=Atom) for i in range(6)]
        for i, v in enumerate(mock_vertices, 1):
            v.neighbours = None
            v.index = i
            v.__lt__ = mock_atom_lt
            v.__eq__ = mock_atom_eq
            v.in_polyhedra =[] 
        self.coordination_polyhedron = CoordinationPolyhedron(
                                           central_atom=mock_central_atom,
                                           vertices=mock_vertices)

    def test_equal_members(self):
        other_coordination_polyhedron = copy.deepcopy(self.coordination_polyhedron)
        other_coordination_polyhedron.vertices[0].neighbours = {0: [1, 2, 3]}
        other_coordination_polyhedron.vertices[4].neighbours = {4: [1, 3, 5]}
        self.assertTrue(self.coordination_polyhedron.equal_members(
                         other_coordination_polyhedron))

    def test_vertex_vectors(self):
        vectors = [np.array([1.0, 0.0, 0.0]),
                   np.array([0.0, 1.0, 2.0])]
        self.coordination_polyhedron._abstract_geometry = Mock(spec=AbstractGeometry)
        self.coordination_polyhedron.abstract_geometry.points_wocs_ctwocc.return_value = vectors
        returned_vectors = self.coordination_polyhedron.vertex_vectors
        for v1, v2 in zip(vectors, returned_vectors):
            np.testing.assert_equal(v1, v2)

    def test_angles(self):
        vertex_vectors = np.array( [ [ 1.0, 0.0, 0.0 ],
                                     [ 0.0, 1.0, 0.0 ],
                                     [ 0.0, -1.0, 0.0 ] ] )
        with patch('polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors', new_callable=PropertyMock) as mock_vertex_vectors:
            mock_vertex_vectors.return_value = vertex_vectors
            angles = self.coordination_polyhedron.angles()
            np.testing.assert_equal( angles, [ 90.0, 90.0, 180.0 ] )

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
             
if __name__ == '__main__':
    unittest.main()

