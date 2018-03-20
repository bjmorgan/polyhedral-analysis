import unittest
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from unittest.mock import Mock, patch

class TestCoordinationPolyhedronInit( unittest.TestCase ):

    def test_coordination_polyhedron_is_initialised( self ):
        mock_central_atom = Mock( spec=Atom )
        mock_central_atom.in_polyhedra = []
        mock_central_atom.index = 10
        mock_central_atom.label = 'Li'
        mock_vertices = [ Mock( spec=Atom ) for i in range(6) ]
        for v in mock_vertices:
            v.neighbours = None
        with patch( 'polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.construct_edge_graph' ) as mock_construct_edge_graph:
            mock_construct_edge_graph.return_value = { 0: [ 1, 2, 3, 4 ],
                                                       1: [ 0, 2, 3, 5 ],
                                                       2: [ 0, 1, 3, 5 ],
                                                       3: [ 0, 2, 4, 5 ],
                                                       4: [ 0, 1, 3, 5 ],
                                                       5: [ 1, 2, 3, 4 ] }
            CoordinationPolyhedron( central_atom=mock_central_atom, 
                                    vertices=mock_vertices )

if __name__ == '__main__':
    unittest.main()

