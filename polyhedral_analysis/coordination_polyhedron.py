from .atom import Atom
from .symmetry_measure import symmetry_measures_from_coordination
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry, symmetry_measure
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.util.coord import pbc_shortest_vectors
import numpy as np
from scipy.spatial import ConvexHull
from itertools import permutations

class CoordinationPolyhedron:

    def __init__( self, central_atom, vertices, label=None ):
        """
        Initialise a CoordinationPolyhedron object.
 
        Args:
            central_atom (Atom): the central atom.
            vertices (list(Atom)): A list of atoms that define the coordination environment.
            label (:obj:Str, optional): An optional string used to label this coordination polyhedron.
                if the label is not defined, the label of the central atom will be used.
        Returns:
            None
        """
        self.central_atom = central_atom
        self.central_atom.in_polyhedra.append( self )
        self.vertices = vertices
        if label:
            self.label = label
        else:
            self.label = central_atom.label
        for vertex, neighbour_list in zip( self.vertices, self.construct_edge_graph().values() ):
            if not vertex.neighbours:
                vertex.neighbours = {}
            vertex.neighbours[ self.index ] = neighbour_list

    def __repr__( self ):
        if self.label:
            to_return = 'Coordination Polyhedron {}\n'.format( self.label ) 
        else:
            to_return = 'Coordination Polyhedron\n'
        to_return += '{}\n'.format( self.central_atom )
        to_return += '----------\n'
        for v in self.vertices:
            to_return += '{}\n'.format( v )
        return to_return

    def intersection( self, other_polyhedron ):
        return list( set( self.vertex_indices ) & set( other_polyhedron.vertex_indices ) )
       
    @property
    def vertex_indices( self ):
        return [ v.index for v in self.vertices ]
 
    @property
    def coordination_number( self ):
        return len( self.vertices )

    @property
    def index( self ):
        return self.central_atom.index

    @property
    def edge_graph( self ):
        to_return = {}
        for v in self.vertices:
            to_return[ v.index ] = v.neighbours[ self.index ]
        return to_return
   
    # should be calculated once upon initialisation and stored 
    @property
    def abstract_geometry( self ):
        """
        Returns the polyhedron as a Pymatgen AbstractGeometry object.
        """
        return AbstractGeometry( central_site=self.central_atom.coords, bare_coords=self.minimum_image_vertex_coordinates(), 
                                 include_central_site_in_centroid=False )

    @property
    def symmetry_measure( self ):
        if self.coordination_number not in symmetry_measures_from_coordination:
            raise ValueError( 'No symmetry measure objects for coordination number of {}'.format( self.coordination_number ) )
        msm = {}
        for string, sm in symmetry_measures_from_coordination[ self.coordination_number ].items():
            msm[ string ] = sm.minimum_symmetry_measure( self.abstract_geometry )
        return msm 

    @property
    def best_fit_geometry( self ):
        psm = self.symmetry_measure
        best_fit = min( psm, key=psm.get )
        return { 'geometry': best_fit, 'symmetry_measure': psm[ best_fit ] }

    def minimum_image_vertex_coordinates( self ):
        vertex_frac_coords = [ v.frac_coords for v in self.vertices ]
        pbc_vectors = pbc_shortest_vectors( self.central_atom.lattice, self.central_atom.frac_coords, vertex_frac_coords )[0]
        vertex_minimum_image_coords = [ self.central_atom.coords + v for v in pbc_vectors ]
        return vertex_minimum_image_coords

    def faces( self ):
        """
        Args:
            None

        Returns:
            (dict(int:list[int])
        """
        return [ [ self.vertex_indices[v] for v in simplex ] 
                 for simplex in merge_coplanar_simplices( self.convex_hull() ) ]

    def convex_hull( self ):
        return ConvexHull( self.minimum_image_vertex_coordinates() )

    def construct_edge_graph( self ):
        connected_vertices = { i : set() for i in range( self.coordination_number ) }
        if self.coordination_number > 3:
            convex_hull = self.convex_hull()
            for m in merge_coplanar_simplices( convex_hull ):
                if len(m) == 3:
                    for r in range(3):
                        rotated_simplex = np.roll(m,r)
                        connected_vertices[ rotated_simplex[0] ].add( rotated_simplex[1] )
                        connected_vertices[ rotated_simplex[0] ].add( rotated_simplex[2] )
                else: # non-triangular face with > 4 vertices. This will be a composite of more than one simplex.
                    component_simplices = []
                    for s in convex_hull.simplices:
                        if np.all( [ i in m for i in s ] ):
                            component_simplices.append( s )
                    # common elements are linked along an internal edge.
                    for s_roll in range( len( component_simplices ) ):
                        rotated_component_simplices = np.roll( component_simplices, s_roll, axis=0 )
                        this_simplex = rotated_component_simplices[0]
                        other_simplices = rotated_component_simplices[1:]
                        for roll in range(3):
                            internal_edge = True
                            rotated_simplex = np.roll( this_simplex, roll )
                            edge = [ rotated_simplex[0], rotated_simplex[1] ]
                            if not np.all( [ i in np.unique( other_simplices ) for i in edge ] ):
                                connected_vertices[ edge[0] ].add( edge[1] )
                            edge = [ rotated_simplex[0], rotated_simplex[2] ]
                            if not np.all( [ i in np.unique( other_simplices ) for i in edge ] ):
                                connected_vertices[ edge[0] ].add( edge[1] )
        else:
            for roll in range( self.coordination_number ):
                rotated_list = np.roll( list( range( self.coordination_number ) ), roll )
                for i in rotated_list[1:]:
                    connected_vertices[ rotated_list[0] ].add( i )
        edge_list = {}
        for i in range( self.coordination_number ):
            edge_list[ self.vertex_indices[i] ] = [ self.vertex_indices[v] for v in connected_vertices[i] ] 
        return edge_list

    @property
    def vertex_distances( self ):
        """
        Returns a list of distances from the central atom to the vertex atoms.

        Args:
            None

        Returns:
            list(float): a list of atomic separations.
        """ 
        return [ self.central_atom.site.distance( v.site ) for v in self.vertices ]

    def equal_vertices( self, other ):
        """
        Test whether this CoordinationPolyhedron has vertices with the same labels as
        another CoordinationPolyhedron.

        Args:
            other(CoordinationPolyhedron): The other CoordinationPolyhedron.

        Returns:
            (bool): True / False.
        """
        return self.vertex_indices == other.vertex_indices

    def equal_edge_graph( self, other ):
        """
        Test whether this CoordinationPolyhedron has the same edge graph as
        another CoordinationPolyhedron.

        Args:
            other(CoordinationPolyhedraon): The other CoordinationPolyhedrom.

        Returns:
            (bool): True or False.
        """
        return self.edge_graph == other.edge_graph

    def __eq__( self, other ):
        """
        Two CoordinationPolyhedron objects are considered equal if they
        have equal edge graphs.

        Args:
            other(CoordinationPolyhedraon): The other CoordinationPolyhedrom.

        Returns:
            (bool): True or False.
        """
        return self.equal_edge_graph( other )
      
def merge_coplanar_simplices( complex_hull, tolerance=0.1 ):
    triangles_to_merge = []
    # TODO: there has to be a better way of doing this pairwise loop, e.g. using itertools.permutations
    for i, e1 in enumerate( complex_hull.equations ):
        for j, e2 in enumerate( complex_hull.equations[i+1:], i+1 ):
            if np.all( e1 == e2 ):
                continue
            dr = e1 - e2
            distance = np.dot( dr, dr )
            if distance < tolerance:
                triangles_to_merge.append( [ i, j ] )
    merged_simplices = []
    for i, s in enumerate( complex_hull.simplices ):
        if i not in np.unique( triangles_to_merge ):
            merged_simplices.append( s )
    for i, j in triangles_to_merge:
        merged_simplices.append( np.unique( [ complex_hull.simplices[i], complex_hull.simplices[j] ] ) ) 
    return merged_simplices # note: this simplex is not ordered: should not be used to construct the edge_graph.
