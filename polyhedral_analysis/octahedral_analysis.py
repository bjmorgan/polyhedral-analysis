# Functions for analysing octahedra

def check_octahedra( polyhedron ):
    """
    Check whether a polyhedron is considered an octahedron.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The :obj:`CoordinationPolyhedron`
            to be tested.

    Returns:
        None

    Raises:
        ValueError: If the "best fit" geometry (lowest continuous symmetry
            measure) is not an octahedron, a ValueError is raised.

    """
    if not polyhedron.best_fit_geometry['geometry'] is 'Octahedron':
        raise ValueError( 'This polyhedron is not recognised a an octahedron' )

def opposite_vertex_pairs( polyhedron, check=True ):
    """
    For an octahedral polyhedron, find the pairs of vertices opposite each other.
   
    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.

    Returns:
        (list): a list of 3 pairse of vertices.

    """
    if check:
        check_octahedra( polyhedron )
    vertex_pairs = []
    seen_indices = []
    for v1 in polyhedron.vertices:
        if v1.index in seen_indices:
            continue
        v2 = [ v for v in polyhedron.vertices 
            if v not in [ v1 ] + polyhedron.vertices_by_indices( v1.neighbours[ polyhedron.index ] ) ][0]
        vertex_pairs.append( [ v1, v2 ] )
        seen_indices.extend( [ v1.index, v2.index ] )
    return vertex_pairs

def opposite_vertex_distances( polyhedron, check=True ):
    if check:
        check_octahedra( polyhedron )
    distances = [ vp[0].site.distance( vp[1].site ) for vp in opposite_vertex_pairs( polyhedron ) ]
    return distances
