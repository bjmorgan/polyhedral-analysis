# Functions for analysing octahedra

def opposite_vertex_pairs( polyhedron ):
    """
    For an octahedral polyhedron, find the pairs of vertices opposite each other.
   
    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.

    Returns:
        (list): a list of 3 pairse of vertices.

    """
    assert polyhedron.best_fit_geometry['geometry'] is 'Octahedron'
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
