from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from typing import List, Tuple
import numpy as np # type: ignore

# Functions for analysing octahedra

VertexPairs = Tuple[Tuple[Atom, Atom],
                    Tuple[Atom, Atom],
                    Tuple[Atom, Atom]]

def check_octahedra(polyhedron: CoordinationPolyhedron) -> None:
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
    if polyhedron.best_fit_geometry['geometry'] != 'Octahedron':
        raise ValueError( 'This polyhedron is not recognised a an octahedron' )

def opposite_vertex_pairs(polyhedron: CoordinationPolyhedron,
                          check: bool = True) -> VertexPairs:
    """For an octahedral polyhedron, find the pairs of vertices opposite each other.
   
    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this
            polyhedron is an octahedron. Default is `True`.

    Returns:
        (tuple): 3 pairs of vertex atoms.

    """
    if check:
        check_octahedra(polyhedron)
    assert {len(n) for n in polyhedron.edge_graph.values()} == {4}, "Edge graph does not describe an octahedron."
    vertex_pairs = []
    seen_indices: List[int] = []
    for v1 in polyhedron.vertices:
        if v1.index in seen_indices:
            continue
        v1_neighbours = polyhedron.vertices_by_indices(v1.neighbours[polyhedron.index])
        v2 = next(v for v in polyhedron.vertices if v not in [v1, *v1_neighbours])
        vertex_pairs.append((v1, v2))
        seen_indices.extend([v1.index, v2.index])
    return (vertex_pairs[0], vertex_pairs[1], vertex_pairs[2])

def opposite_vertex_distances(polyhedron: CoordinationPolyhedron,
                              check: bool = True) -> Tuple[float, float, float]:
    """For an octahedral polyhedron, return the distances between pairs of cis- vertices.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Returns:
        (tuple): a length-3 tuple of distances.

    """
    distances = tuple(vp[0].distance(vp[1]) 
                      for vp in opposite_vertex_pairs(polyhedron, check=check))
    return (distances[0], distances[1], distances[2])

def trans_vertex_vectors(polyhedron: CoordinationPolyhedron,
                         check: bool = True) -> List[np.ndarray]:
    """For an octahedral polyhedron, return the vectors between pairs of trans vertices.

    Args:
        polyhedron (:obj: `CoordinationPolyhedron`): The polyhedron to be analysed.
        check (Optional, :obj: `bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Returns:
        List[np.ndarray]: A list of three numpy arrays, each representing a vector
                          between a pair of trans vertices.

    Raises:
        ValueError: If the polyhedron is not recognized as an octahedron.
    """
    if check:
        check_octahedra(polyhedron)
    vertex_pairs = opposite_vertex_pairs(polyhedron, check=False)
    vectors = []
    for v1, v2 in vertex_pairs:
        vector = polyhedron.vertex_vectors()[polyhedron.vertex_internal_index_from_global_index(v2.index)] - \
                 polyhedron.vertex_vectors()[polyhedron.vertex_internal_index_from_global_index(v1.index)]
        vectors.append(vector)
    return vectors

def isomer_is_trans(polyhedron: CoordinationPolyhedron,
             check: bool = True) -> bool:
    """For an octahedral polyhedron with 2+4 coordination, return True for trans coordination,
    and False for cis coordination.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Return:
        (bool): `True` for trans, `False` for cis.

    """
    coord = list(polyhedron.vertex_count.values())
    if sorted(coord) != [2, 4]:
        raise ValueError(f'cis/trans not defined for {polyhedron.vertex_count.values()} coordination.')
    to_return = bool(np.all([vp[0].label == vp[1].label
                     for vp in opposite_vertex_pairs(polyhedron, check=check)]))
    return to_return

def isomer_is_cis(polyhedron: CoordinationPolyhedron,
                  check: bool = True) -> bool:
    """For an octahedral polyhedron with 2+4 coordination, return True for cis coordination,
    and False for trans coordination.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Return:
        (bool): `True` for cis, `False` for trans.

    """
    return not isomer_is_trans(polyhedron, check=check)

def isomer_is_fac(polyhedron: CoordinationPolyhedron,
                  check: bool = True) -> bool:
    """For an octahedral polyhedron with 3+3 coordination, return True for fac coordination,
    and False for mer coordination.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Return:
        (bool): `True` for fac, `False` for mer.

    """
    coord = list(polyhedron.vertex_count.values())
    if sorted(coord) != [3, 3]:
        raise ValueError(f'fac/mer not defined for {polyhedron.vertex_count.values()} coordination.')
    to_return = bool(np.all([vp[0].label != vp[1].label
                     for vp in opposite_vertex_pairs(polyhedron, check=check)]))
    return to_return

def isomer_is_mer(polyhedron: CoordinationPolyhedron,
                  check: bool = True) -> bool:
    """For an octahedral polyhedron with 3+3 coordination, return True for mer coordination,
    and False for fac coordination.

    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is `True`.

    Return:
        (bool): `True` for fac, `False` for mer.

    """
    return not isomer_is_fac(polyhedron, check=check)

def trans_vector_orthogonality(polyhedron: CoordinationPolyhedron,
                               check: bool = True) -> Tuple[float, float, float]:
    """
    For each trans pair vector, calculate the smallest angle between the vector and 
    the plane defined by the other two vectors.

    Args:
        polyhedron (CoordinationPolyhedron): The polyhedron to be analysed.
        check (bool, optional): Optional flag to set whether to check that this polyhedron
            is an octahedron. Default is True.

    Returns:
        Tuple[float, float, float]: A tuple containing the three smallest angles (in degrees) 
        between each vector and the plane defined by the other two vectors.

    Raises:
        ValueError: If the polyhedron is not recognized as an octahedron.

    """
    if check:
        check_octahedra(polyhedron)
    
    trans_vectors = trans_vertex_vectors(polyhedron, check=False)
    
    angles = []
    for i in range(3):
        v1 = trans_vectors[i]
        v2 = trans_vectors[(i+1)%3]
        v3 = trans_vectors[(i+2)%3]
        
        # Calculate the normal vector to the plane containing v2 and v3
        normal = np.cross(v2, v3)
        normal = normal / np.linalg.norm(normal)  # Normalize the normal vector
        # Calculate the angle between v1 and the normal
        v1_normalized = v1 / np.linalg.norm(v1)
        dot_product = np.dot(normal, v1_normalized)
        angle = np.arccos(np.clip(abs(dot_product), -1.0, 1.0))
        angle_degrees = np.degrees(angle)
        angle_degrees = np.degrees(angle)
        
        angles.append(angle_degrees)
    
    return tuple(angles)
