from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom
from typing import List, Tuple

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
    if not polyhedron.best_fit_geometry['geometry'] is 'Octahedron':
        raise ValueError( 'This polyhedron is not recognised a an octahedron' )

def opposite_vertex_pairs(polyhedron: CoordinationPolyhedron,
                          check: bool = True) -> VertexPairs:
    """
    For an octahedral polyhedron, find the pairs of vertices opposite each other.
   
    Args:
        polyhedron (:obj:`CoordinationPolyhedron`): The polyhedron to be analysed.
        check: (Optional, :obj:`bool`): Optional flag to set whether to check that this
            polyhedron is an octahedron. Defulat is `True`.

    Returns:
        (tuple): 3 pairs of vertex atoms.

    """
    if check:
        check_octahedra(polyhedron)
    vertex_pairs = []
    seen_indices: List[int] = []
    for v1 in polyhedron.vertices:
        v1_neighbours = [polyhedron.vertices[v] for v in v1.neighbours[polyhedron.index]]
        if v1.index in seen_indices:
            continue
        v2 = next(v for v in polyhedron.vertices if v not in [v1, *v1_neighbours])
        vertex_pairs.append((v1, v2))
        seen_indices.extend([v1.index, v2.index])
    return (vertex_pairs[0], vertex_pairs[1], vertex_pairs[2])

def opposite_vertex_distances(polyhedron: CoordinationPolyhedron,
                              check: bool = True) -> Tuple[float, float, float]:
    if check:
        check_octahedra(polyhedron)
    distances = tuple(vp[0].distance(vp[1]) for vp in opposite_vertex_pairs(polyhedron))
    return (distances[0], distances[1], distances[2])
