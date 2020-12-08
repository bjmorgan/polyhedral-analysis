from typing import List, TypeVar, Sequence, Dict, Tuple
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron

T = TypeVar('T')

"""
Utility functions
"""


def flatten(this_list: Sequence[Sequence[T]]) -> List[T]:
    """Flattens a nested list.

    Args:
        (list): A list of lists.

    Returns:
        (list): The flattened list.

    """
    return [item for sublist in this_list for item in sublist]


def lattice_mc_string(polyhedron: CoordinationPolyhedron,
                      neighbour_list: Dict[int, Tuple[int, ...]]) -> str:
    """Returns a string representation of a polyhedron as a `lattice_mc`_
    site-input formatted site.

    .. _lattice_mc: https://github.com/bjmorgan/lattice_mc

    Args:
        polyhedron (CoordinationPolyhedron): The coordination polyhedron.
        neighbour_list (dict): Neighbour list dictionary.

    Returns:
        str

    Example:

        >>> nlist = {1: (3, 5, 8), ...}
        >>> lattice_mc_string(polyhedron, neighbour_list=nlist
        site: 1
        centre: 2.94651000 1.70116834 1.20290767
        neighbours: 3 5 8
        label: oct

    """
    string = f'site: {polyhedron.index}\n'
    string += f'centre: {" ".join(f"{c:.8f}" for c in polyhedron.central_atom.coords)}\n'
    string += f'neighbours: {" ".join(str(i) for i in neighbour_list[polyhedron.index])}\n'
    string += f'label: {polyhedron.label}\n'
    return string


def prune_neighbour_list(neighbours: Dict[int, Tuple[int, ...]], 
                         indices: List[int]) -> Dict[int, Tuple[int, ...]]:
    """TODO"""
    pruned_neighbours = {}
    for k, v in neighbours.items():
        if k in indices:
            pruned_neighbours[k] = tuple(set(v).intersection(set(indices)))
    return pruned_neighbours
