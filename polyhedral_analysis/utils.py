from collections.abc import Sequence
from typing import TypeVar
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron

T = TypeVar('T')

"""
Utility functions
"""


def flatten(this_list: Sequence[Sequence[T]]) -> list[T]:
    """Flattens a nested list.

    Args:
        (list): A list of lists.

    Returns:
        (list): The flattened list.

    """
    return [item for sublist in this_list for item in sublist]


def lattice_mc_string(polyhedron: CoordinationPolyhedron,
                      neighbour_list: dict[int, tuple[int, ...]]) -> str:
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


def prune_neighbour_list(neighbours: dict[int, tuple[int, ...]],
                         indices: list[int]) -> dict[int, tuple[int, ...]]:
    """Prune a neighbour list to include only the specified indices.

    Filters both the keys and the values of the neighbour dictionary,
    keeping only entries whose keys appear in ``indices`` and removing
    any neighbour references to indices not in the set.

    Args:
        neighbours: A neighbour list mapping each polyhedron index
            to a tuple of its neighbour indices.
        indices: The subset of indices to retain.

    Returns:
        The pruned neighbour list.
    """
    pruned_neighbours = {}
    for k, v in neighbours.items():
        if k in indices:
            pruned_neighbours[k] = tuple(set(v).intersection(set(indices)))
    return pruned_neighbours
