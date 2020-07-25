from typing import List, TypeVar, Sequence, Dict
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
                      neighbour_list: Dict[int, List[int]]):
    """TODO"""
    string = f'site: {polyhedron.index}\n'
    string += f'centre: {" ".join(f"{c:.8f}" for c in polyhedron.central_atom.coords)}\n'
    string += f'neighbours: {" ".join(str(i) for i in neighbour_list[polyhedron.index])}\n'
    string += f'label: {polyhedron.label}\n'
    return string
