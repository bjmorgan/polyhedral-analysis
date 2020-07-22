from typing import List, TypeVar, Sequence

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

