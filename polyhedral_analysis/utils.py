from typing import List

"""
Utility functions
"""

def flatten(this_list: List[List]) -> List:
    """Flattens a nested list.

    Args:
        (list): A list of lists.

    Returns:
        (list): The flattened list.

    """
    return [item for sublist in this_list for item in sublist]
