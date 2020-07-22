from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from typing import Optional, List

class PolyhedronTrajectory:
    """
    Collects a sequence of polyhedra from different configurations into a trajectory.
    """

    def __init__(self, 
                 polyhedra: List[CoordinationPolyhedron],
                 config_numbers: Optional[List[int]] = None) -> None:
        """
        Initialise a PolyhedronTrajectory object.

        Args:
            polyhedra (list(CoordinationPolyhedron)): A list of CoordinationPolyhedron objects.
            config_numbers (:obj:list(int), optional): An optional list of configuration numbers 
                for each of the polyhedra.
       
        Returns:
            None
        """
        self.polyhedra = polyhedra
        if config_numbers:
            if len(config_numbers) != len(polyhedra):
                raise ValueError('polyhedra and config_numbers lists must be the same length.')
            self.config_numbers = config_numbers
        else:
            self.config_numbers = list(range(1, len(self.polyhedra) + 1))
