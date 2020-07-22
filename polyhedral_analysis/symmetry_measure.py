from __future__ import annotations

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import symmetry_measure # type: ignore
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries # type: ignore
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from itertools import permutations
import numpy as np # type: ignore
from typing import Dict, List

class SymmetryMeasure:

    def __init__(self, 
                 reference_points: np.ndarray,
                 string: str) -> None:
        self.string = string
        self.reference_points = reference_points

    def minimum_symmetry_measure(self,
                                 ag: AbstractGeometry) -> float:
        distorted_points = ag.points_wocs_csc()
        to_return =  min(symmetry_measure(np.array(p), self.reference_points )['symmetry_measure'] 
                         for p in permutations(distorted_points))
        assert isinstance(to_return, float)
        return to_return

    @classmethod
    def from_name(cls, name: str) -> SymmetryMeasure:
        reference_points = AllCoordinationGeometries().get_geometry_from_name(name).points
        return cls(reference_points, name)

symmetry_measures_to_construct = {4: ['Tetrahedron'],
                                  5: ['Trigonal bipyramid', 
                                      'Square pyramid'],
                                  6: ['Octahedron',
                                      'Trigonal prism'],
                                  7: ['Pentagonal bipyramid',
                                      'Square-face capped trigonal prism',
                                      'Face-capped octahedron'],
                                  8: ['Cube',
                                      'Square antiprism',
                                      'Square-face bicapped trigonal prism',
                                      'Triangular-face bicapped trigonal prism',
                                      'Dodecahedron with triangular faces',
                                      'Dodecahedron with triangular faces',
                                      'Hexagonal bipyramid',
                                      'Bicapped octahedron (opposed cap faces)', 
                                      'Bicapped octahedron (cap faces with one atom in common)', 
                                      'Bicapped octahedron (cap faces with one edge in common)']}

symmetry_measures_from_coordination: Dict[int, Dict[str, SymmetryMeasure]] = {}
for coordination_number, strings in symmetry_measures_to_construct.items():
    symmetry_measures_from_coordination[coordination_number] = {}
    for string in strings:
        symmetry_measures_from_coordination[coordination_number][string] = SymmetryMeasure.from_name(string)
