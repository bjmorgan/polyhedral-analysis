from __future__ import annotations

import numpy as np

from polyhedral_analysis.csm import continuous_symmetry_measure
from polyhedral_analysis.reference_geometries import get_reference_geometry


def _compute_reduced_permutations(reference_points: np.ndarray) -> np.ndarray:
    """Compute symmetry-inequivalent permutation representatives for reference points.

    Uses bsym to find the point group of the reference geometry and enumerate
    one representative permutation per orbit under the symmetry group action.

    Args:
        reference_points: An Nx3 array of reference coordinates.

    Returns:
        An MxN integer array where each row is a permutation of [0, ..., N-1]
        and M = N! / |G|, where |G| is the size of the symmetry group.

    Note:
        The point group detection relies on pymatgen's PointGroupAnalyzer
        (via bsym). A known bug (materialsproject/pymatgen#4596) causes
        under-detection of symmetry for the Square-face capped trigonal prism,
        resulting in more representatives than strictly necessary. This does
        not affect correctness.
    """
    from pymatgen.core import Molecule
    from bsym.interface.pymatgen import point_group_from_molecule  # type: ignore[import]
    from bsym import ConfigurationSpace  # type: ignore[import]

    n = len(reference_points)
    mol = Molecule(['H'] * n, reference_points)
    point_group = point_group_from_molecule(mol)
    config_space = ConfigurationSpace(
        objects=list(range(n)),
        symmetry_group=point_group,
    )
    site_distribution = {i: 1 for i in range(n)}
    unique_configs = config_space.unique_configurations(
        site_distribution=site_distribution,
    )
    return np.array(
        [c.tolist() for c in unique_configs],
        dtype=np.intp,
    )


class SymmetryMeasure:

    def __init__(self,
                 reference_points: np.ndarray,
                 string: str) -> None:
        self.string = string
        self.reference_points = reference_points
        self._permutations: np.ndarray | None = None

    @property
    def permutations(self) -> np.ndarray:
        """Symmetry-inequivalent index permutations for this geometry.

        Computed lazily on first access from the reference points using bsym.
        """
        if self._permutations is None:
            self._permutations = _compute_reduced_permutations(self.reference_points)
        return self._permutations

    def minimum_symmetry_measure(self,
                                 distorted_points: np.ndarray) -> float:
        return min(
            continuous_symmetry_measure(distorted_points[perm], self.reference_points).symmetry_measure
            for perm in self.permutations
        )

    @classmethod
    def from_name(cls, name: str) -> SymmetryMeasure:
        reference_points = get_reference_geometry(name)
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
                                      'Hexagonal bipyramid',
                                      'Bicapped octahedron (opposed cap faces)',
                                      'Bicapped octahedron (cap faces with one atom in common)',
                                      'Bicapped octahedron (cap faces with one edge in common)']}

symmetry_measures_from_coordination: dict[int, dict[str, SymmetryMeasure]] = {}
for coordination_number, strings in symmetry_measures_to_construct.items():
    symmetry_measures_from_coordination[coordination_number] = {}
    for string in strings:
        symmetry_measures_from_coordination[coordination_number][string] = SymmetryMeasure.from_name(string)
