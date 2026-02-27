from __future__ import annotations

from polyhedral_analysis.csm import continuous_symmetry_measure
import numpy as np
from typing import TypedDict
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron


def _analyse_point_group(reference_points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Analyse the point group of a reference geometry.

    Performs a single symmetry analysis and returns both the
    symmetry-reduced permutation indices and the proper rotation matrices.

    Args:
        reference_points: An Nx3 array of reference vertex coordinates.

    Returns:
        A tuple of (reduced_permutations, proper_rotations) where:

        - reduced_permutations is an MxN integer array of symmetry-inequivalent
          permutation representatives.
        - proper_rotations is a Kx3x3 array of proper (det > 0) rotation
          matrices of the point group.
    """
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    from pymatgen.util.coord import coord_list_mapping
    from bsym import ConfigurationSpace  # type: ignore[import]
    from bsym import PointGroup as BsymPointGroup  # type: ignore[import]
    from bsym.symmetry_operation import SymmetryOperation  # type: ignore[import]

    n = len(reference_points)
    mol = Molecule(['H'] * n, reference_points - reference_points.mean(axis=0))
    pga = PointGroupAnalyzer(mol)
    symm_ops = pga.get_symmetry_operations()

    # Extract proper 3x3 rotation matrices
    proper_rotations = np.array([
        op.rotation_matrix for op in symm_ops
        if np.linalg.det(op.rotation_matrix) > 0
    ])

    # Convert symmetry operations to permutation representations for bsym
    coords = mol.cart_coords
    seen: set[tuple[int, ...]] = set()
    mappings: list[list[int]] = []
    for op in symm_ops:
        mapped_coords = op.operate_multi(coords)
        new_mol = Molecule(mol.species, mapped_coords)
        mapping = [x + 1 for x in coord_list_mapping(
            new_mol.cart_coords, coords)]
        key = tuple(mapping)
        if key not in seen:
            seen.add(key)
            mappings.append(mapping)

    sym_ops = [SymmetryOperation.from_vector(m) for m in mappings]
    point_group = BsymPointGroup(symmetry_operations=sym_ops)
    config_space = ConfigurationSpace(
        objects=list(range(n)),
        symmetry_group=point_group,
    )
    unique_configs = config_space.unique_configurations(
        site_distribution={i: 1 for i in range(n)},
    )
    reduced_permutations = np.array(
        [c.tolist() for c in unique_configs],
        dtype=np.intp,
    )

    return reduced_permutations, proper_rotations


class OrientationDict(TypedDict):
    orientation_index: int
    reference_geometry_index: int
    rotational_distance: float
    symmetry_measure: float
    all_rotational_distances: np.ndarray

class RotationAnalyser:
    """Class for analysing rotational orientation of polyhedra.

    Attributes:
        reference_points (np.ndarray): A Mx3xN numpy array of points around the
            origin ([0.0, 0.0, 0.0]) that define sets of centre-vertex
            vectors for specific reference orientations, e.g. for a single
            orientation of a tetrahedron, reference points will be a size 1x3x4
            array, e.g.::

                [[[ 1.0, -1.0,  1.0],
                  [-1.0, -1.0, -1.0],
                  [ 1.0,  1.0, -1.0],
                  [-1.0,  1.0,  1.0]]]

    """

    def __init__(self,
                 reference_points: np.ndarray) -> None:
        """Initialise a RotationAnalyser instance.

        Args:
            reference_points: Either a 3xN or Mx3xN numpy array of points
                around the origin ([0.0, 0.0, 0.0]) that define sets of
                centre-vertex vectors, and are used to classify the rotational
                orientation of polyhedra.

                If a 3xN array is passed in, this will be converted to a
                1x3xN array.
        """
        if len(reference_points.shape) == 2:
            reference_points = np.array([reference_points])
        if len(reference_points.shape) == 1:
            raise ValueError(
                'Reference points must all contain the same number of coordinates')
        self.reference_points = reference_points
        self._reduced_permutations: np.ndarray | None = None
        self._proper_rotations: np.ndarray | None = None

    def _ensure_symmetry_analysis(self) -> None:
        """Run point group analysis if not already done.

        Populates both ``_reduced_permutations`` and ``_proper_rotations``
        from a single symmetry analysis.
        """
        if self._reduced_permutations is None:
            self._reduced_permutations, self._proper_rotations = (
                _analyse_point_group(self.reference_points[0]))

    @property
    def reduced_permutations(self) -> np.ndarray:
        """Symmetry-inequivalent permutation indices for the reference geometry.

        Computed lazily on first access.
        """
        self._ensure_symmetry_analysis()
        return self._reduced_permutations  # type: ignore[return-value]

    @property
    def proper_rotations(self) -> np.ndarray:
        """Proper rotation matrices (Kx3x3) of the reference geometry's point group.

        Computed lazily on first access.
        """
        self._ensure_symmetry_analysis()
        return self._proper_rotations  # type: ignore[return-value]

    def discrete_orientation(self,
                             points: np.ndarray) -> OrientationDict:
        """Find the discrete "closest orientation" for an input polyhedron of points.

        For example, a tetrahedron has 12 proper rotation symmetry operations.
        A distorted tetrahedron with vertices approximately aligned along the
        unordered vectors::

            [[ 1.0, -1.0,  1.0],
             [-1.0, -1.0, -1.0],
             [ 1.0,  1.0, -1.0],
             [-1.0,  1.0,  1.0]]

        can be in one of 12 orientations. This method compares the input
        points to the reference orientations and returns the orientation that
        minimises the rotational distance.

        The algorithm uses a two-phase approach to avoid the N! cost of
        evaluating all vertex permutations:

        1. Find the best vertex alignment using symmetry-reduced permutations
           of the reference points (N!/|G| CSM evaluations instead of N!).
        2. Compose the best-fit Procrustes rotation R_best with all proper
           rotation matrices S_i of the reference geometry's point group to
           give the rotation matrix for each orientation: R_i = S_i @ R_best.
        3. Calculate rotational distances from the composed matrices.
        4. Return the orientation with the minimum rotational distance.

        The composed rotation matrices in step 2 are exact for undistorted
        geometries. For distorted geometries they are approximate, since the
        true Procrustes rotation for each symmetry-equivalent permutation
        differs slightly from the composed matrix. In practice the results
        are numerically indistinguishable from the brute-force approach for
        moderate distortions.

        If exact rotation matrices are needed for highly distorted
        geometries, an alternative is a two-phase orbit-expansion approach:
        use reduced permutations to identify the minimum-CSM orbit (phase 1),
        then expand that orbit to all |G| equivalent permutations and compute
        CSM for each to obtain exact rotation matrices (phase 2). This costs
        N!/|G| + |G| CSM evaluations per reference orientation rather than
        N!/|G|, but avoids the composition approximation.

        Args:
            points: Nx3 numpy array of points describing the coordinates of
                the input polyhedron, centred around zero.

        Returns:
            Dictionary describing the orientation, with keys:

            - ``orientation_index`` (int): Index of this particular orientation.
            - ``reference_geometry_index`` (int): Index of the reference geometry
              the closest discrete orientation is equivalent to.
            - ``rotational_distance`` (float): Angle of rotation from the
              relevant reference orientation, in radians.
            - ``all_rotational_distances`` (np.ndarray): Array of angles of
              rotation from each reference orientation.
            - ``symmetry_measure`` (float): The continuous symmetry measure
              (CSM) for this polyhedron.
        """
        points = points - np.mean(points, axis=0, dtype=float)
        # Phase 1: find best alignment using reduced permutations
        best_csm = float('inf')
        best_rotation = np.eye(3)
        best_ref_index = 0
        for i, rp in enumerate(self.reference_points):
            for perm in self.reduced_permutations:
                result = continuous_symmetry_measure(points, rp[perm])
                if result.symmetry_measure < best_csm:
                    best_csm = result.symmetry_measure
                    best_rotation = result.rotation_matrix
                    best_ref_index = i
        # Phase 2: compose best-fit rotation with proper symmetry operations
        n_proper = len(self.proper_rotations)
        all_rot_distances = np.empty(len(self.reference_points) * n_proper)
        for i, rp in enumerate(self.reference_points):
            # Find best rotation for this specific reference orientation
            if i == best_ref_index:
                R_i = best_rotation
            else:
                R_i = min(
                    (continuous_symmetry_measure(points, rp[perm])
                     for perm in self.reduced_permutations),
                    key=lambda r: r.symmetry_measure,
                ).rotation_matrix
            composed = self.proper_rotations @ R_i
            traces = np.trace(composed, axis1=1, axis2=2)
            start = i * n_proper
            all_rot_distances[start:start + n_proper] = np.arccos(
                np.clip((traces - 1.0) / 2.0, -1.0, 1.0))
        index = int(np.argmin(all_rot_distances))
        reference_geometry_index = index // n_proper
        return OrientationDict(
            orientation_index=index,
            reference_geometry_index=reference_geometry_index,
            rotational_distance=all_rot_distances[index],
            symmetry_measure=best_csm,
            all_rotational_distances=all_rot_distances,
        )

    def polyhedron_orientation(self,
                               polyhedron: CoordinationPolyhedron) -> OrientationDict:
        points = polyhedron.vertex_vectors(reference='central_atom')
        return self.discrete_orientation(points)
