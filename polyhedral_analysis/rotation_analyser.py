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
        if np.linalg.det(op.rotation_matrix) > 0.5
    ])
    if len(proper_rotations) == 0:
        raise RuntimeError(
            'No proper rotations found for the reference geometry.')

    # Convert symmetry operations to permutation representations for bsym
    coords = mol.cart_coords
    seen: set[tuple[int, ...]] = set()
    mappings: list[list[int]] = []
    for op in symm_ops:
        mapped_coords = op.operate_multi(coords)
        new_mol = Molecule(mol.species, mapped_coords)
        mapping_zero = coord_list_mapping(new_mol.cart_coords, coords)
        if len(set(mapping_zero)) != n:
            raise RuntimeError(
                'coord_list_mapping did not produce a bijection; '
                'symmetry operation maps multiple vertices to the same target.')
        mapping = [x + 1 for x in mapping_zero]
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

    All reference orientations must be rotated copies of the same geometry
    (i.e. share the same point group). The symmetry analysis is performed
    once from the first reference orientation and reused for all others.

    Attributes:
        reference_points (np.ndarray): An MxNx3 numpy array of points around
            the origin ([0.0, 0.0, 0.0]) that define sets of centre-vertex
            vectors for specific reference orientations, e.g. for a single
            orientation of a tetrahedron, reference points will be a size
            1x4x3 array, e.g.::

                [[[ 1.0, -1.0,  1.0],
                  [-1.0, -1.0, -1.0],
                  [ 1.0,  1.0, -1.0],
                  [-1.0,  1.0,  1.0]]]

    """

    def __init__(self,
                 reference_points: np.ndarray) -> None:
        """Initialise a RotationAnalyser instance.

        Args:
            reference_points: Either an Nx3 or MxNx3 numpy array of points
                around the origin ([0.0, 0.0, 0.0]) that define sets of
                centre-vertex vectors, and are used to classify the rotational
                orientation of polyhedra.

                If an Nx3 array is passed in, this will be converted to a
                1xNx3 array.
        """
        if len(reference_points.shape) == 2:
            reference_points = np.array([reference_points])
        if len(reference_points.shape) == 1:
            raise ValueError(
                'Reference points must all contain the same number of coordinates')
        self.reference_points = reference_points
        self._symmetry_data: tuple[np.ndarray, np.ndarray] | None = None

    def _ensure_symmetry_analysis(self) -> tuple[np.ndarray, np.ndarray]:
        """Run point group analysis if not already done.

        Returns:
            A tuple of (reduced_permutations, proper_rotations) from
            a single symmetry analysis.
        """
        if self._symmetry_data is None:
            self._symmetry_data = _analyse_point_group(self.reference_points[0])
        return self._symmetry_data

    @property
    def reduced_permutations(self) -> np.ndarray:
        """Symmetry-inequivalent permutation indices for the reference geometry.

        Computed lazily on first access.
        """
        return self._ensure_symmetry_analysis()[0]

    @property
    def proper_rotations(self) -> np.ndarray:
        """Proper rotation matrices (Kx3x3) of the reference geometry's point group.

        Computed lazily on first access.
        """
        return self._ensure_symmetry_analysis()[1]

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
           of the reference points (N! / ``|G|`` CSM evaluations instead of N!).
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
        then expand that orbit to all ``|G|`` equivalent permutations and compute
        CSM for each to obtain exact rotation matrices (phase 2). This costs
        N! / ``|G|`` + ``|G|`` CSM evaluations per reference orientation rather than
        N! / ``|G|``, but avoids the composition approximation.

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
        expected_n = self.reference_points.shape[1]
        if points.shape != (expected_n, 3):
            raise ValueError(
                f'Expected points with shape ({expected_n}, 3), '
                f'got {points.shape}.')
        points = points - np.mean(points, axis=0, dtype=float)
        n_refs = len(self.reference_points)
        n_proper = len(self.proper_rotations)
        # Phase 1: find best alignment per reference orientation
        best_csm_per_ref = np.full(n_refs, float('inf'))
        best_rotation_per_ref: list[np.ndarray] = [np.eye(3) for _ in range(n_refs)]
        for i, rp in enumerate(self.reference_points):
            for perm in self.reduced_permutations:
                result = continuous_symmetry_measure(points, rp[perm])
                if result.symmetry_measure < best_csm_per_ref[i]:
                    best_csm_per_ref[i] = result.symmetry_measure
                    best_rotation_per_ref[i] = result.rotation_matrix
        # Phase 2: compose best-fit rotation with proper symmetry operations
        all_rot_distances = np.empty(n_refs * n_proper)
        for i in range(n_refs):
            composed = self.proper_rotations @ best_rotation_per_ref[i]
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
            symmetry_measure=float(best_csm_per_ref[reference_geometry_index]),
            all_rotational_distances=all_rot_distances,
        )

    def polyhedron_orientation(self,
                               polyhedron: CoordinationPolyhedron) -> OrientationDict:
        """Find the discrete orientation of a coordination polyhedron.

        Convenience wrapper around :meth:`discrete_orientation` that
        extracts vertex vectors from the polyhedron relative to its
        central atom.

        Args:
            polyhedron: The coordination polyhedron to analyse.

        Returns:
            Dictionary describing the orientation.
            See :meth:`discrete_orientation` for details.
        """
        points = polyhedron.vertex_vectors(reference='central_atom')
        return self.discrete_orientation(points)
