"""Continuous symmetry measure (CSM) via SVD-based Procrustes analysis.

References:
    Pinsky & Avnir, Inorganic Chemistry 37, 5575 (1998).
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class SymmetryMeasureResult(NamedTuple):
    """Result of a continuous symmetry measure calculation."""

    #: The CSM value as a percentage (0--100).
    symmetry_measure: float
    #: The optimal 3x3 rotation matrix.
    rotation_matrix: np.ndarray
    #: The optimal scaling factor.
    scaling_factor: float


def continuous_symmetry_measure(
    points_distorted: np.ndarray,
    points_perfect: np.ndarray,
) -> SymmetryMeasureResult:
    """Compute the continuous symmetry measure (CSM) of a set of points.

    Uses SVD-based orthogonal Procrustes analysis to find the optimal
    rotation and scaling that align the distorted points to the perfect
    reference, then computes the normalised residual as a percentage.

    Args:
        points_distorted: An Nx3 array of distorted vertex coordinates.
        points_perfect: An Nx3 array of ideal vertex coordinates.

    Returns:
        A ``SymmetryMeasureResult`` with fields ``symmetry_measure``
        (float, percentage 0--100), ``rotation_matrix`` (3x3 ndarray),
        and ``scaling_factor`` (float).
    """
    if len(points_distorted) == 1:
        return SymmetryMeasureResult(
            symmetry_measure=0.0,
            rotation_matrix=np.eye(3),
            scaling_factor=1.0,
        )
    if np.allclose(points_distorted, 0.0):
        raise ValueError(
            "All distorted points are at the origin; CSM is undefined."
        )
    # Orthogonal Procrustes via SVD
    H = points_distorted.T @ points_perfect
    U, _S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # Optimal scaling factor
    rotated = (R @ points_distorted.T).T
    scaling = float(
        np.tensordot(rotated, points_perfect, axes=2)
        / np.tensordot(rotated, rotated, axes=2)
    )
    # CSM as normalised residual (percentage)
    diff = points_perfect - scaling * rotated
    csm = float(
        np.tensordot(diff, diff, axes=2)
        / np.tensordot(points_perfect, points_perfect, axes=2)
        * 100.0
    )
    return SymmetryMeasureResult(
        symmetry_measure=csm,
        rotation_matrix=R,
        scaling_factor=scaling,
    )
