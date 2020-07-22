from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry # type: ignore
import numpy as np # type: ignore

def cos_theta(a: float, b: float) -> float:
    """TODO"""
    to_return = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    assert isinstance(to_return, float)
    return to_return

def projection_xyz(vec_in: np.ndarray) -> float:
    """TODO"""
    scores = []
    for vec in np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]]):
        scores.append(2.0 * cos_theta(vec_in, vec)**2 - 1.0)
    return max( scores )

def oct_rotational_order_parameter(ag: AbstractGeometry) -> float:
    """Order parameter for rotational oriention of an octahedron with respect to <100> vectors.

    Gives 1.0 for a perfectly aligned octahedron.
    Gives 0.33 for a 45° rotated octahedron around one axis.

    """
    return sum(projection_xyz(point) for point in ag.points_wocs_ctwocc()) / 6.0
