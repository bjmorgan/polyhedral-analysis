import numpy as np

def cos_theta(a: np.ndarray, 
              b: np.ndarray) -> float:
    """Calculate the cosine of the angle between two vectors.

    Args:
        a (np.ndarray): Vector a.
        b (np.ndarray): Vector b.
    
    Returns:
        float

    """
    to_return = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return to_return

def projection_xyz(vec_in: np.ndarray) -> float:
    """TODO"""
    scores = []
    for vec in np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]]):
        scores.append(2.0 * cos_theta(vec_in, vec)**2 - 1.0)
    return max(scores)

def oct_rotational_order_parameter(points: np.ndarray) -> float:
    """Order parameter for rotational orientation of an octahedron with respect to <100> vectors.

    Args:
        points: A 6x3 array of centroid-centred vertex vectors.

    Returns:
        1.0 for a perfectly aligned octahedron.
        0.33 for a 45 degree rotated octahedron around one axis.
    """
    return sum(projection_xyz(point) for point in points) / 6.0
