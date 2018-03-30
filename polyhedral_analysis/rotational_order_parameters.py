import numpy as np

def cos_theta( a, b ):
    return np.dot( a, b ) / ( np.linalg.norm(a) * np.linalg.norm(b) )

def projection( vec_in ):
    scores = []
    for vec in np.array( [ [ 1.0, 0.0, 0.0 ],
                           [ 0.0, 1.0, 0.0 ],
                           [ 0.0, 0.0, 1.0 ] ] ):
        scores.append( 2.0 * cos_theta( vec_in, vec )**2 - 1.0 )
    return max( scores )

def oct_rotational_order_parameter( ag ):
    """
    Order parameter for rotational oriention of an octahedron with respect to <100> vectors.
    Gives 1.0 for a perfectly aligned octahedron.
    Gives 0.33 for a 45° rotated octahedron around one axis.
    """
    return sum( projection( point ) for point in ag.points_wocs_ctwocc() ) / 6.0

def orientation( ag ):
    to_return = []
    for point in ag.points_wocs_ctwocc():
        for vec in np.array( [ [ 1.0, 0.0, 0.0 ],
                               [ 0.0, 1.0, 0.0 ],
                               [ 0.0, 0.0, 1.0 ] ] ):
            to_return.append( cos_theta( vec, point ) )
    return np.array( to_return ).reshape(6,3)
