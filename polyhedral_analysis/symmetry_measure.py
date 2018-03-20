from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import symmetry_measure
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from itertools import permutations
import numpy as np

class SymmetryMeasure:

    def __init__( self, reference_points, string ):
        self.string = string
        self.reference_points = reference_points

    def minimum_symmetry_measure( self, ag ):
        distorted_points = ag.points_wocs_csc()
        return min( symmetry_measure( np.array(p), self.reference_points )['symmetry_measure'] 
                    for p in permutations( distorted_points ) )

    @classmethod
    def from_name( cls, name ):
        reference_points = AllCoordinationGeometries().get_geometry_from_name( name ).points
        return cls( reference_points, name )

symmetry_measures_to_construct = { 4: [ 'Tetrahedron' ],
                                   5: [ 'Trigonal bipyramid', 'Square pyramid' ],
                                   6: [ 'Octahedron', 'Trigonal prism' ],
                                   7: [ 'Pentagonal bipyramid', 'Square-face capped trigonal prism', 'Face-capped octahedron' ] }

symmetry_measures_from_coordination = {}
for coordination_number, strings in symmetry_measures_to_construct.items():
    symmetry_measures_from_coordination[ coordination_number ] = {}
    for string in strings:
        symmetry_measures_from_coordination[ coordination_number ][ string ] = SymmetryMeasure.from_name( string )
