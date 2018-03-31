class Atom:
    """Atom class"""

    def __init__( self, index, site, label=None ):
        """
        Initialise an Atom object.

        Args:
            index (int): Numerical index identifying this atom.
            site  (pymatgen.Site): The pymatgen Site (or PeriodicSite) object describing this atom.
            label (:obj:`str`, optional): An optional string labelling this atom.

        Attributes:
            in_polyhedra (list): List of polyhedra that this atom is part of.
            neighbours: (None|dict(list)): Dictionary of lists of neighbours for each 
                polyhedron this atom is part of. Dictionary keys are the indexes of 
                each respective polyhedron.
        """
        self.index = index
        self.site = site
        self.label = label
        self.in_polyhedra = []
        self.neighbours = None

    def __repr__( self ):
        return "{} {}".format( self.index, self.site )

#    def construct_neighbour_list( self, cutoff, vertex_atoms, poly_index ):
#        vertex_indices = [ v.index for v in vertex_atoms ]
#        if not self.neighbours:
#            self.neighbours = {}
#        self.neighbours[ poly_index ] = []
#        self.neighbours[ poly_index ] = [ v.index for v in vertex_atoms if self.site.distance( v.site ) < cutoff and v is not self ]

    @property
    def frac_coords( self ):
        return self.site.frac_coords

    @property
    def coords( self ):
        return self.site.coords

    @property
    def lattice( self ):
        return self.site.lattice

    def __lt__( self, other ):
        return self.index < other.index

    def __eq__( self, other ):
        return self.index == other.index

    def __hash__( self ):
        return self.index
