class Atom:
    """
    Class definition for Atom instances.
    """

    def __init__( self, index, site, label=None ):
        """
        Initialise an Atom index.

        Args:
            index:
            site:
            label (:obj:(Str), optional):
        """
        self.index = index
        self.site = site
        self.label = label
        self.in_polyhedra = []
        self.neighbours = None

    def __repr__( self ):
        return "{} {}".format( self.index, self.site )

    def construct_neighbour_list( self, structure, cutoff, vertex_atoms, poly_index ):
        vertex_indices = [ v.index for v in vertex_atoms ]
        if not self.neighbours:
            self.neighbours = {}
        self.neighbours[ poly_index ] = []
        self.neighbours[ poly_index ] = [ v.index for v in vertex_atoms if self.site.distance( v.site ) < cutoff and v is not self ]

    @property
    def frac_coords( self ):
        return self.site.frac_coords

    @property
    def coords( self ):
        return self.site.coords

    @property
    def lattice( self ):
        return self.site.lattice
