from fnmatch import fnmatch
from monty.io import zopen
import json
import os

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

    def as_dict( self ):
        """
        json-serializable :obj:`dict` representation of Atom.
        """
        d = { 'index': self.index,
              'site': self.site.as_dict(),
              'label': self.label,
              'in_polyhedra': [ p.index for p in self.in_polyhedra ],
              'neighbours': self.neighbours,
              '@module': self.__class__.__module__,
              '@class': self.__class__.__name__ }
        return d

    def __repr__( self ):
        return "{} {}".format( self.index, self.site )

    def to( self, fmt=None, filename=None ):
        """
        Outputs the structure to a file or string.
    
        Args:
            fmt (:obj:`str`, optional): Format to ouput to. Defaults to JSON unlees the filename
                is provided. If fmt is specified this overrides the filename extension
                in the filename.
            filename (:obj:`str`, optional): If provided, the output will be written to a file.
                If `fmt` is not specified, the format will be determined from the filename
                extension.
    
        Returns:
            (str) if filename is `None`. `None` otherwise.
        """
        if not filename:
            filename = ""
        if not fmt:
            fmt = "json"
        else:
            fmt = fmt.lower()
        fname = os.path.basename( filename )
        if fmt == "json" or fnmatch( fname.lower(), "*.json" ):
            s = json.dumps( self.as_dict() )
            if filename:
                with zopen( filename, 'wt' ) as f:
                    f.write( '{}'.format(s) )
            else:
                return s
        
    @property
    def frac_coords( self ):
        """Fractional coordinates"""
        return self.site.frac_coords

    @property
    def coords( self ):
        """Cartesian coordinates"""
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
