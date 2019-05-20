from functools import partial
import numpy as np
from .coordination_polyhedron import CoordinationPolyhedron
from .utils import flatten

def matching_sites( structure, reference_structure, species=None ):
    """
    Returns a subset of sites from structure (as a list) where each site is the closest to one 
    site in the reference structure.
    
    Args:
        structure (Structure): The structure being analysed.
        reference_structure (Structure): A Structure object containing a set of reference sites.
        species (:obj:`list(str)`, optional): Optional list of strings specifying a subset of species
            to return. Default is None, which specifies all species matching species are returned.
        
    Returns:
        (list([Site,index]))
    """
    matched_sites = []
    for ref_site in reference_structure:
        dr = [ site.distance( ref_site ) for site in structure ]
        i = np.argmin( dr )
        matched_sites.append( [ structure[ i ], i ] )
    if species:
        matched_sites = [ s for s in matched_sites if str( s[0].specie ) in species ]
    return matched_sites

def matching_site_indices( structure, reference_structure, species=None ):
    """
    Returns a subset of site indices from structure (as a list) where each site is the closest to one 
    site in the reference structure.
    
    Args:
        structure (Structure): The structure being analysed.
        reference_structure (Structure): A Structure object containing a set of reference sites.
        species (:obj:`list[str]`, optional): A list of species labels. If this is set, only matching
            sites will be included in the returned set.
        
    Returns:
        (list[int])
    """
    return [ m[1] for m in matching_sites( structure, reference_structure, species=species ) ]

def create_matching_site_generator( reference_structure, species=None ):
    """
    Args:
        reference_structure (Structure):
        species (:obj:`list[str]`, optional):

    Returns:
        (functools.partial): 
    """
    return partial( matching_site_indices, reference_structure=reference_structure, species=species )

def generator_from_atom_argument( arg ):
    """
    Returns a generator function for selecting a subset of sites from a pymatgen :obj:`Structure` object.

    Args:
        arg (various): Argument used to construct the generator function.

    Returns:
        (func): Generator function that takes one argument (:obj:`Structure`) and
            returns an appropriate list of site indices.

    """
    if type( arg ) is str:
        return partial( _get_indices_from_str, arg=arg )
    elif type( arg ) in [ list, tuple ]:
        if type( arg[0] ) is str:
            return partial( _get_indices_from_list_str, arg=arg )
        elif type( arg[0] ) is int:
            return partial( _get_indices_from_list_int, arg=arg )
    elif callable( arg ):
        return arg

def _get_indices_from_str( structure, arg ):
    return structure.indices_from_symbol( arg )

def _get_indices_from_list_str( structure, arg ):
    return flatten( [ s.indices_from_symbol( sp ) for sp in arg ] )

def _get_indices_from_list_int( structure, arg ):
    return arg

class PolyhedraRecipe:

    allowed_methods = [ 'distance cutoff', 'closest centre' ]

    def __init__( self, method, central_atoms, vertex_atoms, 
                  coordination_cutoff=None, vertex_graph_cutoff=None,
                  label=None ):
        """
        Create a :obj:`PolyhedraRecipe` object.

        Args:
            method (str): Method used for constructing coordination polyhedra. 
                Implemented options are:

                    - `distance cutoff`: include all coordinating ions within a cutoff distance.
                    - `closest centre`: include all coordinating ions that share a common closest centre atom.
            coordination_cutoff (:obj:`float`, optional): Cutoff distance for vertex atoms to
                be considered as coordinated to a central atom.
            central_atoms (various): Defines the set of atoms considered polyhedron centres.
            vertex_atoms  (various): Defines the set of atoms considered polyhedron vertices.
            vertex_graph_cutoff (:obj:`float`, optional): TODO.
            label (str): Label for this recipe.

        Returns:
            None

        Notes:
            `central_atoms` and `vertex_atoms` both define how to select a subset of
            atoms from a single :obj:`Structure`. How this is implemented depends on
            the argument type provided when initialising a :obj:`PolyhedraRecipe` instance.

                    - (:obj:`str`): Select all atoms with matching species strings, e.g. ``'Ti'``.
                    - (:obj:`list` of :obj:`str`): Select all atoms with species strings that match
                        one of the list entries, e.g. ``[ 'Ti', 'Nb' ]``.
                    - (:obj:`list(int)`): Select matching atom indices, e.g. ``[ 1, 2, 3 ]``.
                    - (:obj:`func`): Generator function. This can be any function that 
                        accepts a pymatgen `Structure` as its only argument, and returns
                        a list of selected atom indices, e.g.
                        ``lambda s: s.indices_from_symbol('Ti')``
                        gives the same result as passing ``'Ti'``. 

        """
        if method not in PolyhedraRecipe.allowed_methods:
            methods_string = '\n'.join( [ "    - '{}'".format( m ) 
                for m in PolyhedraRecipe.allowed_methods ] )
            raise ValueError( "\n'{}' is not a recognised string for selecting the recipe method.\n Valid methods are:\n{}".format( 
                              method, methods_string ) )
        self.method = method
        self._central_atom_list_generator = generator_from_atom_argument( central_atoms )
        self._vertex_atom_list_generator = generator_from_atom_argument( vertex_atoms )
        self._central_atom_list = None
        self._vertex_atom_list = None
        self.coordination_cutoff = coordination_cutoff
        self.vertex_graph_cutoff = vertex_graph_cutoff
        self.label = label

    def central_atom_list( self, structure=None, recalculate=False ):
        if structure:
            if recalculate or not self._central_atom_list:
                self._central_atom_list = self._central_atom_list_generator( structure )
        elif not self._central_atom_list or recalculate:
            raise ValueError( 'Needs structure argument' )
        return self._central_atom_list

    def vertex_atom_list( self, structure=None, recalculate=False ):
        if structure:
            if recalculate or not self._vertex_atom_list:
                self._vertex_atom_list = self._vertex_atom_list_generator( structure )
        elif not self._vertex_atom_list or recalculate:
            raise ValueError( 'Needs structure argument' )
        return self._vertex_atom_list

    def find_polyhedra( self, atoms, structure=None ):
        polyhedra_method = { 'distance cutoff': partial( polyhedra_from_distance_cutoff, cutoff=self.coordination_cutoff ),
                             'closest centre': polyhedra_from_closest_centre }
        central_atom_list = self.central_atom_list( structure )
        vertex_atom_list = self.vertex_atom_list( structure )
        central_atoms = [ atom for atom in atoms if atom.index in central_atom_list ]
        vertex_atoms = [ atom for atom in atoms if atom.index in vertex_atom_list ]
        return polyhedra_method[ self.method ]( central_atoms=central_atoms, 
                                                vertex_atoms=vertex_atoms, 
                                                label=self.label )

def polyhedra_from_distance_cutoff( central_atoms, vertex_atoms, cutoff, label=None ):
    polyhedra = []
    for c_atom in central_atoms:
        vertices = [ a for a in vertex_atoms if a.site.distance( c_atom.site ) <= cutoff ]
        polyhedra.append( CoordinationPolyhedron( central_atom=c_atom, 
                                                  vertices=vertices, 
                                                  label=label ) )
    return polyhedra

def polyhedra_from_closest_centre( central_atoms, vertex_atoms, label=None ):
    coordination_coords = [ co_atom.coords for co_atom in vertex_atoms ]
    central_coords = [ c_atom.coords for c_atom in central_atoms ]
    closest_site_index = [ np.argmin( [ a.site.distance( c_atom.site ) 
                           for c_atom in central_atoms ] ) 
                           for a in vertex_atoms ]
    polyhedra = []
    for i, c_atom in enumerate( central_atoms ):
        vertices = [ co for c, co in zip( closest_site_index, vertex_atoms ) if c == i ]
        polyhedra.append( CoordinationPolyhedron( central_atom=c_atom, 
                                                  vertices=vertices, 
                                                  label=label ) )
    return polyhedra

