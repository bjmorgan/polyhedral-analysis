from functools import partial
import numpy as np
from .coordination_polyhedron import CoordinationPolyhedron

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

class PolyhedraRecipe:
    """
    Defines a "recipe" for defining coordination polyhedra. 

    The algorithm for defining the coordination polyhedra must be specified using the keyword `method`.

    Available options are:
        - `distance cutoff`: include all coordinating ions within a cutoff distance.
        - `closest centre`: include all coordinating ions that share a common closest centre atom.
    """

    allowed_methods = [ 'distance cutoff', 'closest centre' ]

    def __init__( self, method, coordination_cutoff=None, vertex_graph_cutoff=None,
                  central_atom_list=None, coordination_atom_list=None,
                  central_atom_list_generator=None, coordination_atom_list_generator=None,
                  label=None ):
        if method not in PolyhedraRecipe.allowed_methods:
            raise ValueError( '{} is not a valid recipe method:\n valid methods: {}'.format( 
                              method, PolyhedraRecipe.allowed_methods ) )
        self.method = method
        if not (central_atom_list or central_atom_list_generator):
            raise ValueError( 'central_atom_list or central_atom_list_generator must be specified' )
        self._central_atom_list = central_atom_list
        self._central_atom_list_generator = central_atom_list_generator
        if not (coordination_atom_list or coordination_atom_list_generator):
            raise ValueError( 'coordination_atom_list or coordination_atom_list_generator must be specified' )
        self._coordination_atom_list = coordination_atom_list
        self._coordination_atom_list_generator = coordination_atom_list_generator
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

    def coordination_atom_list( self, structure=None, recalculate=False ):
        if structure:
            if recalculate or not self._coordination_atom_list:
                self._coordination_atom_list = self._coordination_atom_list_generator( structure )
        elif not self._coordination_atom_list or recalculate:
            raise ValueError( 'Needs structure argument' )
        return self._coordination_atom_list

    def find_polyhedra( self, atoms, structure=None ):
        polyhedra_method = { 'distance cutoff': partial( polyhedra_from_distance_cutoff, cutoff=self.coordination_cutoff ),
                             'closest centre': polyhedra_from_closest_centre }
        central_atom_list = self.central_atom_list( structure )
        coordination_atom_list = self.coordination_atom_list( structure )
        central_atoms = [ atom for atom in atoms if atom.index in central_atom_list ]
        coordination_atoms = [ atom for atom in atoms if atom.index in coordination_atom_list ]
        return polyhedra_method[ self.method ]( central_atoms=central_atoms, 
                                                coordination_atoms=coordination_atoms, 
                                                label=self.label )

def polyhedra_from_distance_cutoff( central_atoms, coordination_atoms, cutoff, label=None ):
    polyhedra = []
    for c_atom in central_atoms:
        vertices = [ a for a in coordination_atoms if a.site.distance( c_atom.site ) <= cutoff ]
        polyhedra.append( CoordinationPolyhedron( central_atom=c_atom, 
                                                  vertices=vertices, 
                                                  label=label ) )
    return polyhedra

def polyhedra_from_closest_centre( central_atoms, coordination_atoms, label=None ):
    coordination_coords = [ co_atom.coords for co_atom in coordination_atoms ]
    central_coords = [ c_atom.coords for c_atom in central_atoms ]
    closest_site_index = [ np.argmin( [ a.site.distance( c_atom.site ) 
                           for c_atom in central_atoms ] ) 
                           for a in coordination_atoms ]
    polyhedra = []
    for i, c_atom in enumerate( central_atoms ):
        vertices = [ co for c, co in zip( closest_site_index, coordination_atoms ) if c == i ]
        polyhedra.append( CoordinationPolyhedron( central_atom=c_atom, 
                                                  vertices=vertices, 
                                                  label=label ) )
    return polyhedra

