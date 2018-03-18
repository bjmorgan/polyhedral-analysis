from functools import partial
import numpy as np

def matching_site_indices( structure, reference_structure, species=None ):
    """
    Returns a subset of sites from structure (as a list) where each site is the closest to one 
    site in the reference structure.
    
    Args:
        structure (Structure): The structure being analysed.
        reference_structure (Structure): A Structure object containing a set of reference sites.
        species (:obj:`list[str]`, optional): A list of species labels. If this is set, only matching
            sites will be included in the returned set.
        
    Returns:
        (list[Site,index])
    """
    matched_sites = []
    for ref_site in reference_structure:
        dr = [ site.distance( ref_site ) for site in structure ]
        i = np.argmin( dr )
        matched_sites.append( [ structure[ i ], i ] )
    if species:
        matched_sites = [ s for s in matched_sites if str( s[0].specie ) in species ]
    return [ m[1] for m in matched_sites ]

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
    """

    def __init__( self, coordination_cutoff, vertex_graph_cutoff,
                  central_atom_list=None, coordination_atom_list=None,
                  central_atom_list_generator=None, coordination_atom_list_generator=None,
                  label=None ):
        if not (central_atom_list or central_atom_list_generator):
            raise ValueError
        if not (coordination_atom_list or coordination_atom_list_generator):
            raise ValueError
        self._central_atom_list = central_atom_list
        self._coordination_atom_list = coordination_atom_list
        self._central_atom_list_generator = central_atom_list_generator
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

    






