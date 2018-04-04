from .atom import Atom
from .utils import flatten

class Configuration:

    def __init__( self, structure, recipes, config_number=None ):
        self.structure = structure
        self.config_number = config_number
        self.atoms = [ Atom( index=i, site=site, label=site.species_string ) 
                       for i, site in enumerate( self.structure ) ]
        self.polyhedra = []
        for recipe in recipes:
            self.polyhedra.extend( recipe.find_polyhedra( self.atoms, self.structure ) )
        self.central_atoms = sorted( list( set( 
                                 [ p.central_atom for p in self.polyhedra ] ) ) )
        self.coordination_atoms = sorted( list( set( flatten( 
                                      [ p.vertices for p in self.polyhedra ] ) ) ) )

    def coordination_atoms_by_index( self, index ):
        return next( atom for atom in self.coordination_atoms if atom.index == index )

    @property
    def polyhedra_labels( self ):
        return [ p.label for p in self.polyhedra ]

    def polyhedra_by_label( self, label ):
        return [ p for p in self.polyhedra if p.label == label ]
