from .atom import Atom
from .coordination_polyhedron import CoordinationPolyhedron

class Configuration:

    def __init__( self, structure, recipes, config_number=None ):
        self.structure = structure
        self.config_number = config_number
        self.polyhedra = []
        self.central_atoms = []
        self.coordination_atoms = []
        for recipe in recipes:
            # TODO: these get repeated for each recipe!
            central_atoms = []
            coordination_atoms = []
            for i in recipe.central_atom_list( structure ):
                if i in self.central_atom_indices:
                    central_atoms.append( self.central_atoms_by_index_list( [ i ] )[0] )
                else: 
                    new_atom = Atom( i, structure[i], label=recipe.label )
                    central_atoms.append( new_atom )
                    self.central_atoms.append( new_atom )
            for i in recipe.coordination_atom_list( structure ):
                if i in self.coordination_atom_indices:
                    coordination_atoms.append( self.coordination_atoms_by_index_list( [ i ] )[0] )
                else:
                    new_atom = Atom( i, structure[i] )
                    coordination_atoms.append( new_atom )
                    self.coordination_atoms.append( new_atom )
            self.central_atoms += central_atoms
            self.coordination_atoms += coordination_atoms
            for c_atom in central_atoms:
                c_atom.in_polyhedra.append( c_atom.index )
                vertices = get_neighbouring_coordination_atoms( c_atom, coordination_atoms, recipe.coordination_cutoff )
                self.polyhedra.append( CoordinationPolyhedron( central_atom=c_atom,
                                                               vertices=vertices ) )
                                                           
    @property
    def central_atom_indices( self ):
        return [ a.index for a in self.central_atoms ]

    @property
    def coordination_atom_indices( self ):
        return [ a.index for a in self.coordination_atoms ]

    def central_atoms_by_index_list( self, index_list ):
        return [ a for a in self.index_list if a.index in index_list ]

    def coordination_atoms_by_index_list( self, index_list ):
        return [ a for a in self.coordination_atoms if a.index in index_list ]

    @property
    def polyhedra_labels( self ):
        return [ p.label for p in self.polyhedra ]

def get_neighbouring_coordination_atoms( atom, coordination_atoms, cutoff ):
    return [ a for a in coordination_atoms if a.site.distance( atom.site ) <= cutoff ]

    




