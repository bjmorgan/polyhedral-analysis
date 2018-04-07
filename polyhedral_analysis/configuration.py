from .atom import Atom
from .utils import flatten

class Configuration:

    def __init__( self, structure, recipes, config_number=None ):
        """
        A Configuration object describes a single atomic geometry.

        Args:
            structure (pymatgen.Structure): A pymatgen Structure object for this configuration.
            recipes (list(PolyhedraRecipe): A list of PolyhedraRecipe objects used to construct
                polyhedra for this configuration.
            config_number (:obj:`int`): An optional integer value to identify this configuration
                in a sequence (e.g. the frame number in a trajectory).
 
        Attributes:
            atoms (list(Atom)): A list of atoms that make up this configuration.
            polyhedra (list(CoordinationPolyhedron): A list of polyhedra, generated using
                the PolyhedraRecipe definitions passed in as the `recipes` list.
            central_atoms (list(Atom)): A list of atoms that define the centres of the 
                coordination polyhedra.
            coordination_atoms (list(Atom)): A list of atoms that define the vertices of 
                the coordination polyhedra.
         """
        self.structure = structure
        self.config_number = config_number
        self.atoms = [ Atom( index=i, site=site, label=site.species_string ) 
                       for i, site in enumerate( self.structure.sites ) ]
        self.polyhedra = []
        for recipe in recipes:
            self.polyhedra.extend( recipe.find_polyhedra( self.atoms, self.structure ) )
        self.central_atoms = sorted( list( set( 
                                 [ p.central_atom for p in self.polyhedra ] ) ) )
        self.coordination_atoms = sorted( list( set( flatten( 
                                      [ p.vertices for p in self.polyhedra ] ) ) ) )

    def coordination_atom_by_index( self, index ):
        return next( atom for atom in self.coordination_atoms if atom.index == index )

    @property
    def polyhedra_labels( self ):
        return [ p.label for p in self.polyhedra ]

    def polyhedra_by_label( self, label ):
        return [ p for p in self.polyhedra if p.label == label ]
