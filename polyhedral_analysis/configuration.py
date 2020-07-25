from polyhedral_analysis.atom import Atom
from polyhedral_analysis.utils import flatten
from pymatgen.core.structure import Structure
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from typing import Optional, List, Union, Dict

class Configuration:

    def __init__(self,
                 structure: Structure,
                 recipes: List[PolyhedraRecipe]) -> None:
        """
        A Configuration object describes a single atomic geometry.

        Args:
            structure (pymatgen.Structure): A pymatgen Structure object for this configuration.
            recipes (list(PolyhedraRecipe): A list of PolyhedraRecipe objects used to construct
                polyhedra for this configuration.

        Attributes:
            atoms (list(Atom)): A list of atoms that make up this configuration.
            polyhedra (list(CoordinationPolyhedron): A list of polyhedra, generated using
                the PolyhedraRecipe definitions passed in as the `recipes` list.
            central_atoms (list(Atom)): A list of atoms that define the centres of the 
                coordination polyhedra.
            coordination_atoms (list(Atom)): A list of atoms that define the vertices of 
                the coordination polyhedra.
         """
        self.atoms = [Atom(index=i,
                           site=site,
                           label=site.species_string)
                      for i, site in enumerate(structure.sites)]
        self.polyhedra: List[CoordinationPolyhedron] = []
        for recipe in recipes:
            self.polyhedra.extend(recipe.find_polyhedra(self.atoms, structure))
        self.central_atoms = sorted(list(set(
            [p.central_atom for p in self.polyhedra])),
            key=lambda x: x.index)
        self.coordination_atoms = sorted(list(set(flatten(
            [p.vertices for p in self.polyhedra]))),
            key=lambda x: x.index)

    def coordination_atom_by_index(self,
                                   index: int) -> Union[Atom, None]:
        """
        Return the coordination atom with a specific index.

        Args:
            index (int): The atom index to match.

        Returns:
            (Atom|None): The matching coordination atom. If the desired index does not match
                         any of the coordination atoms for this configuration, None is returned.
        """
        coordination_atom_indices = [
            atom.index for atom in self.coordination_atoms]
        if index not in coordination_atom_indices:
            return None
        else:
            return self.coordination_atoms[coordination_atom_indices.index(index)]

    @property
    def polyhedra_labels(self) -> List[Union[str, None]]:
        return [p.label for p in self.polyhedra]

    def polyhedra_by_label(self,
                           label: Union[str, List[str]]) -> List[CoordinationPolyhedron]:
        """Returns a list of polyhedra for this configuration with matching labels.

        Args:
            (:obj:`str` or :obj:`list(str)`): Either a single label string, or a list of
                label strings.

        Returns:
            (list(CoordinationPolyhedron)

        """
        if isinstance(label, str):
            return [p for p in self.polyhedra if p.label == label]
        elif isinstance(label, list):
            return [p for p in self.polyhedra if p.label in label]
        else:
            raise TypeError('Invalid type for label argument')

    def neighbour_list(self,
                       label1: Optional[str] = None,
                       label2: Optional[str] = None) -> Dict[int, List[int]]:
        """
        TODO
        """
        pass 
