from polyhedral_analysis.atom import Atom
from polyhedral_analysis.utils import flatten, prune_neighbour_list, lattice_mc_string
from pymatgen.core.structure import Structure
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from typing import List, Union, Dict, Tuple

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
            polyhedra (list(CoordinationPolyhedron)): A list of polyhedra, generated using
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
        """Return the coordination atom with a specific index.

        Args:
            index (int): The atom index to match.

        Returns:
            Union[Atom, None]: The matching coordination atom. If the desired index does not match
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
            (Union[str, List[str]]): Either a single label string, or a list of
                label strings.

        Returns:
            list(CoordinationPolyhedron)

        """
        if isinstance(label, str):
            return [p for p in self.polyhedra if p.label == label]
        elif isinstance(label, list):
            return [p for p in self.polyhedra if p.label in label]
        else:
            raise TypeError('Invalid type for label argument')

    def to_lattice_mc(self, 
                      filename: str, 
                      labels: List[str],
                      neighbour_list: Dict[int, Tuple[int, ...]]) -> None:
        site_indices = [p.index for p in self.polyhedra_by_label(labels)]
        neighbour_list = prune_neighbour_list(neighbour_list, site_indices)
        with open(filename, 'w') as f:
            n_sites = len(self.polyhedra_by_label(labels))
            f.write(f'{n_sites}\n\n')
            for l in labels:
                for p in self.polyhedra_by_label(l):
                    f.write(f'{lattice_mc_string(p, neighbour_list)}\n')

    def face_sharing_neighbour_list(self, 
                                    labels: List[str]) -> Dict[int, Tuple[int, ...]]:
        return {p.index: p.face_sharing_neighbour_list()
                for p in self.polyhedra_by_label(labels)}
