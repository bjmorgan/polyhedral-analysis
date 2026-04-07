from functools import partial
import numpy as np # type: ignore
from pymatgen.core import Site, Structure
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.utils import flatten
from polyhedral_analysis.atom import Atom
from collections.abc import Callable, Sequence
from typing import Any, Protocol, cast

class IndexGenerator(Protocol):
    def __call__(self, structure: Structure) -> Sequence[int]: ...

AtomSpec = str | list[str] | list[int] | IndexGenerator

def matching_sites(structure: Structure, 
                   reference_structure: Structure,
                   species: list[str] | None = None) -> list[tuple[Site, int]]:
    """
    Returns a subset of sites from structure (as a list) where each site is the closest to one 
    site in the reference structure.
    
    Args:
        structure (Structure): The structure being analysed.
        reference_structure (Structure): A Structure object containing a set of reference sites.
        species (:obj:`list(str)`, optional): Optional list of strings specifying a subset of species
            to return. Default is None, which specifies all species are returned.
        
    Returns:
        (list[tuple(Site,int)]): A list of length-2 tuples for each matching site. Each tuple
            contains the corresponding pymatgen Site, and the site index.

    """
    matched_sites: list[tuple[Site, int]] = []
    for ref_site in reference_structure:
        dr = [site.distance(ref_site) for site in structure]
        i = int(np.argmin(dr))
        matched_sites.append((structure[i], i))
    if species:
        matched_sites = [(s, i) for s, i in matched_sites if s.species_string in species]
    return matched_sites

def matching_site_indices(structure: Structure,
                          reference_structure: Structure,
                          species: list[str] | None = None) -> list[int]:
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
    return [m[1] for m in matching_sites(structure, reference_structure, species=species)]

def create_matching_site_generator(reference_structure: Structure,
                                   species: list[str] | None = None) -> Callable[[Structure], list[int]]:
    """
    Args:
        reference_structure (Structure):
        species (:obj:`list[str]`, optional):

    Returns:
        (func): 

    """
    return partial(matching_site_indices, reference_structure=reference_structure, species=species)

def generator_from_atom_argument(arg: AtomSpec) -> IndexGenerator:
    """
    Returns a generator function for selecting a subset of sites from a pymatgen :obj:`Structure` object.

    Args:
        arg (various): Argument used to construct the generator function.

    Returns:
        (func): Generator function that takes one argument (:obj:`Structure`) and
            returns an appropriate list of site indices.

    """
    if callable(arg):
        return arg
    if isinstance(arg, str):
        return partial(_get_indices_from_str, arg=arg)
    elif isinstance(arg, (list, tuple)):
        if all(isinstance(item, str) for item in arg):
            return partial(_get_indices_from_list_str, arg=cast(Sequence[str], arg))
        elif all(isinstance(item, int) for item in arg):
            return partial(_get_indices_from_list_int, arg=cast(Sequence[int], arg))
        else:
            raise TypeError("List items must be all strings or all integers")
    else:
        raise TypeError("Argument must be callable, string, or list of strings or integers")

def _get_indices_from_str(structure: Structure,
                          arg: str) -> Sequence[int]:
    return structure.indices_from_symbol(arg)

def _get_indices_from_list_str(structure: Structure, 
                               arg: Sequence[str]) -> Sequence[int]:
    return tuple(flatten([structure.indices_from_symbol(sp) for sp in arg]))

def _get_indices_from_list_int(structure: Structure,
                               arg: Sequence[int]) -> Sequence[int]:
    return arg

class PolyhedraRecipe:

    allowed_methods: list[str] = ['distance cutoff',
                                  'closest centre',
                                  'nearest neighbours']

    def __init__(self,
                 method: str,
                 central_atoms: AtomSpec,
                 vertex_atoms: AtomSpec, 
                 coordination_cutoff: float | None = None,
                 vertex_graph_cutoff: float | None = None,
                 label: str | None = None,
                 n_neighbours: int | None = None,
                 recalculate: bool = True) -> None:
        """
        Create a :obj:`PolyhedraRecipe` object.

        Args:
            method (str): Method used for constructing coordination polyhedra. 
                Implemented options are:
                    - `distance cutoff`: include all coordinating ions within a cutoff distance.
                    - `closest centre`: include all coordinating ions that share a common closest centre atom.
                    - `nearest neighbours`: include n nearest neighbours to each centre atom.
            coordination_cutoff (:obj:`float`, optional): Cutoff distance for vertex atoms to
                be considered as coordinated to a central atom.
            central_atoms (various): Defines the set of atoms considered polyhedron centres.
            vertex_atoms  (various): Defines the set of atoms considered polyhedron vertices.
            vertex_graph_cutoff (:obj:`float`, optional): TODO.
            n_neighbours (:obj:`int`, optional): Optionally set the maximum number of neighbours to include.
            label (str): Label for this recipe.
            recalculate (:obj:`bool`, optional): Reassign central atom and vertex atom indices for each new
                configuration parsed. Default is `True`. If `recalulate` is set to `False` then the lists of
                atom indices for central and vertex atoms will be fixed once initialised.

        Returns: None

        Notes:
            `central_atoms` and `vertex_atoms` define how to select 
            a subset of atoms from a single :obj:`Structure`. How this is implemented depends on
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
            raise ValueError(f"\n'{method}' is not a recognised string for selecting the recipe method.\n Valid methods are:\n{methods_string}")
            
        self.method = method
        self._central_atom_list_generator = generator_from_atom_argument(central_atoms)
        self._vertex_atom_list_generator = generator_from_atom_argument(vertex_atoms)
        self._central_atom_list: list[int] | None = None
        self._vertex_atom_list: list[int] | None = None
        self.coordination_cutoff = coordination_cutoff
        self.vertex_graph_cutoff = vertex_graph_cutoff
        self.n_neighbours = n_neighbours
        self.label = label
        self.recalculate = recalculate

    def central_atom_list(self,
                          structure: Structure | None = None) -> list[int]:
        """Return the list of central atom indices for a given structure.

        Args:
            structure: The structure to select atoms from. Required on
                first call and whenever ``recalculate`` is ``True``.

        Returns:
            List of integer site indices for the central atoms.

        Raises:
            ValueError: If ``structure`` is not provided when needed.
        """
        if self._central_atom_list is None:
            if structure:
                self._central_atom_list = list(self._central_atom_list_generator(structure))
            else:
                raise ValueError('Needs structure argument')
        elif self.recalculate:
            if structure:
                self._central_atom_list = list(self._central_atom_list_generator(structure))
            else:
                raise ValueError('Needs structure argument')
        return self._central_atom_list

    def vertex_atom_list(self,
                         structure: Structure | None = None) -> list[int]:
        """Return the list of vertex atom indices for a given structure.

        Args:
            structure: The structure to select atoms from. Required on
                first call and whenever ``recalculate`` is ``True``.

        Returns:
            List of integer site indices for the vertex atoms.

        Raises:
            ValueError: If ``structure`` is not provided when needed.
        """
        if self._vertex_atom_list is None:
            if structure:
                self._vertex_atom_list = list(self._vertex_atom_list_generator(structure))
            else:
                raise ValueError('Needs structure argument')
        elif self.recalculate:
            if structure:
                self._vertex_atom_list = list(self._vertex_atom_list_generator(structure))
            else:
                raise ValueError('Needs structure argument')
        return self._vertex_atom_list
                                             
    def find_polyhedra(self,
                       atoms: list[Atom],
                       structure: Structure | None = None) -> list[CoordinationPolyhedron]:
        """Construct coordination polyhedra from a list of atoms.

        Applies this recipe's method and atom selection rules to build
        a list of :class:`CoordinationPolyhedron` objects.

        Args:
            atoms: All atoms in the configuration.
            structure: The pymatgen Structure, used to resolve atom
                indices when ``recalculate`` is ``True``.

        Returns:
            A list of CoordinationPolyhedron objects.
        """
        central_atom_list = self.central_atom_list(structure)
        vertex_atom_list = self.vertex_atom_list(structure)
        central_atoms = [atom for atom in atoms if atom.index in central_atom_list]
        vertex_atoms = [atom for atom in atoms if atom.index in vertex_atom_list]
        
        if self.method == 'distance cutoff':
            if self.coordination_cutoff is None:
                raise ValueError("coordination_cutoff must be set for 'distance cutoff' method")
            return polyhedra_from_distance_cutoff(central_atoms=central_atoms, 
                                                  vertex_atoms=vertex_atoms, 
                                                  cutoff=self.coordination_cutoff,
                                                  label=self.label)
        elif self.method == 'nearest neighbours':
            if self.n_neighbours is None:
                raise ValueError("n_neighbours must be set for 'nearest neighbours' method")
            return polyhedra_from_nearest_neighbours(central_atoms=central_atoms, 
                                                     vertex_atoms=vertex_atoms, 
                                                     nn=self.n_neighbours,
                                                     label=self.label)
        elif self.method == 'closest centre':
            return polyhedra_from_closest_centre(central_atoms=central_atoms,
                                                 vertex_atoms=vertex_atoms,
                                                 label=self.label)
        else:
            raise ValueError(f"Unsupported method: {self.method}")

def polyhedra_from_distance_cutoff(central_atoms: list[Atom],
                                   vertex_atoms: list[Atom],
                                   cutoff: float,
                                   label: str | None = None) -> list[CoordinationPolyhedron]:
    """Construct polyhedra by including all vertex atoms within a cutoff distance.

    Args:
        central_atoms: Atoms to use as polyhedron centres.
        vertex_atoms: Candidate vertex atoms.
        cutoff: Maximum distance for a vertex atom to be included.
        label: Optional label for the resulting polyhedra.

    Returns:
        A list of CoordinationPolyhedron objects.
    """
    if not central_atoms:
        return []
    polyhedra = []
    lattice = central_atoms[0].site.lattice
    distance_matrix = lattice.get_all_distances([c.frac_coords for c in central_atoms],
                                                [v.frac_coords for v in vertex_atoms])
    for dr, c_atom in zip(distance_matrix, central_atoms):
        indices = np.where(dr <= cutoff)[0]
        vertices = [vertex_atoms[i] for i in indices]
        polyhedra.append(CoordinationPolyhedron(central_atom=c_atom, 
                                                vertices=vertices, 
                                                label=label))
    return polyhedra

def polyhedra_from_nearest_neighbours(central_atoms: list[Atom],
                                      vertex_atoms: list[Atom],
                                      nn: int,
                                      label: str | None = None) -> list[CoordinationPolyhedron]:
    """Construct polyhedra from the *n* nearest vertex atoms to each centre.

    Args:
        central_atoms: Atoms to use as polyhedron centres.
        vertex_atoms: Candidate vertex atoms.
        nn: Number of nearest neighbours to include.
        label: Optional label for the resulting polyhedra.

    Returns:
        A list of CoordinationPolyhedron objects.
    """
    polyhedra = []
    for c_atom in central_atoms:
        vertices = sorted(vertex_atoms, key=lambda atom: atom.distance(c_atom))[:nn]
        polyhedra.append(CoordinationPolyhedron(central_atom=c_atom,
                                                vertices=vertices,
                                                label=label))
    return polyhedra

def polyhedra_from_closest_centre(central_atoms: list[Atom],
                                  vertex_atoms: list[Atom],
                                  label: str | None = None) -> list[CoordinationPolyhedron]:
    """Construct polyhedra by assigning each vertex atom to its closest centre.

    Each vertex atom is assigned to the central atom it is nearest to,
    so every vertex belongs to exactly one polyhedron.

    Args:
        central_atoms: Atoms to use as polyhedron centres.
        vertex_atoms: Candidate vertex atoms.
        label: Optional label for the resulting polyhedra.

    Returns:
        A list of CoordinationPolyhedron objects.
    """
    coordination_coords = [co_atom.coords for co_atom in vertex_atoms]
    central_coords = [c_atom.coords for c_atom in central_atoms]
    closest_site_index = [np.argmin([a.distance(c_atom) 
                                     for c_atom in central_atoms]) 
                                     for a in vertex_atoms]
    polyhedra = []
    for i, c_atom in enumerate(central_atoms):
        vertices = [co for c, co in zip( closest_site_index, vertex_atoms) if c == i]
        polyhedra.append(CoordinationPolyhedron(central_atom=c_atom, 
                                                vertices=vertices, 
                                                label=label))
    return polyhedra

def polyhedra_from_atom_indices(central_atoms: list[Atom],
                                vertex_atoms: list[Atom],
                                central_indices: list[int],
                                vertex_indices: list[list[int]],
                                label: str | None = None) -> list[CoordinationPolyhedron]:
    """Construct a set of polyhedra from lists of atom indices for central and vertex atoms.
c	

    Args:
        central_atoms (list(Atom)): List of Atom objects describing the set of possible centre atoms.
        vertex_atoms (list(Atom)): List of Atom objects describing the set of possible vertex atoms.
        central_indices (list(int)): List of integer indices specifying each central atom.
        vertex_indices (list(list(int)): Nested list of integer indices for the vertex atoms in each polyhedron.

    Returns:
        list(CoordinationPolyhedron)

    Raises:
        ValueError: If the lengths of the central_indices and vertex_indices lists are unequal.

    """
    if len(central_indices) != len(vertex_indices):
        raise ValueError('central_indices and vertex_indices are different lengths: '
                         f'{len(central_indices)} vs. {len(vertex_indices)}.')
    central_atom_map = {atom.index: atom for atom in central_atoms}
    vertex_atom_map = {atom.index: atom for atom in vertex_atoms}
    polyhedra = []
    for ic, iv in zip(central_indices, vertex_indices):
        try:
            central_atom = central_atom_map[ic]
        except KeyError:
            raise ValueError(
                f'Central atom index {ic} not found in central_atoms.') from None
        try:
            vertices = [vertex_atom_map[i] for i in iv]
        except KeyError:
            missing = [i for i in iv if i not in vertex_atom_map]
            raise ValueError(
                f'Vertex atom indices {missing} not found in vertex_atoms.') from None
        polyhedra.append(CoordinationPolyhedron(central_atom=central_atom,
                                                vertices=vertices,
                                                label=label))
    return polyhedra 

