from __future__ import annotations

import copy
import multiprocessing
from collections.abc import Iterable
from typing import Any

from pymatgen.io.vasp import Xdatcar
from pymatgen.core.structure import Structure
from tqdm.auto import tqdm

from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.utils import flatten

class Trajectory:

    def __init__(self,
                 structures: list[Structure],
                 configurations: list[Configuration]) -> None:
        """Initialise a Trajectory from pre-built structures and configurations.

        Args:
            structures: A list of pymatgen Structure objects, one per frame.
            configurations: A list of Configuration objects corresponding
                to each structure.

        Raises:
            TypeError: If any element of ``structures`` is not a
                :class:`~pymatgen.core.structure.Structure`, or any element
                of ``configurations`` is not a
                :class:`~polyhedral_analysis.configuration.Configuration`.
        """
        for s in structures:
            if not isinstance( s, Structure ):
                raise TypeError( "structures argument should contain pymatgen Structure objects" )
        for c in configurations:
            if not isinstance( c, Configuration ):
                raise TypeError( "configurations argument should contain Configuration objects" )
        self._structures = structures
        self._configurations = configurations

    @classmethod
    def _get_configuration(cls,
                           args: dict[str, Any]) -> Configuration:
         return Configuration(structure=args['structure'], 
                              recipes=args['recipes'])
    
    @classmethod
    def from_structures(cls,
                        structures: list[Structure],
                        recipes: list[PolyhedraRecipe],
                        verbose: bool = False,
                        ncores: int | None = None,
                        progress: bool = False) -> Trajectory:
        """Generate a Trajectory by applying recipes to a series of structures.

        Args:
            structures: A list of pymatgen Structure objects.
            recipes: List of PolyhedraRecipe recipes, where each recipe
                defines how to construct a set of CoordinationPolyhedra
                for each configuration.
            verbose: Verbose output while parsing the input structures.
                Default is ``False``.
            ncores: Number of cores for parallel processing. Default
                is ``None`` (serial).
            progress: Show a progress bar. Uses ``tqdm.auto`` so the
                correct widget is chosen automatically for terminals
                and Jupyter notebooks. Default is ``False``.

        Returns:
            A Trajectory object.
        """
        args = [{'structure': structure,
                 'recipes': recipes} for structure in structures]
        iterable: Iterable[Configuration]
        if ncores:
            with multiprocessing.Pool(ncores) as p:
                iterable = p.imap(cls._get_configuration, args)
                if progress:
                    iterable = tqdm(iterable, total=len(args), unit=' configs')
                configurations = list(iterable)
        else:
            iterable = map(cls._get_configuration, args)
            if progress:
                iterable = tqdm(iterable, total=len(args), unit=' configs')
            configurations = list(iterable)
        return cls(structures=structures,
                   configurations=configurations)

    def extend(self, 
               other: Trajectory,
               offset: int = 0) -> None:
        """
        Extend this Trajectory with the data from another Trajectory.

        Args:
            other (`Trajectory`): The trajectory to be appended.
            offset (int): Offset to apply to the frame numbers.

        Returns:
            None
        """
        extended_configurations = copy.deepcopy(other.configurations)
        self.configurations.extend(extended_configurations)

    def __add__(self, other: object) -> Trajectory:
        """Concatenate two trajectories, returning a new Trajectory.

        Args:
            other: Another Trajectory to concatenate with this one.

        Returns:
            A new Trajectory containing the structures and configurations
            from both trajectories.
        """
        if not isinstance(other, Trajectory):
            return NotImplemented
        summed_structures = self.structures + other.structures
        summed_configurations = self.configurations + other.configurations
        new_trajectory = Trajectory(summed_structures, summed_configurations)
        return new_trajectory

    @property
    def structures(self) -> list[Structure]:
        return self._structures

    @property
    def configurations(self) -> list[Configuration]:
        return self._configurations

    @classmethod
    def from_xdatcar(cls,
                     filename: str,
                     recipes: list[PolyhedraRecipe],
                     verbose: bool = False,
                     progress: bool = False,
                     ncores: int | None = None) -> Trajectory:
        """Generate a Trajectory from a VASP XDATCAR file.

        Args:
            filename: Path to a VASP XDATCAR file.
            recipes: List of PolyhedraRecipe recipes defining how to
                construct coordination polyhedra for each configuration.
            verbose: Verbose output while parsing. Default is ``False``.
            progress: Show a progress bar. Default is ``False``.
            ncores: Number of cores for parallel processing. Default
                is ``None`` (serial).

        Returns:
            A Trajectory object.
        """
        xdatcar = Xdatcar( filename )
        structures = xdatcar.structures
        return cls.from_structures(structures, recipes, verbose=verbose, progress=progress,
                                    ncores=ncores)

    @classmethod
    def from_xdatcars(cls,
                      filenames: list[str],
                      recipes: list[PolyhedraRecipe],
                      verbose: bool = False,
                      ncores: int | None = None,
                      progress: bool = False) -> Trajectory:
        """Generate a Trajectory from multiple VASP XDATCAR files.

        This is a convenience method for loading trajectories split across
        several XDATCAR files (e.g. from consecutive VASP runs). The
        structures from all files are concatenated before building
        configurations.

        Args:
            filenames: List of paths to VASP XDATCAR files.
            recipes: List of PolyhedraRecipe recipes defining how to
                construct coordination polyhedra for each configuration.
            verbose: Verbose output while parsing. Default is ``False``.
            ncores: Number of cores for parallel processing. If set,
                both XDATCAR parsing and configuration building are
                parallelised. Default is ``None`` (serial).
            progress: Show a progress bar. Default is ``False``.

        Returns:
            A Trajectory object.
        """
        if ncores and len(filenames) > 1:
            with multiprocessing.Pool(ncores) as p:
                xdatcars = p.map(_get_xdatcar, filenames)
        else:
            xdatcars = [Xdatcar(f) for f in filenames]
        structures = flatten([x.structures for x in xdatcars])
        return cls.from_structures(structures,
                                   recipes,
                                   verbose,
                                   ncores=ncores,
                                   progress=progress) 

    def __len__(self) -> int:
        return len(self.structures)

def _get_xdatcar(filename: str) -> Xdatcar:
    """
    Internal method to support multiprocessing.
    """
    return Xdatcar(filename)

