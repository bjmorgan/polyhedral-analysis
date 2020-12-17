from __future__ import annotations

from pymatgen.io.vasp import Xdatcar # type: ignore
from pymatgen.core.structure import Structure
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.utils import flatten
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.configuration import Configuration
import re
import copy
from monty.io import zopen # type: ignore
import pickle
import multiprocessing
from tqdm import tqdm, tqdm_notebook # type: ignore
from functools import partial
from typing import List, Optional, Union, Any, Dict

class Trajectory:

    def __init__(self,
                 structures: List[Structure],
                 configurations: List[Configuration]) -> None:
        """TODO"""
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
                           args: Dict[str, Any]) -> Configuration:
         return Configuration(structure=args['structure'], 
                              recipes=args['recipes'])
    
    @classmethod
    def from_structures(cls,
                        structures: List[Structure],
                        recipes: List[PolyhedraRecipe],
                        verbose: bool = False,
                        ncores: Optional[int] = None,
                        progress: Union[bool, str] = False) -> Trajectory:
        """
        Generate a :obj:`Trajectory` object by applying one or more polyhedral recipes to
        a series of `pymatgen` :obj:`Structure` objects.,

        Args:
            structures (list(:obj:`pymatgen.Structure`)): A list of `pymatgen` 
                :obj:`Structure`` objects.
            recipes (list(PolyhedraRecipe): List of `PolyhedraRecipe` recipes, where each recipe
                defines how to construct a set of `CoordinationPolyhedra` for each configuration.
            verbose (:obj:`bool`, optional): verbose output while parsing the input structures.
                Default is False.
            ncores (:obj:`int`, optional): TODO.
            progress (:obj:`str`, optional): Show a progress bar.
               Setting to ``True`` gives a simple progress bar.
               Setting to ``"notebook"`` gives a Jupyter notebook compatible progress bar.
               Default is ``False``.

        Returns:
            None
        """
        args = [{'structure': structure, 
                 'recipes': recipes} for structure in structures]
        if progress:
            progress_kwargs = {'total': len(args),
                               'unit': ' configurations'}
            if progress == 'notebook':
                pbar = tqdm_notebook
            else:
                pbar = tqdm
        else:
            pbar = iter
        # generate polyhedra configurations
        if ncores:
            with multiprocessing.Pool(ncores) as p:
                configurations = list(pbar(p.imap(cls._get_configuration, args), 
                                                  **progress_kwargs))
        else:
            configurations = list(pbar(map(cls._get_configuration, args),
                                           **progress_kwargs))
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
        """TODO"""
        if not isinstance(other, Trajectory):
            return NotImplemented
        summed_structures = self.structures + other.structures
        summed_configurations = self.configurations + other.configurations
        new_trajectory = Trajectory(summed_structures, summed_configurations)
        return new_trajectory

    @property
    def structures(self) -> List[Structure]:
        return self._structures

    @property
    def configurations(self) -> List[Configuration]:
        return self._configurations

    @classmethod
    def from_xdatcar(cls,
                     filename: str,
                     recipes: List[PolyhedraRecipe],
                     verbose: bool = False,
                     progress: Union[str, bool] = False,
                     ncores: Optional[int] = None) -> Trajectory:
        """
        Args:
            filename (str): Filename for a VASP XDATCAR file.
            recipes (list(:obj:`PolyhedraRecipe`): List of `PolyhedraRecipe` recipes, 
                where each recipe defines how to construct a set of `CoordinationPolyhedra` 
                for each configuration.

        Returns:
            (:obj:`Trajectory`): a :obj:`Trajectory` object.

        """
        xdatcar = Xdatcar( filename )
        structures = xdatcar.structures
        return cls.from_structures(structures, recipes, verbose=verbose, progress=progress,
                                    ncores=ncores)

    @classmethod
    def from_xdatcars(cls,
                      filenames: List[str],
                      recipes: List[PolyhedraRecipe],
                      verbose: bool = False,
                      ncores: Optional[int] = None,
                      progress: Union[bool, str] = False) -> Trajectory:
        """
        TODO
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

