from __future__ import annotations

from fnmatch import fnmatch
from monty.io import zopen  # type: ignore
import json
import os
from pymatgen.core.sites import Site
from pymatgen.core.lattice import Lattice
from typing import List, Dict, Optional, Union, Any
import numpy as np  # type: ignore

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron


class Atom:
    """Atom class"""

    def __init__(self,
                 index: int,
                 site: Site,
                 label: Optional[str] = None) -> None:
        """
        Initialise an Atom object.

        Args:
            index (int): Numerical index identifying this atom.
            site  (pymatgen.Site): The pymatgen Site (or PeriodicSite) object describing this atom.
            label (:obj:`str`, optional): An optional string labelling this atom. If no label
                is given the site species string will be used.

        Attributes:
            in_polyhedra (list): List of polyhedra that this atom is part of.
            neighbours: (dict(int: list)): Dictionary of lists of neighbours for each 
                polyhedron this atom is part of. Dictionary keys are the indexes of 
                each respective polyhedron.

        """
        self.index = index
        self.site = site
        self.label: str = label if label else site.species_string
        self.in_polyhedra: List[CoordinationPolyhedron] = []
        self._neighbours: Dict[int, List[int]] = {}

    @property
    def neighbours(self) -> Dict[int, List[int]]:
        if len(self._neighbours) != len(self.in_polyhedra):
            for p in self.in_polyhedra:
                p.construct_edge_graph()
                p.update_vertex_neighbours()
        return self._neighbours

    def as_dict(self) -> Dict[str, Any]:
        """
        json-serializable :obj:`dict` representation of Atom.
        """
        d = {'index': self.index,
             'site': self.site.as_dict(),  # type: ignore
             'label': self.label,
             'in_polyhedra': [p.index for p in self.in_polyhedra],
             'neighbours': self.neighbours,
             '@module': self.__class__.__module__,
             '@class': self.__class__.__name__}
        return d

    def __repr__(self) -> str:
        return "{} {}".format(self.index, self.site)

    def to(self,
           fmt: Optional[str] = None,
           filename: Optional[str] = None) -> Union[None, str]:
        """
        Outputs the structure to a file or string.

        Args:
            fmt (:obj:`str`, optional): Format to ouput to. Defaults to JSON unlees the filename
                is provided. If fmt is specified this overrides the filename extension
                in the filename.
            filename (:obj:`str`, optional): If provided, the output will be written to a file.
                If `fmt` is not specified, the format will be determined from the filename
                extension.

        Returns:
            (str) if filename is `None`. `None` otherwise.
        """
        if not filename:
            filename = ""
        if not fmt:
            fmt = "json"
        else:
            fmt = fmt.lower()
        fname = os.path.basename(filename)
        if fmt == "json" or fnmatch(fname.lower(), "*.json"):
            s = json.dumps(self.as_dict())
            if filename:
                with zopen(filename, 'wt') as f:
                    f.write('{}'.format(s))
                return None
            else:
                return s
        else:
            raise ValueError(f'Output format "{fmt}" not recognised')

    @property
    def frac_coords(self) -> np.ndarray:
        """Fractional coordinates"""
        return self.site.frac_coords

    @property
    def coords(self) -> np.ndarray:
        """Cartesian coordinates"""
        return self.site.coords

    @property
    def lattice(self) -> Lattice:
        l = self.site.lattice
        assert isinstance(l, Lattice)
        return l

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Atom):
            return NotImplemented
        return self.index < other.index

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Atom):
            return NotImplemented
        return self.index == other.index

    def __hash__(self) -> int:
        return self.index

    def distance(self, other: Atom) -> float:
        """Shortest distance using periodic boundary conditions to another `Atom`.

        Args:
            other (Atom): The other atom.

        Returns:
            float

        """
        d = self.site.distance(other.site)  # type: ignore
        assert isinstance(d, float)
        return d
