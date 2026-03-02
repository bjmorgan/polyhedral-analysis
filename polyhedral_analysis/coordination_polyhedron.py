from __future__ import annotations

from polyhedral_analysis.atom import Atom
from polyhedral_analysis.symmetry_measure import symmetry_measures_from_coordination
from polyhedral_analysis.orientation_parameters import cos_theta
from pymatgen.core.sites import Site
from pymatgen.util.coord import pbc_shortest_vectors
import numpy as np 
from scipy.spatial import ConvexHull  # type: ignore[import-untyped]
from itertools import combinations
import vg  # type: ignore
from collections import Counter
from typing import Literal
import typing

class CoordinationPolyhedron:
    """A coordination polyhedron defined by a central atom and its vertex atoms.

    Provides geometric analysis including bond distances, angles, volume,
    edge connectivity, continuous symmetry measures, and neighbour-sharing
    relationships with other polyhedra.
    """

    def __init__(self,
                 central_atom: Atom,
                 vertices: list[Atom],
                 label: str | None = None) -> None:
        """
        Initialise a CoordinationPolyhedron object.

        Args:
            central_atom (Atom): the central atom.
            vertices (list(Atom)): A list of atoms that define the coordination environment.
            label (:obj:`str`, optional): An optional string used to label this coordination polyhedron.
                if the label is not defined, the label of the central atom will be used.

        Returns:
            None
        """
        self.central_atom = central_atom
        self.central_atom.in_polyhedra.append(self)
        self.vertices = sorted(vertices)
        for v in self.vertices:
            v.in_polyhedra.append(self)
        self._label: str = label if label else central_atom.label
        self._edge_graph: dict[int, list[int]] | None = None
        self._symmetry_measure: dict[str, float] | None = None

    @property
    def label(self) -> str:
        """The label for this polyhedron."""
        return self._label

    def set_label(self, label: str) -> None:
        """Set the label for this polyhedron.

        Args:
            label: The new label string.
        """
        self._label = label

    def update_vertex_neighbours(self) -> None:
        """Update the neighbour lists on each vertex atom from the edge graph."""
        for vertex, neighbour_list in zip(self.vertices, self.edge_graph.values()):
            vertex._neighbours[self.index] = neighbour_list

    def __repr__(self) -> str:
        """
        String representation of a Polyhedron object.

        The output includes the polyhedron label (if this is set), and
        information about the central and vertex atoms. For each atom
        the output includes::

            atom_index [ x y z ] atom_species

        Examples:

            >>> print(polyhedron)
            Coordination Polyhedron 4c
            255 [12.71362322 17.90999634 12.74490767] S
            ----------
            31 [12.46919306 20.2317206  12.2641591 ] Li
            55 [13.0016308  17.39863735 10.46318072] Li
            71 [10.4034848  18.18407515 12.43873978] Li
            103 [12.17924193 15.66932958 13.34077502] Li
            159 [13.24242002 18.43469275 15.02193658] Li
            175 [15.02830461 17.60091516 12.52079631] Li

        """
        if self.label:
            to_return = 'Coordination Polyhedron {}\n'.format(self.label)
        else:
            to_return = 'Coordination Polyhedron\n'
        to_return += '{}\n'.format(self.central_atom)
        to_return += '----------\n'
        for v in self.vertices:
            to_return += '{}\n'.format(v)
        return to_return

    def intersection(self,
                     other_polyhedron: CoordinationPolyhedron) -> tuple[int, ...]:
        """Returns a tuple of atom indices for vertex atoms shared with another polyhedron.

        Args:
            other_polyhedron (:obj:`CoordinationPolyhedron`): The other coordination polyhedron.

        Returns:
            tuple(int): Tuple of shared vertex indices.

        """
        return tuple(sorted(set(self.vertex_indices) & set(other_polyhedron.vertex_indices)))

    @property
    def vertex_indices(self) -> list[int]:
        """Global indices of the vertex atoms."""
        return [v.index for v in self.vertices]

    def vertex_vectors(self,
                       reference: str='centroid') -> np.ndarray:
        """Returns an array of polyhedron centre -> vertex vectors.
        The polyhedron centre can be set as either the centroid of the polyhedron (default)
        or the position of the central atom.

        Args:
            reference (str, optional): The reference point used to compute the vertex vectors.
                Can either be 'centroid' (default) or 'central_atom'.

        Returns:
            np.ndarray: Array of vectors from the reference point to each vertex.

        Raises:
            ValueError: If an invalid reference point setting is passed.

        """
        if reference == 'centroid':
            reference_point = self.centroid()
        elif reference == 'central_atom':
            reference_point = self.central_atom.coords
        else:
            raise ValueError("Invalid reference point. Use 'centroid' or 'central_atom'.")
        return self.minimum_image_vertex_coordinates() - reference_point

    @property
    def vertex_coords(self) -> np.ndarray:
        """Cartesian coordinates of the vertex atoms as an Nx3 array."""
        return np.array([v.coords for v in self.vertices])

    @property
    def vertex_labels(self) -> list[str | None]:
        """Species labels of the vertex atoms."""
        return [v.label for v in self.vertices]

    @property
    def coordination_number(self) -> int:
        """Number of vertex atoms in this polyhedron."""
        return len(self.vertices)

    @property
    def index(self) -> int:
        """Global index of the central atom."""
        return self.central_atom.index

    @property
    def edge_graph(self) -> dict[int, list[int]]:
        """Edge connectivity graph for this polyhedron.

        A dictionary mapping each vertex index to a list of indices of
        vertices connected to it by an edge. Computed lazily from the
        convex hull on first access.
        """
        if self._edge_graph is None:
            self._edge_graph = self.construct_edge_graph()
            self.update_vertex_neighbours()
        return {k: list(v) for k, v in self._edge_graph.items()}

    def edge_vertex_indices(self) -> tuple[tuple[int, int], ...]:
        """Return all edges as sorted pairs of vertex indices."""
        edge_pairs: set[tuple[int, int]] = set()
        for v1, v2_list in self.edge_graph.items():
            for v2 in v2_list:
                edge = (min(v1, v2), max(v1, v2))
                edge_pairs.add(edge)
        return tuple(sorted(edge_pairs))

    @property
    def symmetry_measure(self) -> dict[str, float]:
        """Continuous symmetry measures against all reference geometries for this coordination number.

        Returns a dictionary mapping geometry names to CSM values.
        Computed lazily on first access.

        Raises:
            ValueError: If no reference geometries are defined for this
                coordination number.
        """
        if self._symmetry_measure is None:
            if self.coordination_number not in symmetry_measures_from_coordination:
                raise ValueError('No symmetry measure objects for coordination number of {}'.format(
                    self.coordination_number))
            msm = {}
            for string, sm in symmetry_measures_from_coordination[self.coordination_number].items():
                msm[string] = sm.minimum_symmetry_measure(self.vertex_vectors(reference='central_atom'))
            self._symmetry_measure = msm
        return dict(self._symmetry_measure)

    @property
    def best_fit_geometry(self) -> dict[str, str | float]:
        """The reference geometry with the lowest symmetry measure.

        Returns a dictionary with keys ``'geometry'`` (the name) and
        ``'symmetry_measure'`` (the CSM value).
        """
        psm = self.symmetry_measure
        best_fit = min(psm, key=psm.get) # type: ignore
        return {'geometry': best_fit, 'symmetry_measure': psm[best_fit]}

    def minimum_image_vertex_coordinates(self) -> np.ndarray:
        """Vertex coordinates under the minimum-image convention.

        Returns Cartesian coordinates for each vertex, shifted to be the
        closest periodic image to the central atom.

        Returns:
            np.ndarray: An Nx3 array of minimum-image vertex coordinates.
        """
        vertex_frac_coords = [v.frac_coords for v in self.vertices]
        pbc_vectors = pbc_shortest_vectors(self.central_atom.lattice,
                                           self.central_atom.frac_coords,
                                           vertex_frac_coords)[0]
        vertex_minimum_image_coords = [
            self.central_atom.coords + v for v in pbc_vectors]
        return np.array(vertex_minimum_image_coords)

    def faces(self) -> tuple[tuple[int, ...], ...]:
        """
        Returns a nested set of tuples for each face of the polyhedron. Each 
        per-face tuple returned contains the indices of the vertices that define
        that face, in sorted numerical order.

        Args:
            None

        Returns:
            tuple[tuple[int]]

        Examples: 

            >>> print(polyhedron)
            Coordination Polyhedron 4c
            255 [12.71362322 17.90999634 12.74490767] S
            ----------
            31 [12.46919306 20.2317206  12.2641591 ] Li
            55 [13.0016308  17.39863735 10.46318072] Li
            71 [10.4034848  18.18407515 12.43873978] Li
            103 [12.17924193 15.66932958 13.34077502] Li
            159 [13.24242002 18.43469275 15.02193658] Li
            175 [15.02830461 17.60091516 12.52079631] Li
            >>> polyhedron.faces()
            ((31, 159, 175), (31, 55, 71), (31, 55, 175), (31, 71, 159)
             (55, 71, 103), (55, 103, 175), (71, 103, 159), (103, 159, 175))

        """
        return tuple(
            tuple(
                sorted([self.vertex_indices[v] for v in simplex])
            ) for simplex in merge_coplanar_simplices(self.convex_hull())
        )

    def convex_hull(self) -> ConvexHull:
        """Compute the convex hull of the minimum-image vertex coordinates."""
        return ConvexHull(self.minimum_image_vertex_coordinates())

    def construct_edge_graph(self) -> dict[int, list[int]]:
        """Build the edge connectivity graph from the convex hull.

        For polyhedra with more than three vertices, edges are derived
        from the convex hull simplices, merging coplanar faces. For
        three or fewer vertices, all vertices are connected.

        Returns:
            A dictionary mapping each vertex index to the indices of
            its edge-connected neighbours.
        """
        connected_vertices: dict[int, set[int]] = {
            i: set() for i in range(self.coordination_number)}
        if self.coordination_number > 3:
            convex_hull = self.convex_hull()
            for m in merge_coplanar_simplices(convex_hull):
                if len(m) == 3:
                    for r in range(3):
                        rotated_simplex = np.roll(m, r)
                        connected_vertices[rotated_simplex[0]].add(
                            rotated_simplex[1])
                        connected_vertices[rotated_simplex[0]].add(
                            rotated_simplex[2])
                # non-triangular face with > 4 vertices. This will be a composite of more than one simplex.
                else:
                    component_simplices = []
                    for s in convex_hull.simplices:
                        if np.all([i in m for i in s]):
                            component_simplices.append(s)
                    # common elements are linked along an internal edge.
                    for s_roll in range(len(component_simplices)):
                        rotated_component_simplices = np.roll(
                            component_simplices, s_roll, axis=0)
                        this_simplex = rotated_component_simplices[0]
                        other_simplices = rotated_component_simplices[1:]
                        for roll in range(3):
                            rotated_simplex = np.roll(this_simplex, roll)
                            edge = [rotated_simplex[0], rotated_simplex[1]]
                            if not np.all([i in np.unique(other_simplices) for i in edge]):
                                connected_vertices[edge[0]].add(edge[1])
                            edge = [rotated_simplex[0], rotated_simplex[2]]
                            if not np.all([i in np.unique(other_simplices) for i in edge]):
                                connected_vertices[edge[0]].add(edge[1])
        else:
            for roll in range(self.coordination_number):
                rotated_list = np.roll(
                    list(range(self.coordination_number)), roll)
                for i in rotated_list[1:]:
                    connected_vertices[rotated_list[0]].add(i)
        edge_list = {}
        for i in range(self.coordination_number):
            edge_list[self.vertex_indices[i]] = [self.vertex_indices[v]
                                                 for v in connected_vertices[i]]
        return edge_list

    @property
    def vertex_count(self) -> typing.Counter[str | None]:
        """Count of vertex atoms grouped by species label."""
        return Counter([v.label for v in self.vertices])

    def vertex_distances(self,
                         reference: Literal['central_atom', 'centroid'] = 'central_atom'
                        ) -> np.ndarray:
        """
        Returns an array of distances from either the central atom or the centroid to the vertex atoms.
 
        Args:
            reference (str, optional): The reference point for distance calculations. 
                Can be either 'central_atom' (default) or 'centroid'.

        Returns:
            np.ndarray[float, ...]: An array of distances between each vertex and the reference point.

        Raises:
            ValueError: If an invalid reference point is provided.

        """
        if reference not in ['central_atom', 'centroid']:
            raise ValueError("Invalid reference point. Use 'central_atom' or 'centroid'.")
        vectors = self.vertex_vectors(reference=reference)
        distances = np.linalg.norm(vectors, axis=1)
        return distances

    def vertex_distances_and_labels(self,
                                    reference: Literal['central_atom', 'centroid'] = 'central_atom'
                                    ) -> tuple[tuple[float, str | None], ...]:
        """
        Returns a tuple of distances and species labels from the reference point to the vertex atoms.

        Args:
            reference (str, optional): The reference point for distance calculations. 
                Can be either 'central_atom' (default) or 'centroid'.

        Returns:
            tuple[tuple[float, str | None], ...]: A tuple of length-2 tuples, containing for each vertex the
                distance from the reference point and the species label of the vertex atom.

        Raises:
            ValueError: If an invalid reference point is provided.

        """
        if reference not in ['central_atom', 'centroid']:
            raise ValueError("Invalid reference point. Use 'central_atom' or 'centroid'.")
    
        return tuple(zip(self.vertex_distances(reference=reference),
                         self.vertex_labels))

    def equal_vertices(self, other: object) -> bool:
        """
        Test whether this :obj:`CoordinationPolyhedron` has vertices with the same labels as
        another :obj:`CoordinationPolyhedron`.

        Args:
            other (:obj:`CoordinationPolyhedron`): The other :obj:`CoordinationPolyhedron`.

        Returns:
            bool: True / False.

        """
        if not isinstance(other, CoordinationPolyhedron):
            raise TypeError
        return self.vertex_indices == other.vertex_indices

    def equal_edge_graph(self, other: object) -> bool:
        """
        Test whether this :obj:`CoordinationPolyhedron` has the same edge graph as
        another :obj:`CoordinationPolyhedron`.

        Args:
            other (:obj:`CoordinationPolyhedron`): The other :obj:`CoordinationPolyhedron`.

        Returns:
            bool: True or False.

        """
        if not isinstance(other, CoordinationPolyhedron):
            raise TypeError
        return self.edge_graph == other._edge_graph

    def equal_members(self, other: object) -> bool:
        """
        Test whether this :obj:`CoordinationPolyhedron` has the same member atoms
        as another :obj:`CoordinationPolyhedron`.

        Args:
            other (:obj:`CoordinationPolyhedron`): The other :obj:`CoordinationPolyhedron`.

        Returns:
            bool: True or False.

        """
        if not isinstance(other, CoordinationPolyhedron):
            raise TypeError
        equal_central_atom = self.central_atom == other.central_atom
        equal_vertex_atoms = self.vertices == other.vertices
        return equal_central_atom & equal_vertex_atoms

    def neighbours(self) -> tuple[CoordinationPolyhedron, ...]:
        """Returns a tuple of neighbouring polyhedra.
        Two polyhedra are considered to be neighbours if they 
        have one or more vertex atoms in common.

        Args:
            none

        Returns:
            tuple(CoordinationPolyhedron): Tuple of neighbouring polyhedra.

        """
        neighbours = []
        for v in self.vertices:
            for p in v.in_polyhedra:
                if p.index != self.index and p not in neighbours:
                    neighbours.append(p)
        return tuple(sorted(neighbours, key=lambda p: p.index))

    def neighbours_by_index_and_shared_vertices(self) -> dict[int, tuple[int, ...]]:
        """Returns a dictionary of vertex atoms shared with neighbouring polyhedra.
        Two polyhedra are considered to be neighbours if they have one of more
        vertex atoms in common.

        Args:
            None

        Returns:
            dict(int, tuple(int, ...)): Dictionary where the keys are polyhedra
                indices, and the values are the indices of common vertex atoms.

        """
        return {p.index: self.intersection(p) for p in self.neighbours()}

    def corner_sharing_neighbour_list(self) -> tuple[int, ...]:
        """Indices of neighbouring polyhedra that share exactly one vertex (corner-sharing)."""
        return tuple(k for k, v in self.neighbours_by_index_and_shared_vertices().items()
                     if len(v) == 1)

    def edge_sharing_neighbour_list(self) -> tuple[int, ...]:
        """Indices of neighbouring polyhedra that share exactly two vertices (edge-sharing)."""
        return tuple(k for k, v in self.neighbours_by_index_and_shared_vertices().items()
                     if len(v) == 2)

    def face_sharing_neighbour_list(self) -> tuple[int, ...]:
        """Indices of neighbouring polyhedra that share three or more vertices (face-sharing)."""
        return tuple(k for k, v in self.neighbours_by_index_and_shared_vertices().items()
                     if len(v) >= 3)

    def __eq__(self, other: object) -> bool:
        """
        Two :obj:`CoordinationPolyhedron` objects are considered equal if they
        have equal edge graphs.

        Args:
            other (:obj:`CoordinationPolyhedron`): The other :obj:`CoordinationPolyhedron`.

        Returns:
            bool: True or False.

        """
        if not isinstance(other, CoordinationPolyhedron):
            return NotImplemented
        return self.equal_edge_graph(other)

    def vertex_vector_projections(self, vectors: np.ndarray, reference: str = 'centroid') -> np.ndarray:
        """
        Calculate the projection of each vertex vector on one or more input vectors.

        Args:
            vectors (np.ndarray): A Nx3 numpy array, where each row is a vector used
                                  to calculate the projection. `vectors` can also be
                                  a single length 3 array.
            reference (str, optional): The reference point for vertex vectors. 
                                       Can be either 'centroid' (default) or 'central_atom'.

        Returns:
            np.ndarray: A (N_vertex x N_vector) dimension numpy array.

        """
        vectors = np.atleast_2d(vectors)
        vertex_vecs = self.vertex_vectors(reference=reference)
        vectors_normalized = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
        projections = np.dot(vertex_vecs, vectors_normalized.T)
        return projections

    def coordination_distances(self) -> list[float]:
        """
        Returns a list of distances from the central atom to the vertex atoms.

        This method is a wrapper around vertex_distances() that maintains backward compatibility
        and returns a list instead of a tuple.

        Args:
            None

        Returns:
            list[float]: A list of distances between each vertex and the central atom.
        """
        return list(self.vertex_distances(reference='central_atom'))

    def angles(self) -> list[float]:
        """
        List of all vertex--centroid--vertex angles.

        Args:
            None

        Returns:
             list: List of angles.

        """
        return [vg.angle(p1, p2) for p1, p2 in combinations(self.vertex_vectors(reference='centroid'), 2)]

    def vertex_vector_orientations(self,
                                   units: Literal['degrees', 'radians'] = 'degrees',
                                   return_distance: bool = False,
                                   reference: Literal['centroid', 'central_atom'] = 'centroid') -> list[tuple[float, ...]]:
        """Returns the angular orientations of each vertex vector.
    
        The orientation is defined by two angles, theta and phi. 
        Theta is the angle with respect to [0, 0, 1] and ranges from 0 to 180 degrees. 
        Phi is the angle with respect to [1, 0, 0] and ranges from -180 to +180 degrees.
    
        Args:
            units (:obj:`str`, optional): Optionally select the units for the calculated angles.
                                          Options are `degrees` or `radians`. 
                                          Default is `degrees`.
            return_distance (:obj:`bool`, optional): Optionally also return the distance. 
            reference (:obj:`str`, optional): The reference point for vertex vectors. 
                                              Options are `centroid` or `central_atom`.
                                              Default is `centroid`.
    
        Returns:
            list(tuple): A list of `(theta,phi)` tuple pairs, or `(theta,phi,distance)` if
                         `return_distance` is True.
    
        Raises:
            ValueError: If an invalid reference point is provided.
        """
        if reference not in ['centroid', 'central_atom']:
            raise ValueError("Invalid reference point. Use 'centroid' or 'central_atom'.")
        
        vertex_vectors = self.vertex_vectors(reference=reference)
        vg_units = {'degrees': 'deg',
                    'radians': 'rad'}
        
        result: list[tuple[float, ...]] = []
        for point in vertex_vectors:
            theta = vg.angle(np.array([0.0, 0.0, 1.0]), point, units=vg_units[units])
            
            # Handle the edge case where theta is 0 or 180 degrees (or equivalent in radians)
            if np.isclose(theta, 0) or np.isclose(theta, np.pi if units == 'radians' else 180):
                phi = 0.0  # phi is undefined in this case, so we set it to 0
            else:
                phi = vg.signed_angle(np.array([1.0, 0.0, 0.0]), point,
                                      look=np.array([0.0, 0.0, 1.0]), units=vg_units[units])
            
            if return_distance:
                distance = vg.magnitude(point)
                result.append((theta, phi, distance))
            else:
                result.append((theta, phi))
        
        return result

    @property
    def volume(self) -> float:
        """
        Volume of this polyhedron.

        Args:
            None

        Returns:
            float: The volume.

        """
        return self.convex_hull().volume

    @classmethod
    def from_sites(cls,
                   central_site: Site,
                   vertex_sites: list[Site],
                   label: str | None = None) -> CoordinationPolyhedron:
        """
        Create a :obj:`CoordinationPolyhedron` from a set of `pymatgen` :obj:`PeriodicSite` objects.

        Args:
            central_site (:obj:`pymatgen.PeriodicSite`): A `pymatgen` :obj:`PeriodicSite` object describing an atom at the nominal centre of the polyhedron.
            vertex_sites (list[:obj:`pymatgen.PeriodicSite`): A list of `pymatgen` :obj:`PeriodicSite` objects describing the atoms at the vertices.
            label (:obj:`str`, optional): An optional string used to label this coordination polyhedron.

        Returns:
            :obj:`CoordinationPolyhedron`: The :obj:`CoordinationPolyhedron` object.

        """
        vertices = [Atom(i, s) for i, s in enumerate(vertex_sites)]
        central_atom = Atom(-1, central_site)
        return cls(central_atom=central_atom, vertices=vertices, label=label)

    def vertices_by_indices(self,
                            vertex_indices: list[int]) -> list[Atom]:
        """
        Select a subset of vertices from this polyhedron with a list of vertex indices.

        Args:
            vertex_indices (list): List of vertex indices (`int`).

        Returns:
            list: A list of :obj:`Atom` objects containing the matching vertices.

        """
        return [v for v in self.vertices if v.index in vertex_indices]

    def vertex_internal_index_from_global_index(self,
                                                vertex_global_index: int) -> int:
        """
        Returns the internal index for the vertex with a given global index.

        Args:
            vertex_global_index (int): The global index for a given vertex.

        Returns:
            int: the internal index for the corresponding vertex.

        Raises:
            ValueError: If this polyhedron does not have a vertex with a global
                index equal to `vertex_global_index`.

        """
        try:
            return self.vertex_indices.index(vertex_global_index)
        except ValueError:
            raise ValueError(f"This polyhedron does not have a vertex with global index {vertex_global_index}")

    def centroid(self) -> np.ndarray:
        """
        Calculate the centroid of the vertices, accounting for periodic boundary conditions.

        Returns:
            np.ndarray: The 3D coordinates of the centroid.
        """
        return np.mean(self.minimum_image_vertex_coordinates(), axis=0)

    def centroid_to_central_atom_vector(self) -> np.ndarray:
        """
        Calculate the vector displacement from the centroid of the vertices to the central atom,
        accounting for periodic boundary conditions.

        Returns:
            np.ndarray: A 3D vector representing the displacement from the vertex centroid to the central atom.

        """
        return self.central_atom.coords - self.centroid()

    @property
    def off_centre_displacement(self) -> float:
        """
        Returns the displacment from the centroid of the polyhedron vertices to the central atom,
        accounting for periodic boundary conditions.

        Returns:
            float

        """
        return float(np.linalg.norm(self.centroid_to_central_atom_vector()))

    def radial_distortion_parameter(self, 
                                    reference: Literal['central_atom', 'centroid'] = 'centroid', 
                                    normalize: bool = True,
                                    method: Literal['MSD', 'MAD'] = 'MSD') -> float:
        """
        Calculate the radial distortion parameter for the coordination polyhedron.

        The radial distortion parameter is defined as:

        - For MSD: Delta = mean(((d_i - d_mean) / (d_mean if normalize else 1))^2)
        - For MAD: Delta = mean(abs(d_i - d_mean) / (d_mean if normalize else 1))

        where d_i are the N distances from the reference point to each vertex,
        and d_mean is their mean value.

        Args:
            reference (Literal['central_atom', 'centroid'], optional): The reference point for calculating distances. 
                Can be either 'centroid' (default) or 'central_atom'.
            normalize (bool, optional): Whether to normalize the distortion parameter
                with respect to the mean vertex distance. Default is True.
            method (Literal['MSD', 'MAD'], optional): The method to use for calculating distortion.
                Can be either 'MSD' (Mean Squared Deviation, default) or 
                'MAD' (Mean Absolute Deviation).

        Returns:
            float: The radial distortion parameter.

        Raises:
            ValueError: If an invalid reference point or method is provided.
        """
        if method not in ['MSD', 'MAD']:
            raise ValueError("Invalid method. Use 'MSD' or 'MAD'.")
        distances = self.vertex_distances(reference=reference)
        d_mean = np.mean(distances)
        if normalize:
            relative_deviations = (distances - d_mean) / d_mean
        else:
            relative_deviations = distances - d_mean
        if method == 'MAD':
            distortion = np.mean(np.abs(relative_deviations))
        else:  # MSD
            distortion = np.mean(relative_deviations**2)
        return distortion

    def vertex_angles(self, 
                      vertex_pairs: tuple[tuple[int, int], ...],
                      reference: Literal['central_atom', 'centroid'] = 'central_atom') -> np.ndarray:
        """
        Calculate angles between specified pairs of vertices.

        Args:
            vertex_pairs (tuple[tuple[int, int], ...]): Pairs of vertex global indices to calculate angles for.
            reference (Literal['central_atom', 'centroid'], optional): The reference point for vertex vectors. 
                Defaults to 'central_atom'.

        Returns:
            np.ndarray: Array of calculated angles in degrees.

        Raises:
            ValueError: If an invalid vertex index is provided.

        """
        vertex_vectors = self.vertex_vectors(reference=reference)
        angles = []

        for i1, i2 in vertex_pairs:
            try:
                v1 = vertex_vectors[self.vertex_internal_index_from_global_index(i1)]
                v2 = vertex_vectors[self.vertex_internal_index_from_global_index(i2)]
            except ValueError as e:
                raise ValueError(f"Invalid vertex index in pair ({i1}, {i2}): {str(e)}")

            # Calculate the angle between v1 and v2
            dot_product = np.dot(v1, v2)
            magnitudes = np.linalg.norm(v1) * np.linalg.norm(v2)
            cos_angle = dot_product / magnitudes
            angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))  # clip to handle floating point errors
            angles.append(np.degrees(angle))

        return np.array(angles)

def merge_coplanar_simplices(convex_hull: ConvexHull,
                             tolerance: float = 0.1) -> list[list[int]]:
    triangles_to_merge = []
    # TODO: there has to be a better way of doing this pairwise loop, e.g. using itertools.permutations
    for i, e1 in enumerate(convex_hull.equations):
        for j, e2 in enumerate(convex_hull.equations[i + 1:], i + 1):
            if np.all(e1 == e2):
                continue
            dr = e1 - e2
            distance = np.dot(dr, dr)
            if distance < tolerance:
                triangles_to_merge.append([i, j])
    merged_simplices = []
    for i, s in enumerate(convex_hull.simplices):
        if i not in np.unique(triangles_to_merge):
            merged_simplices.append(s)
    for i, j in triangles_to_merge:
        merged_simplices.append(
            np.unique([convex_hull.simplices[i], convex_hull.simplices[j]]))
    # note: this simplex is not ordered: should not be used to construct the edge_graph.
    return merged_simplices
