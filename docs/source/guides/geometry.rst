Geometric analysis
==================

The :doc:`/core-concepts` page covers the basics of working with
coordination polyhedra: bond distances, angles, volume, edge graphs,
symmetry measures, and neighbour sharing. This guide covers the
remaining geometric analysis tools available on
:class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`.

All examples assume you have a
:class:`~polyhedral_analysis.configuration.Configuration` with
polyhedra already constructed:

.. code-block:: python

    poly = config.polyhedra[0]


Vertex species composition
--------------------------

In mixed-anion systems (e.g. TiOF2, where Ti is coordinated by both O
and F), it is useful to know which species occupy the vertex sites.

The
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_labels`
property returns the species label for each vertex:

.. code-block:: python

    poly.vertex_labels
    # ['F', 'F', 'O', 'O', 'O', 'O']

To get a count of vertices grouped by species, use
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_count`:

.. code-block:: python

    poly.vertex_count
    # Counter({'O': 4, 'F': 2})

To pair each vertex distance with its species label, use
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_distances_and_labels`:

.. code-block:: python

    poly.vertex_distances_and_labels()
    # ((1.98, 'F'), (2.01, 'F'), (1.95, 'O'), (1.96, 'O'), (1.97, 'O'), (1.97, 'O'))

This accepts an optional ``reference`` argument (``'central_atom'`` or
``'centroid'``) to choose the reference point for the distance
calculation.


Distances and vectors
---------------------

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_distances`
method returns an array of distances from either the central atom or
the centroid to each vertex:

.. code-block:: python

    poly.vertex_distances()                       # from central atom (default)
    poly.vertex_distances(reference='centroid')    # from centroid

The existing
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.coordination_distances`
method is a convenience wrapper that always measures from the central
atom and returns a list rather than a numpy array.

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vectors`
method returns the vectors from a reference point to each vertex:

.. code-block:: python

    poly.vertex_vectors()                         # from centroid (default)
    poly.vertex_vectors(reference='central_atom') # from central atom

These vectors are the basis for most other geometric calculations
(angles, orientations, symmetry measures).

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.centroid`
method returns the centroid of the vertex positions, correctly
accounting for periodic boundary conditions:

.. code-block:: python

    poly.centroid()  # np.array([x, y, z])


Off-centre displacement
-----------------------

In ferroelectric and related materials, the displacement of the central
atom away from the geometric centre of its coordination shell is
physically significant.

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.centroid_to_central_atom_vector`
method returns the 3D displacement vector:

.. code-block:: python

    poly.centroid_to_central_atom_vector()  # np.array([dx, dy, dz])

The scalar magnitude is available as the
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.off_centre_displacement`
property:

.. code-block:: python

    poly.off_centre_displacement  # 0.042


Radial distortion
-----------------

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.radial_distortion_parameter`
method quantifies how much the vertex distances deviate from their
mean value. By default it computes the normalised mean squared
deviation (MSD):

.. code-block:: python

    poly.radial_distortion_parameter()

The ``method`` argument selects between MSD and mean absolute deviation
(MAD):

.. code-block:: python

    poly.radial_distortion_parameter(method='MAD')

Other options control the reference point and normalisation:

.. code-block:: python

    # Measure from the central atom instead of the centroid
    poly.radial_distortion_parameter(reference='central_atom')

    # Un-normalised (absolute deviations, not divided by the mean distance)
    poly.radial_distortion_parameter(normalize=False)


Angular analysis
----------------

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_angles`
method calculates angles between specific pairs of vertices, identified
by their global atom indices:

.. code-block:: python

    pairs = ((3, 7), (3, 12))
    poly.vertex_angles(pairs)  # np.array([90.1, 88.5])

This complements the existing
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.angles`
method, which returns *all* pairwise vertex--centroid--vertex angles.

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vector_orientations`
method returns the angular orientation of each vertex vector as
``(theta, phi)`` pairs. Theta is the angle from [0, 0, 1] (0--180
degrees) and phi is the angle from [1, 0, 0] (-180 to +180 degrees):

.. code-block:: python

    poly.vertex_vector_orientations()  # from centroid (default)
    # [(45.0, 90.0), (135.0, -90.0), ...]

    # In radians
    poly.vertex_vector_orientations(units='radians')

    # Also return the distance for each vertex
    poly.vertex_vector_orientations(return_distance=True)
    # [(45.0, 90.0, 1.98), (135.0, -90.0, 2.01), ...]


Vector projections
------------------

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vector_projections`
method projects each vertex vector onto one or more input vectors. This
is useful for measuring the component of each bond along a specific
crystallographic direction:

.. code-block:: python

    import numpy as np

    # Project onto the c-axis
    c_axis = np.array([0.0, 0.0, 1.0])
    poly.vertex_vector_projections(c_axis)
    # np.array of shape (n_vertices, 1)

    # Project onto multiple directions at once
    axes = np.array([[1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0]])
    poly.vertex_vector_projections(axes)
    # np.array of shape (n_vertices, 3)


Topology: faces, edges, and convex hull
----------------------------------------

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.faces`
method returns the faces of the polyhedron as tuples of vertex atom
indices:

.. code-block:: python

    poly.faces()
    # ((3, 7, 12), (3, 7, 16), (3, 12, 19), ...)

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.edge_vertex_indices`
method returns all edges as sorted pairs:

.. code-block:: python

    poly.edge_vertex_indices()
    # ((3, 7), (3, 12), (3, 16), (7, 12), ...)

This is an alternative representation of the same connectivity as the
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.edge_graph`
property (documented in :doc:`/core-concepts`), which uses an adjacency
list format instead.

For custom analysis, the underlying
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.convex_hull`
method returns a scipy
:class:`~scipy.spatial.ConvexHull` object:

.. code-block:: python

    hull = poly.convex_hull()


Identity and utility
--------------------

**Labels and indices:**

Each polyhedron has a
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.label`
(defaults to the central atom species) and a method to change it:

.. code-block:: python

    poly.label              # 'Ti'
    poly.set_label('Ti_1')

The
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.index`
property returns the global atom index of the central atom, and
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_indices`
returns the global indices of the vertex atoms:

.. code-block:: python

    poly.index           # 0
    poly.vertex_indices  # [3, 7, 12, 16, 19, 23]

The Cartesian coordinates of the vertex atoms are available as an Nx3
array:

.. code-block:: python

    poly.vertex_coords  # np.array of shape (6, 3)

**Periodic boundary conditions:**

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.minimum_image_vertex_coordinates`
method returns vertex coordinates where each vertex is the closest
periodic image to the central atom:

.. code-block:: python

    poly.minimum_image_vertex_coordinates()  # np.array of shape (6, 3)

**Shared vertices:**

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.intersection`
method returns the atom indices shared with another polyhedron:

.. code-block:: python

    poly_a = config.polyhedra[0]
    poly_b = config.polyhedra[1]
    poly_a.intersection(poly_b)  # (12, 16)

This is the underlying operation used by the neighbour analysis methods
(see :doc:`neighbour-analysis`).

**Comparing polyhedra:**

Three methods provide different levels of equality testing:

.. code-block:: python

    # Same vertex atom indices?
    poly_a.equal_vertices(poly_b)

    # Same edge connectivity graph?
    poly_a.equal_edge_graph(poly_b)

    # Same central atom AND same vertices?
    poly_a.equal_members(poly_b)

The default ``==`` operator uses ``equal_edge_graph``.

**Index lookups:**

To retrieve a subset of vertex atoms by their global indices:

.. code-block:: python

    poly.vertices_by_indices([12, 16])  # [Atom, Atom]

To convert a global atom index to the internal (0-based) vertex index:

.. code-block:: python

    poly.vertex_internal_index_from_global_index(12)  # 2

**Creating polyhedra directly:**

The
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.from_sites`
class method creates a polyhedron directly from pymatgen
:class:`~pymatgen.core.sites.Site` objects, without going through a
:class:`~polyhedral_analysis.configuration.Configuration`:

.. code-block:: python

    from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron

    poly = CoordinationPolyhedron.from_sites(
        central_site=central_site,
        vertex_sites=vertex_sites,
        label='Ti',  # optional
    )
