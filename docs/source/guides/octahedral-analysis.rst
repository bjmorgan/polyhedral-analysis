Octahedral analysis
===================

The :mod:`polyhedral_analysis.octahedral_analysis` module provides
specialised functions for analysing octahedral coordination
environments. These are particularly useful for mixed-anion systems
such as TiOF2, where octahedra have vertices occupied by different
species.

All functions take a
:class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
as their first argument. By default they verify that the polyhedron has
octahedral symmetry (lowest CSM) and raise ``ValueError`` if not. This
check can be disabled with ``check=False``.

.. code-block:: python

    from polyhedral_analysis.octahedral_analysis import (
        opposite_vertex_pairs,
        adjacent_vertex_pairs,
        opposite_vertex_distances,
        trans_vertex_vectors,
        trans_vector_orthogonality,
        isomer_is_trans,
        isomer_is_cis,
        isomer_is_fac,
        isomer_is_mer,
    )

    poly = config.polyhedra[0]


Opposite and adjacent vertex pairs
-----------------------------------

An octahedron has three pairs of vertices that sit opposite each other
(trans pairs) and twelve pairs of vertices connected by an edge
(adjacent or cis pairs).

:func:`~polyhedral_analysis.octahedral_analysis.opposite_vertex_pairs`
returns the three trans pairs as tuples of
:class:`~polyhedral_analysis.atom.Atom` objects:

.. code-block:: python

    opposite_vertex_pairs(poly)
    # [(Atom_3, Atom_19), (Atom_7, Atom_23), (Atom_12, Atom_16)]

:func:`~polyhedral_analysis.octahedral_analysis.adjacent_vertex_pairs`
returns the twelve cis pairs:

.. code-block:: python

    adjacent_vertex_pairs(poly)
    # [(Atom_3, Atom_7), (Atom_3, Atom_12), ...]

The distances across the octahedron between opposite vertices are
returned by
:func:`~polyhedral_analysis.octahedral_analysis.opposite_vertex_distances`:

.. code-block:: python

    opposite_vertex_distances(poly)
    # (3.96, 3.94, 3.98)

These measure the elongation or compression along each octahedral axis.


Cis/trans isomer identification
-------------------------------

For octahedra with a 2+4 vertex species split (e.g. TiO4F2), the
minority species can occupy either a trans pair (opposite each other) or
a cis pair (adjacent):

.. code-block:: python

    # Check the species split
    poly.vertex_count
    # Counter({'O': 4, 'F': 2})

    isomer_is_trans(poly)  # True if the two F are opposite
    isomer_is_cis(poly)    # True if the two F are adjacent

These raise ``ValueError`` if the vertex species split is not 2+4.

To classify all octahedra in a configuration:

.. code-block:: python

    for poly in config.polyhedra:
        if isomer_is_trans(poly):
            print(f'Polyhedron {poly.index}: trans')
        else:
            print(f'Polyhedron {poly.index}: cis')


Fac/mer isomer identification
-----------------------------

For octahedra with a 3+3 vertex species split (e.g. TiO3F3), the
two species can be arranged as fac (each species occupies a triangular
face -- all three on the same face of the octahedron) or mer (each
species forms a meridian -- one pair opposite with the third adjacent):

.. code-block:: python

    isomer_is_fac(poly)  # True for fac arrangement
    isomer_is_mer(poly)  # True for mer arrangement

These raise ``ValueError`` if the vertex species split is not 3+3.


Trans-axis vectors and orthogonality
-------------------------------------

The three vectors connecting opposite vertex pairs define the octahedral
axes.
:func:`~polyhedral_analysis.octahedral_analysis.trans_vertex_vectors`
returns these as a list of numpy arrays:

.. code-block:: python

    trans_vertex_vectors(poly)
    # [np.array([...]), np.array([...]), np.array([...])]

For a perfect octahedron, these three axes are mutually perpendicular.
:func:`~polyhedral_analysis.octahedral_analysis.trans_vector_orthogonality`
measures the deviation from this ideal: for each axis, it computes the
angle between that axis and the normal to the plane defined by the other
two. A perfect octahedron gives (0.0, 0.0, 0.0) degrees:

.. code-block:: python

    trans_vector_orthogonality(poly)
    # (0.5, 0.3, 0.8)
