Neighbour analysis
==================

Coordination polyhedra that share vertex atoms are neighbours. The type
of sharing -- corner, edge, or face -- is determined by the number of
shared vertices.

All examples assume you have a
:class:`~polyhedral_analysis.configuration.Configuration` with
polyhedra already constructed:

.. code-block:: python

    poly = config.polyhedra[0]


Finding neighbours
------------------

The :meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.neighbours`
method returns all neighbouring polyhedra (those sharing at least one
vertex):

.. code-block:: python

    poly.neighbours()  # tuple of CoordinationPolyhedron objects


Shared vertices
---------------

To see which vertices are shared with each neighbour, use
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.neighbours_by_index_and_shared_vertices`:

.. code-block:: python

    poly.neighbours_by_index_and_shared_vertices()
    # {3: (8, 12), 5: (8,), 7: (12, 14, 16), ...}

The keys are neighbouring polyhedron indices and the values are tuples
of shared vertex atom indices.


Corner-sharing neighbours
-------------------------

Polyhedra that share exactly one vertex:

.. code-block:: python

    poly.corner_sharing_neighbour_list()
    # (5, 9, ...)

Returns a tuple of indices for corner-sharing neighbours.


Edge-sharing neighbours
-----------------------

Polyhedra that share exactly two vertices:

.. code-block:: python

    poly.edge_sharing_neighbour_list()
    # (3, ...)

Returns a tuple of indices for edge-sharing neighbours.


Face-sharing neighbours
-----------------------

Polyhedra that share three or more vertices:

.. code-block:: python

    poly.face_sharing_neighbour_list()
    # (7, ...)

Returns a tuple of indices for face-sharing neighbours.


Configuration-level neighbour lists
------------------------------------

The :class:`~polyhedral_analysis.configuration.Configuration` class
provides a convenience method to compute face-sharing neighbours for
all polyhedra at once:

.. code-block:: python

    config.face_sharing_neighbour_list(labels=['Ti'])
    # {0: (7,), 1: (), 3: (7,), ...}

This returns a dictionary mapping each polyhedron index to a tuple of
its face-sharing neighbour indices.
