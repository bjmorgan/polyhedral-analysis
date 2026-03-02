Core concepts
=============

This page introduces the main classes in polyhedral-analysis and how
they relate to each other.

Overview
--------

The library is built around a three-layer model:

1. A **structure** (a pymatgen :class:`~pymatgen.core.structure.Structure`)
   contains atoms with positions in a periodic unit cell.
2. A **recipe** (:class:`~polyhedral_analysis.polyhedra_recipe.PolyhedraRecipe`)
   defines how to identify coordination polyhedra in a structure.
3. A **configuration** (:class:`~polyhedral_analysis.configuration.Configuration`)
   applies one or more recipes to a structure and holds the resulting
   polyhedra.

Structures
----------

Structures are standard pymatgen
:class:`~pymatgen.core.structure.Structure` objects. You can create
them from any file format that pymatgen supports:

.. code-block:: python

    from pymatgen.io.vasp import Poscar

    structure = Poscar.from_file('POSCAR').structure

Recipes
-------

A :class:`~polyhedral_analysis.polyhedra_recipe.PolyhedraRecipe`
describes *how* to find polyhedra: which atoms are centres, which are
vertices, and what method to use for assigning coordination.

.. code-block:: python

    from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

Three coordination methods are available:

- ``'distance cutoff'`` -- include all vertex atoms within a cutoff distance.
- ``'nearest neighbours'`` -- include the *n* nearest vertex atoms.
- ``'closest centre'`` -- assign each vertex atom to its nearest centre.

See :doc:`guides/recipes` for full details on each method and how to
specify atoms.

Configurations
--------------

A :class:`~polyhedral_analysis.configuration.Configuration` applies
recipes to a structure to build a list of
:class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
objects:

.. code-block:: python

    from polyhedral_analysis.configuration import Configuration

    config = Configuration(structure=structure, recipes=[recipe])

    # All polyhedra found
    config.polyhedra

You can pass multiple recipes to find different types of polyhedra in
the same structure (e.g. octahedra and tetrahedra).

Coordination polyhedra
----------------------

A :class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
is the main object for geometric analysis. It holds a central atom and
its coordinating vertex atoms, and provides properties and methods for
inspecting the local coordination environment.

**Accessing atoms:**

.. code-block:: python

    poly = config.polyhedra[0]

    poly.central_atom          # the central Atom
    poly.vertices              # list of vertex Atoms
    poly.coordination_number   # number of vertices

**Bond distances and angles:**

.. code-block:: python

    poly.coordination_distances()  # distances from centre to each vertex
    poly.angles()                  # all vertex--centroid--vertex angles
    poly.volume                    # polyhedral volume

**Edge connectivity:**

.. code-block:: python

    poly.edge_graph  # {vertex_index: [connected_vertex_indices, ...], ...}

**Symmetry measures:**

.. code-block:: python

    poly.symmetry_measure    # {'Octahedron': 0.00, 'Trigonal prism': 16.7, ...}
    poly.best_fit_geometry   # {'geometry': 'Octahedron', 'symmetry_measure': ...}

See :doc:`guides/symmetry-measures` for more on continuous symmetry
measures.

**Neighbour analysis:**

.. code-block:: python

    poly.neighbours()                     # neighbouring polyhedra
    poly.corner_sharing_neighbour_list()  # share 1 vertex
    poly.edge_sharing_neighbour_list()    # share 2 vertices
    poly.face_sharing_neighbour_list()    # share 3+ vertices

See :doc:`guides/neighbour-analysis` for details.

What next
---------

- :doc:`guides/recipes` -- all the ways to define polyhedra recipes.
- :doc:`guides/symmetry-measures` -- continuous symmetry measures and
  best-fit geometry identification.
- :doc:`guides/neighbour-analysis` -- corner-, edge-, and face-sharing
  analysis.
- :doc:`guides/trajectories` -- analysing molecular dynamics trajectories.
