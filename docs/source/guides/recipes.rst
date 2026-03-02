Defining polyhedra recipes
==========================

A :class:`~polyhedral_analysis.polyhedra_recipe.PolyhedraRecipe`
describes how to identify coordination polyhedra in a structure. This
guide covers the different ways to select atoms and the available
coordination methods.

All examples use TiOF2, which adopts the ReO3 structure with anion
disorder. We want to find the TiX6 octahedra in a 2x2x2 supercell.

.. code-block:: python

    from pymatgen.io.vasp import Poscar
    from polyhedral_analysis.configuration import Configuration
    from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

    structure = Poscar.from_file('POSCAR').structure


Selecting atoms by species
--------------------------

The simplest approach is to pass species strings. A single string
selects all atoms of that species; a list of strings selects atoms
matching any of the listed species.

.. code-block:: python

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

    config = Configuration(structure=structure, recipes=[recipe])
    config.polyhedra  # list of TiX6 polyhedra


Selecting atoms by index
-------------------------

You can pass explicit lists of integer site indices:

.. code-block:: python

    ti_indices = list(structure.indices_from_symbol('Ti'))
    of_indices = (list(structure.indices_from_symbol('O'))
                  + list(structure.indices_from_symbol('F')))

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms=ti_indices,
        vertex_atoms=of_indices,
    )

This gives the same result but is useful when you need to select a
specific subset of sites.


Selecting atoms with functions
------------------------------

For maximum flexibility, pass a callable that takes a pymatgen
:class:`~pymatgen.core.structure.Structure` and returns a sequence of
site indices:

.. code-block:: python

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms=lambda s: s.indices_from_symbol('Ti'),
        vertex_atoms=lambda s: (s.indices_from_symbol('O')
                                + s.indices_from_symbol('F')),
    )

This is particularly useful for trajectory analysis where site indices
may need to be recalculated for each frame. The helper function
:func:`~polyhedral_analysis.polyhedra_recipe.create_matching_site_generator`
creates a generator that matches sites to a reference structure,
which is useful when atom ordering may change between frames.


Distance cutoff method
----------------------

Include all vertex atoms within a specified distance of each centre:

.. code-block:: python

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

This is the most common method. The ``coordination_cutoff`` parameter
is required.


Nearest neighbours method
--------------------------

Include the *n* nearest vertex atoms to each centre:

.. code-block:: python

    recipe = PolyhedraRecipe(
        method='nearest neighbours',
        n_neighbours=6,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

The ``n_neighbours`` parameter is required. This guarantees a fixed
coordination number for every polyhedron.


Closest centre method
---------------------

Assign each vertex atom to whichever centre it is nearest to:

.. code-block:: python

    recipe = PolyhedraRecipe(
        method='closest centre',
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

Every vertex atom is assigned to exactly one polyhedron. This method
does not require a cutoff or neighbour count, but the coordination
number may vary between polyhedra.


Combining multiple recipes
--------------------------

Pass a list of recipes to
:class:`~polyhedral_analysis.configuration.Configuration` to find
different types of polyhedra in the same structure:

.. code-block:: python

    oct_recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
        label='oct',
    )
    tet_recipe = PolyhedraRecipe(
        method='nearest neighbours',
        n_neighbours=4,
        central_atoms='Si',
        vertex_atoms='O',
        label='tet',
    )

    config = Configuration(
        structure=structure,
        recipes=[oct_recipe, tet_recipe],
    )

The ``label`` parameter lets you distinguish polyhedra types when
inspecting the results.
