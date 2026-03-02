Getting started
===============

Installation
------------

Install from PyPI:

.. code-block:: bash

    pip install polyhedral-analysis

Requirements
~~~~~~~~~~~~

- Python 3.11+
- numpy
- pymatgen >= 2024.7.18
- scipy
- bsym

All dependencies are installed automatically by pip.

Quick start
-----------

The core workflow has three steps: load a structure, define a recipe for
finding polyhedra, and build a configuration.

.. code-block:: python

    from pymatgen.io.vasp import Poscar
    from polyhedral_analysis.configuration import Configuration
    from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

    # 1. Load a structure
    structure = Poscar.from_file('POSCAR').structure

    # 2. Define a recipe
    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

    # 3. Build the configuration
    config = Configuration(structure=structure, recipes=[recipe])

    # Inspect a polyhedron
    poly = config.polyhedra[0]
    print(poly.coordination_number)       # 6
    print(poly.best_fit_geometry)         # {'geometry': 'Octahedron', ...}
    print(poly.coordination_distances())  # list of bond lengths

Next steps
----------

- :doc:`core-concepts` explains the key classes and how they fit together.
- :doc:`guides/recipes` covers all the ways to define polyhedra recipes.
- :doc:`guides/symmetry-measures` introduces continuous symmetry measures.
- :doc:`guides/trajectories` shows how to analyse molecular dynamics trajectories.
