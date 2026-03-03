polyhedral-analysis
===================

**polyhedral-analysis** is a Python library for analysing coordination
polyhedra in crystal structures and molecular dynamics trajectories.
Built on top of `pymatgen <https://pymatgen.org>`_, it works with any
structure that pymatgen can read.

Features include coordination polyhedron construction from distance
cutoffs or nearest-neighbour rules, continuous symmetry measures (CSM)
against reference geometries, bond lengths, angles, volumes, edge
connectivity, neighbour-sharing analysis, and trajectory tracking.

Quick example
-------------

.. code-block:: python

    from pymatgen.io.vasp import Poscar
    from polyhedral_analysis.configuration import Configuration
    from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

    structure = Poscar.from_file('POSCAR').structure
    config = Configuration(structure=structure, recipes=[recipe])

    poly = config.polyhedra[0]
    print(poly.coordination_number)   # 6
    print(poly.best_fit_geometry)     # {'geometry': 'Octahedron', ...}

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting started

   getting-started
   core-concepts

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Guides

   guides/recipes
   guides/symmetry-measures
   guides/neighbour-analysis
   guides/geometry
   guides/octahedral-analysis
   guides/orientations
   guides/trajectories

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference

   reference/api
   reference/citing
   reference/changelog
