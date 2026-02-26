# polyhedral-analysis

![Build Status](https://github.com/bjmorgan/polyhedral-analysis/actions/workflows/build.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/bjmorgan/polyhedral-analysis/badge.svg?branch=main)](https://coveralls.io/github/bjmorgan/polyhedral-analysis?branch=main)
[![Documentation Status](https://readthedocs.org/projects/polyhedral-analysis/badge/?version=latest)](http://polyhedral-analysis.readthedocs.io/en/latest/?badge=latest)

`polyhedral-analysis` is a Python module for analysing coordination polyhedra in crystal structures and molecular dynamics trajectories. It identifies coordination environments around central atoms and computes geometric properties such as distortion measures, bond lengths, angles, volumes, and connectivity.

Built on top of [pymatgen](https://pymatgen.org), the package works with any structure that pymatgen can read.

## Features

- Coordination polyhedron construction from distance cutoffs, nearest neighbours, or closest-centre assignment
- Continuous symmetry measures (CSM) against reference geometries (tetrahedron, octahedron, cube, etc.)
- Best-fit geometry identification across all reference polyhedra for a given coordination number
- Bond lengths, angles, volumes, and edge connectivity graphs
- Corner-, edge-, and face-sharing neighbour analysis
- Off-centre displacement and radial distortion parameters
- Vertex vector orientation analysis with stereographic projection plotting
- Trajectory analysis for tracking polyhedral distortions over molecular dynamics runs

## Installation

```bash
pip install polyhedral-analysis
```

### Requirements

- Python 3.11+
- numpy
- pymatgen >= 2024.7.18
- scipy
- bsym

## Quick start

```python
from pymatgen.io.vasp import Poscar
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

# Define a recipe for octahedral coordination
recipe = PolyhedraRecipe(
    method='distance cutoff',
    coordination_cutoff=3.0,
    central_atoms='Ti',
    vertex_atoms=['O', 'F'],
)

# Load a structure and build polyhedra
structure = Poscar.from_file('POSCAR').structure
config = Configuration(structure=structure, recipes=[recipe])

# Inspect a polyhedron
poly = config.polyhedra[0]
print(poly.coordination_number)     # 6
print(poly.best_fit_geometry)       # {'geometry': 'Octahedron', 'symmetry_measure': ...}
print(poly.volume)                  # polyhedral volume
print(poly.coordination_distances())  # list of bond lengths
```

For more examples, see the [example notebooks](./docs/examples).

## Documentation

API documentation is available at [polyhedral-analysis.readthedocs.io](https://polyhedral-analysis.readthedocs.io/en/latest/).
