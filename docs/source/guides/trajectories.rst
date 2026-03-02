Analysing trajectories
======================

The :class:`~polyhedral_analysis.trajectory.Trajectory` class extends
the single-structure analysis to molecular dynamics trajectories,
applying polyhedra recipes to every frame.


Loading from XDATCAR
--------------------

For VASP trajectories stored as XDATCAR files:

.. code-block:: python

    from polyhedral_analysis.trajectory import Trajectory
    from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

    recipe = PolyhedraRecipe(
        method='distance cutoff',
        coordination_cutoff=3.0,
        central_atoms='Ti',
        vertex_atoms=['O', 'F'],
    )

    trajectory = Trajectory.from_xdatcar(
        'XDATCAR',
        recipes=[recipe],
    )

If your trajectory is split across multiple XDATCAR files (e.g. from
consecutive VASP runs), use
:meth:`~polyhedral_analysis.trajectory.Trajectory.from_xdatcars`:

.. code-block:: python

    trajectory = Trajectory.from_xdatcars(
        ['XDATCAR.1', 'XDATCAR.2', 'XDATCAR.3'],
        recipes=[recipe],
    )


Loading from structures
-----------------------

You can also build a trajectory from any list of pymatgen
:class:`~pymatgen.core.structure.Structure` objects:

.. code-block:: python

    trajectory = Trajectory.from_structures(
        structures,
        recipes=[recipe],
    )

This works with structures from any source that pymatgen supports.


Progress bars and parallelism
-----------------------------

For long trajectories, enable a progress bar. The correct widget is
chosen automatically for terminals and Jupyter notebooks (via
``tqdm.auto``):

.. code-block:: python

    trajectory = Trajectory.from_xdatcar(
        'XDATCAR',
        recipes=[recipe],
        progress=True,
    )

To parallelise configuration building across multiple CPU cores:

.. code-block:: python

    trajectory = Trajectory.from_xdatcar(
        'XDATCAR',
        recipes=[recipe],
        ncores=4,
    )


Accessing configurations
------------------------

A :class:`~polyhedral_analysis.trajectory.Trajectory` contains a list
of :class:`~polyhedral_analysis.configuration.Configuration` objects,
one per frame:

.. code-block:: python

    len(trajectory)                   # number of frames
    trajectory.configurations[0]      # first frame
    trajectory.configurations[0].polyhedra  # polyhedra in first frame

    # Iterate over all frames
    for config in trajectory.configurations:
        for poly in config.polyhedra:
            print(poly.best_fit_geometry)


Polyhedron trajectories
-----------------------

A :class:`~polyhedral_analysis.polyhedron_trajectory.PolyhedronTrajectory`
tracks a single polyhedron (identified by its central atom) across all
frames of a trajectory:

.. code-block:: python

    from polyhedral_analysis.polyhedron_trajectory import PolyhedronTrajectory

    # Track the first polyhedron across all frames
    poly_traj = PolyhedronTrajectory(
        polyhedra=[config.polyhedra[0] for config in trajectory.configurations],
    )


Combining trajectories
----------------------

You can concatenate trajectories with the ``+`` operator:

.. code-block:: python

    combined = trajectory_1 + trajectory_2

Or extend an existing trajectory in place:

.. code-block:: python

    trajectory_1.extend(trajectory_2)
