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


Working with configurations
---------------------------

Each :class:`~polyhedral_analysis.configuration.Configuration` provides
convenience attributes and methods for querying its polyhedra and atoms.

**Filtering polyhedra by label:**

When you define a recipe with a ``label`` argument (see
:doc:`recipes`), you can retrieve matching polyhedra with
:meth:`~polyhedral_analysis.configuration.Configuration.polyhedra_by_label`:

.. code-block:: python

    config.polyhedra_by_label('Ti')          # all polyhedra labelled 'Ti'
    config.polyhedra_by_label(['Ti', 'Nb'])  # multiple labels

The
:attr:`~polyhedral_analysis.configuration.Configuration.polyhedra_labels`
property lists the labels of all polyhedra:

.. code-block:: python

    config.polyhedra_labels  # ['Ti', 'Ti', 'Ti', ...]

**Accessing atoms:**

The
:attr:`~polyhedral_analysis.configuration.Configuration.central_atoms`
and
:attr:`~polyhedral_analysis.configuration.Configuration.coordination_atoms`
properties return the central and vertex atoms respectively:

.. code-block:: python

    config.central_atoms       # list of central Atom objects
    config.coordination_atoms  # list of vertex Atom objects

To look up a specific coordination atom by its global index:

.. code-block:: python

    config.coordination_atom_by_index(12)  # Atom or None


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

The ``polyhedra`` attribute gives access to the individual
:class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
objects. You can extract time-series data for any polyhedron property
(see :doc:`geometry` for all available properties):

.. code-block:: python

    # Off-centre displacement over time
    displacements = [p.off_centre_displacement for p in poly_traj.polyhedra]

    # Symmetry measure over time
    csm_values = [p.best_fit_geometry['symmetry_measure'] for p in poly_traj.polyhedra]

    # Coordination number over time
    cn_values = [p.coordination_number for p in poly_traj.polyhedra]


Combining trajectories
----------------------

You can concatenate trajectories with the ``+`` operator:

.. code-block:: python

    combined = trajectory_1 + trajectory_2

Or extend an existing trajectory in place:

.. code-block:: python

    trajectory_1.extend(trajectory_2)


Exporting to lattice_mc
------------------------

The
:meth:`~polyhedral_analysis.configuration.Configuration.to_lattice_mc`
method exports a connectivity description for use with the
`lattice_mc <https://github.com/bjmorgan/lattice_mc>`_ Monte Carlo code.
It takes a filename, a list of polyhedra labels, and a neighbour list
(typically the face-sharing neighbour list):

.. code-block:: python

    neighbour_list = config.face_sharing_neighbour_list(labels=['Ti'])
    config.to_lattice_mc('lattice.dat', labels=['Ti'], neighbour_list=neighbour_list)
