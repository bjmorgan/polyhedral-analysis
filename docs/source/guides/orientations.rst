Rotational orientation analysis
===============================

In many crystalline materials, coordination polyhedra can adopt discrete
orientational states with respect to the crystal lattice. Tracking which
orientation each polyhedron has -- and how orientations change in a
molecular dynamics trajectory -- gives insight into rotational dynamics
and phase transitions.


Defining reference orientations
-------------------------------

The
:class:`~polyhedral_analysis.rotation_analyser.RotationAnalyser` class
classifies polyhedra against one or more reference orientations. Each
reference orientation is defined as a set of centre-to-vertex vectors
(an Nx3 array). All reference orientations must be rotated copies of the
same geometry (they share the same point group).

For a single reference orientation, pass an Nx3 array:

.. code-block:: python

    import numpy as np
    from polyhedral_analysis.rotation_analyser import RotationAnalyser

    # A single reference orientation for a tetrahedron
    reference = np.array([[ 1.0, -1.0,  1.0],
                          [-1.0, -1.0, -1.0],
                          [ 1.0,  1.0, -1.0],
                          [-1.0,  1.0,  1.0]])

    analyser = RotationAnalyser(reference)

For multiple reference orientations (e.g. two distinct octahedral tilts
in a perovskite), pass an MxNx3 array:

.. code-block:: python

    # Two reference orientations for an octahedron
    references = np.array([
        ref_orientation_1,  # 6x3 array
        ref_orientation_2,  # 6x3 array
    ])

    analyser = RotationAnalyser(references)


Classifying polyhedron orientations
------------------------------------

The
:meth:`~polyhedral_analysis.rotation_analyser.RotationAnalyser.polyhedron_orientation`
method takes a
:class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
and returns a dictionary describing its closest discrete orientation:

.. code-block:: python

    result = analyser.polyhedron_orientation(poly)

    result['orientation_index']         # index of the closest orientation
    result['reference_geometry_index']  # which reference geometry it matches
    result['rotational_distance']       # angle from nearest orientation (radians)
    result['symmetry_measure']          # CSM for this polyhedron
    result['all_rotational_distances']  # distances to every orientation

The ``orientation_index`` identifies one of the ``M * n_G``
orientations, where M is the number of reference geometries and n_G
is the order of the point group. The
``reference_geometry_index`` maps back to the input reference array.

For working directly with vertex vector arrays (without a
``CoordinationPolyhedron`` object), use
:meth:`~polyhedral_analysis.rotation_analyser.RotationAnalyser.discrete_orientation`:

.. code-block:: python

    points = poly.vertex_vectors(reference='central_atom')
    result = analyser.discrete_orientation(points)


Orientation order parameters
----------------------------

The :mod:`polyhedral_analysis.orientation_parameters` module provides
functions for quantifying how well a polyhedron is aligned with the
Cartesian axes.

:func:`~polyhedral_analysis.orientation_parameters.oct_rotational_order_parameter`
computes an order parameter for an octahedron's alignment with the
<100> directions. It takes a 6x3 array of centroid-centred vertex
vectors:

.. code-block:: python

    from polyhedral_analysis.orientation_parameters import oct_rotational_order_parameter

    points = poly.vertex_vectors(reference='centroid')
    oct_rotational_order_parameter(points)
    # 1.0 for a perfectly aligned octahedron
    # 0.33 for a 45-degree rotated octahedron

The underlying utility functions are also available:

- :func:`~polyhedral_analysis.orientation_parameters.projection_xyz` --
  maximum projection score of a vector onto the Cartesian axes (1.0 when
  parallel to an axis, -1/3 when equally spaced between all three).
- :func:`~polyhedral_analysis.orientation_parameters.cos_theta` --
  cosine of the angle between two vectors.


Visualising orientation distributions
--------------------------------------

The
:func:`~polyhedral_analysis.plotting.plot_orientation_distribution`
function creates a contour map of the probability density for a set of
vertex orientations. This is useful for visualising the distribution of
orientations across a trajectory or set of configurations:

.. code-block:: python

    from polyhedral_analysis.plotting import plot_orientation_distribution

    # Collect (theta, phi) orientations from all polyhedra in a trajectory
    orientations = []
    for config in trajectory.configurations:
        for poly in config.polyhedra:
            orientations.extend(poly.vertex_vector_orientations())

    fig = plot_orientation_distribution(orientations)

The vertex orientations are obtained using
:meth:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.vertex_vector_orientations`
(see :doc:`geometry`).

Customisation options:

.. code-block:: python

    from cmcrameri import cm

    fig = plot_orientation_distribution(
        orientations,
        title='Vertex orientations',
        figsize=(12, 10),
        cmap=cm.batlow,  # any matplotlib-compatible colourmap
        fontsize=14,
        plot=True,        # display immediately
    )


Tracking orientations across a trajectory
------------------------------------------

Combining the
:class:`~polyhedral_analysis.rotation_analyser.RotationAnalyser` with a
:class:`~polyhedral_analysis.trajectory.Trajectory` lets you track how
individual polyhedra reorient over the course of a simulation:

.. code-block:: python

    from polyhedral_analysis.polyhedron_trajectory import PolyhedronTrajectory

    # Track the first polyhedron across all frames
    poly_traj = PolyhedronTrajectory(
        polyhedra=[config.polyhedra[0] for config in trajectory.configurations],
    )

    # Classify the orientation in each frame
    orientation_indices = [
        analyser.polyhedron_orientation(p)['orientation_index']
        for p in poly_traj.polyhedra
    ]

See :doc:`trajectories` for more on working with trajectories and
polyhedron trajectories.
