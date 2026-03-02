Continuous symmetry measures
============================

A continuous symmetry measure (CSM) quantifies how much a coordination
polyhedron deviates from an ideal reference geometry. A CSM of zero
indicates perfect agreement; larger values indicate greater distortion.

The CSM is computed using the approach of Pinsky and Avnir [1]_, which
finds the optimal rotation, scaling, and vertex permutation to minimise
the distance between the observed and reference vertex positions via
SVD-based Procrustes analysis.


Computing symmetry measures
---------------------------

Every :class:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron`
has a :attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.symmetry_measure`
property that returns the CSM against all reference geometries defined
for its coordination number:

.. code-block:: python

    poly = config.polyhedra[0]
    poly.symmetry_measure
    # {'Octahedron': 1.2e-30, 'Trigonal prism': 16.7, ...}

The keys are geometry names and the values are CSM values.


Best-fit geometry
-----------------

The
:attr:`~polyhedral_analysis.coordination_polyhedron.CoordinationPolyhedron.best_fit_geometry`
property returns the reference geometry with the smallest CSM:

.. code-block:: python

    poly.best_fit_geometry
    # {'geometry': 'Octahedron', 'symmetry_measure': 1.2e-30}


Available reference geometries
------------------------------

The following reference geometries are available, grouped by
coordination number:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - CN
     - Geometries
   * - 4
     - Tetrahedron
   * - 5
     - Trigonal bipyramid, Square pyramid
   * - 6
     - Octahedron, Trigonal prism
   * - 7
     - Pentagonal bipyramid, Square-face capped trigonal prism,
       Face-capped octahedron
   * - 8
     - Cube, Square antiprism, Square-face bicapped trigonal prism,
       Triangular-face bicapped trigonal prism,
       Dodecahedron with triangular faces, Hexagonal bipyramid,
       Bicapped octahedron (opposed cap faces),
       Bicapped octahedron (cap faces with one atom in common),
       Bicapped octahedron (cap faces with one edge in common)

Reference coordinates are defined in
:mod:`polyhedral_analysis.reference_geometries`.


Using SymmetryMeasure directly
------------------------------

For lower-level control, you can construct a
:class:`~polyhedral_analysis.symmetry_measure.SymmetryMeasure` and
compute the CSM against custom vertex coordinates:

.. code-block:: python

    from polyhedral_analysis.symmetry_measure import SymmetryMeasure

    sm = SymmetryMeasure.from_name('Octahedron')
    csm_value = sm.minimum_symmetry_measure(vertex_vectors)

The :meth:`~polyhedral_analysis.symmetry_measure.SymmetryMeasure.minimum_symmetry_measure`
method searches over all symmetry-inequivalent vertex permutations,
using the point-group symmetry of the reference geometry to reduce the
number of permutations that must be tested.


References
----------

.. [1] M. Pinsky and D. Avnir, "Continuous Symmetry Measures. 5. The
   Classical Polyhedra", *Inorganic Chemistry*, **37**, 5575-5582 (1998).
