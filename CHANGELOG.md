# Changelog

## 0.2.0

### Performance

- Reduced CSM (continuous symmetry measure) permutations using symmetry equivalence via bsym. Only symmetry-inequivalent vertex permutations are evaluated, giving 12-49x speedups depending on geometry symmetry.

### Bug fixes

- Removed duplicate "Dodecahedron with triangular faces" entry from CN=8 geometry list.

### Breaking changes

- Now requires Python >=3.11 (previously >=3.7).
- Removed numpy<2.0 version cap; now compatible with numpy 2.x.
- New runtime dependency: bsym.

### Other changes

- Added off-centre displacement calculation.
- Revised vertex distance and vertex pair angle calculations.
- Added octahedral orthogonality metric.
- Added radial distortion parameter.
- Added polyhedron centroid methods.
- Added mapping from global to local vertex indices.
- CI updated to test Python 3.11, 3.12, 3.13, 3.14.

## 0.1

Initial release.
