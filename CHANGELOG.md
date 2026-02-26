# Changelog

## 0.3.0

### Other changes

- Migrated from setup.py to pyproject.toml.
- Removed requirements.txt, .coveragerc, and .mypy.ini (consolidated into pyproject.toml).
- Version now read at runtime via importlib.metadata.
- Updated type hints to modern Python 3.11+ syntax.
- Updated README with description, features, installation, and quick start example.

## 0.2.0

### Performance

- Reduced CSM (continuous symmetry measure) permutations using symmetry equivalence via bsym. Only symmetry-inequivalent vertex permutations are evaluated, giving 12-49x speedups depending on geometry symmetry.

### Bug fixes

- Removed duplicate "Dodecahedron with triangular faces" entry from CN=8 geometry list.
- Fixed mypy errors across codebase.

### Breaking changes

- Now requires Python >=3.11 (previously >=3.7).
- Removed numpy<2.0 version cap; now compatible with numpy 2.x.
- Minimum pymatgen version bumped to >=2024.7.18.
- New runtime dependencies: bsym, cmcrameri.

### New features

- Added orientation plotting.
- Added off-centre displacement calculation.
- Added octahedral orthogonality metric.
- Added radial distortion parameter.
- Added polyhedron centroid methods.
- Added mapping from global to local vertex indices.
- Revised vertex distance and vertex pair angle calculations.

### Other changes

- CI updated to test Python 3.11, 3.12, 3.13, 3.14.
- CI now triggers only on pull requests.

## 0.1

Initial release.
