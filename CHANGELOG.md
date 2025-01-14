# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyfastani/compare/v0.6.1...HEAD


## [v0.6.1] - 2024-01-12
[v0.6.1]: https://github.com/althonos/pyfastani/compare/v0.6.0...v0.6.1

### Fixed
- Deployment to AUR from GitHub Actions.
- Incorrect compile actions in CMake files causing compilation failures with older CMake versions.

### Changed
- Reduce source distribution size by excluding unneeded vendored files.


## [v0.6.0] - 2024-10-22
[v0.6.0]: https://github.com/althonos/pyfastani/compare/v0.5.1...v0.6.0

### Changed
- Rewrite package build using `scikit-build-core`.
- Update documentation to use the PyData theme.
- Bump embedded `cpu-features` library to `v0.9.0`.

### Added
- Support for Python 3.13.

### Fixed
- Null pointer dereference in FASTA parser deallocation.

### Removed
- Support for Python 3.6.


## [v0.5.1] - 2024-02-20
[v0.5.1]: https://github.com/althonos/pyfastani/compare/v0.5.0...v0.5.1

### Added
- Pre-compiled wheels for Python 3.12 and PyPy 3.10.

### Fixed
- Use NEON and SSE2 unconditionally on platforms where they are available by default.


## [v0.5.0] - 2023-07-22
[v0.5.0]: https://github.com/althonos/pyfastani/compare/v0.4.1...v0.5.0

### Changed
- Bump Cython dependency to `v3.0`.

### Fixed
- Add missing signature of `Sketch.index` to type stubs.
- PyPy builds failing on missing `PyInterpreterState_GetID` function.


## [v0.4.1] - 2022-08-20
[v0.4.1]: https://github.com/althonos/pyfastani/compare/v0.4.0...v0.4.1

### Changed
- Use lists instead of arrays to store state of `Minimizers` when pickling.
- Skip serialization of minimizer index to reduce total state size.
- Make `Mapper.query_draft` and `Mapper.query_genome` return hits sorted by identity in descending order.


## [v0.4.0] - 2022-08-04
[v0.4.0]: https://github.com/althonos/pyfastani/compare/v0.3.1...v0.4.0

### Added
- Multi-threaded computation of fragment mapping in `Mapper.query_draft` and `Mapper.query_genome`.

### Fixed
- NEON-specific compile flags in `setup.py` for Aarch64 target platforms.
- Broken compilation for Arm platforms because of missing header guards.


## [v0.3.1] - 2022-07-22
[v0.3.1]: https://github.com/althonos/pyfastani/compare/v0.3.0...v0.3.1

### Added
- Slightly improve documentation in some classes.
- Sphinx documentation for the project hosted on ReadTheDocs.
- Links and instructions to install package from the Bioconda channel.


## [v0.3.0] - 2022-07-17
[v0.3.0]: https://github.com/althonos/pyfastani/compare/v0.2.1...v0.3.0

### Added
- `pickle` protocol implementation to `Mapper` and `Sketch` via `__getstate__` and `__setstate__`.
- `Minimizers` class to access the minimizers of a `Sketch` or `Mapper` object.

### Changed
- Make `Sketcher` and `Mapper` final.
- Prevent direct instantiation of `Mapper` objects.
- Update `Mapper._query_draft` to recycle memory between fragments.
- Vendor `Boost::math` headers (`v1.79`) to allow compiling without depdendencies.

### Fixed
- Broken compilation of `_fastani` extension module as `universal2` binaries on MacOS.


## [v0.2.1] - 2021-06-20
[v0.2.1]: https://github.com/althonos/pyfastani/compare/v0.2.0...v0.2.1

### Fixed
- Missing header files in the `tar.gz` distribution, preventing compilation of the wheel from source.


## [v0.2.0] - 2021-06-20
[v0.2.0]: https://github.com/althonos/pyfastani/compare/v0.1.2...v0.2.0

### Added
- `Sketch.clear` method to remove all sequences currently in a `Sketch` and reset the list of minimizers.
- SIMD code to read and reverse-complement the input sequences efficiently on supported platforms (x86-64 with SSE2 or SSSE3, and ARM with NEON).
### Changed
- Split the `Sketch` type in two depending on whether the object is at the sketching stage (`Sketch`) or at the querying stage (`Mapper`).
- `Sketch.add_genome`, `Sketch.add_draft`, `Mapper.query_genome` and `Mapper.query_draft` can now be passed a Unicode string for the sequence.
### Fixed
- Integer underflow causing minimizers out of the block window to be added to the final minimizers list.


## [v0.1.2] - 2021-06-15
[v0.1.2]: https://github.com/althonos/pyfastani/compare/v0.1.1...v0.1.2

### Changed
- Querying functions now release GIL to allow efficient parallel querying.


## [v0.1.1] - 2021-06-13
[v0.1.1]: https://github.com/althonos/pyfastani/compare/v0.1.0...v0.1.1

### Fixed
- Source distribution missing Cython and C++ sources, thus preventing compilation.


## [v0.1.0] - 2021-06-13
[v0.1.0]: https://github.com/althonos/pyfastani/compare/4bd3017...v0.1.0

Initial release.
