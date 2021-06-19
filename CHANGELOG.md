# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v0.2.0...HEAD


## [v0.2.0] - 2021-06-20
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v0.1.2...v0.2.0

### Added
- `Sketch.clear` method to remove all sequences currently in a `Sketch` and reset the list of minimizers.
- SIMD code to read and reverse-complement the input sequences efficiently on supported platforms (x86-64 with SSE2 or SSSE3, and ARM with NEON).
### Changed
- Split the `Sketch` type in two depending on whether the object is at the sketching stage (`Sketch`) or at the querying stage (`Mapper`).
- `Sketch.add_genome`, `Sketch.add_draft`, `Mapper.query_genome` and `Mapper.query_draft` can now be passed a Unicode string for the sequence.
### Fixed
- Integer underflow causing minimizers out of the block window to be added to the final minimizers list.


## [v0.1.2] - 2021-06-15
[v0.1.2]: https://github.com/althonos/pyrodigal/compare/v0.1.1...v0.1.2

### Changed
- Querying functions now release GIL to allow efficient parallel querying.


## [v0.1.1] - 2021-06-13
[v0.1.1]: https://github.com/althonos/pyrodigal/compare/v0.1.0...v0.1.1

### Fixed
- Source distribution missing Cython and C++ sources, thus preventing compilation.


## [v0.1.0] - 2021-06-13
[v0.1.0]: https://github.com/althonos/pyrodigal/compare/4bd3017...v0.1.0

Initial release.
