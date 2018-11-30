# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

--------------------------------------------------------------------------------------------------------------------------------

## [Unreleased]

### Added
- README.md - Link to JOSS paper

## [v4.0.0] - 2018-11-29 - Final Release

This release only changes documentation files and does not contain modifications to the source code.

### Added
- CONTRIBUTING.md - New file with instructions for how others can contribute to the project
- README.md - Link to new CONTRIBUTING file
- CHANGELOG.md - New file detailing the changes for each release
- README.md - Link to new CHANGELOG file
- papers/v4_paper/paper.bib - Several missing reference entries
- papers/v4_paper/paper.md - Acknowledgement of Dr. Andrew A. Herzing for TEM tomo data

### Changed
- README.md - Current release status to note that v4.0 is now stable
- papers/v4_paper/paper.md - Image alt text to be more descriptive
- papers/v4_paper/paper.md - Formating of plus signs using \texttt command
- Version.cpp - Current_version namespace variable to v4.0.0
- Doxyfile - updated the version number to v4.0.0
- docs - generated updated Doxygen documentation

### Fixed
- README.md - Instructions for running the MPI tests
- papers/v4_paper/paper.md - Problem with underscore by adding a backslash
