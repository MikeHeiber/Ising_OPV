# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

--------------------------------------------------------------------------------------------------------------------------------

## [Unreleased]

## [v4.0.2] - 2019-04-17 - Tortuosity Histogram Bugfix

### Added
- Utils (calculateProbabilityHist) - Added new overloaded function definition that allows users to specify the bin positions

### Changed
- .travis.yml - Extended copyright statement to 2019
- Doxyfile - Updated project version to v4.0.2
- docs - Updated documentation using updated Doxygen 1.8.15
- LICENSE - Extended copyright statement to 2019
- Lattice - Extended copyright statements to 2019
- main.cpp - Extended copyright statement to 2019
- main.cpp - Reduced the tortuosity histogram bin size from 0.02 to 0.01 to produce higher resolution data
- Morphology - Extended copyright statements to 2019
- Parameters - Extended copyright statements to 2019
- Utils - Extended copyright statements to 2019
- Utils (calculateProbabilityHist) - Refactored functions to use data vector size for normalization instead of unnecessarily counting 
- Version - Extended copyright statements to 2019
- Version - Updated Current_version namespace variable to v4.0.02
- test/test.cpp - Extended copyright statement to 2019
- test/test_mpi.cpp - Extended copyright statement to 2019

### Fixed
- main.cpp - Tortuosity histogram output so that the data from the first bin is now correctly output to the data file
- README.md - Current status section to show latest stable release badge

## [v4.0.1] - 2018-12-12 - Tomogram Import Bugfix

### Added
- README.md - Link to JOSS paper
- Morphology (executeIsingSwapping) - Cast return of distance function to an int to prevent compiler warning
- Morphology (importTomogramMorphologyFile) - Cast return of ifstream's tellg function to an int to prevent compiler warning
- Utils (vector_which_median) - Cast return of distance function to int to prevent compiler warning

### Changed
- Doxyfile - Updated file input and output paths to relative paths
- Version - Updated Current_version namespace variable to v4.0.1
- Updated Doxygen documentation
- Morphology (calculatePathDistances) - Moved declaration of local z variable near it's first use to prevent potential conflict with other z variable used in a subsequent loop
- main.cpp (main) - Imported morphology filenames now use main's reusable filename string variable

### Fixed
- Morphology (importTomogramMorphologyFile) - Corrected error loading xml metadata files on Windows by telling the fopen function to open the file as a binary file
- Morphology (importTomogramMorphologyFile) - Corrected miscalculation of the extracted sublattice dimensions that was causing some tomogram files not to be imported
- Parameters (importParameters) - Corrected problem importing the Mixed_frac parameter from the parameter file as an integer instead of a double

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
