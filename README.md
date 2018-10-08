Ising_OPV
=========

This C++ software package can be used to create and analyze bulk heterojunction morphologies for further use in kinetic Monte Carlo simulation tools, such as [Excimontec](https://github.com/MikeHeiber/Excimontec). 
This package implements an Ising-based model to quickly and efficiently generate three-dimensional bulk heterojunction morphologies on a cubic lattice in a parallel computing environment. 
In addition, morphologies derived from other simulations models or experimental measurements can be imported into the tool. 
Generated or imported morphologies are then rigorously analyzed to determine important morphological features such as the domain size, tortuosity, interfacial area to volume ratio, and more. 

#### Major Features:
- Create controlled binary phase separated blends with tunable interaction energies and duration of phase separation
- Use smoothing to create idealized pure phase morphologies
- Impart controlled interfacial mixing between domains
- Create anisotropic morphologies using directionally dependent interaction energies
- Create morphology sets from imported experimental three-dimensional tomograms
- Calculate detailed structural analysis of the morphology set

<img src="./papers/v4_paper/tomogram_analysis.png" alt="Tomography data import and analysis workflow" width="650">

## Current Status

Latest release:  [![GitHub (pre-)release](https://img.shields.io/github/release/MikeHeiber/Ising_OPV/all.svg?style=flat-square)](https://github.com/MikeHeiber/Ising_OPV/releases)

All major planned features for v4.0 are now implemented and have undergone significant testing. 
However, there may still be bugs that need to be fixed. 
Please report any bugs or submit feature requests in the [Issues](https://github.com/MikeHeiber/Ising_OPV/issues) section. 

#### Continuous Integration and Testing Status:

Ising_OPV is currently being tested on [Ubuntu](https://www.ubuntu.com/) v14.04 with the [GCC compiler](https://gcc.gnu.org/) (versions 4.7, 4.8, 4.9, 5, 6, 7, and 8) and on both [Open MPI](http://www.open-mpi.org/) v1.6.5 and [MPICH](http://www.mpich.org/) v3.04 using [Travis CI](https://travis-ci.com/).

| Branch | Status |
| :------: | ------ |
| Master | [![Build Status](https://img.shields.io/travis/MikeHeiber/Ising_OPV/master.svg?style=for-the-badge)](https://travis-ci.org/MikeHeiber/Ising_OPV) |
| Development | [![Build Status](https://img.shields.io/travis/MikeHeiber/Ising_OPV/development.svg?style=for-the-badge)](https://travis-ci.org/MikeHeiber/Ising_OPV) |

Code is being tested using [googletest](https://github.com/google/googletest) with test coverage assessment by [Coveralls](https://coveralls.io/).

| Branch | Status |
| :------: | ------ |
| Master | [![Coveralls Github branch](https://img.shields.io/coveralls/github/MikeHeiber/Ising_OPV/master.svg?style=for-the-badge)](https://coveralls.io/github/MikeHeiber/Ising_OPV?branch=master) |
| Development | [![Coveralls Github branch](https://img.shields.io/coveralls/github/MikeHeiber/Ising_OPV/development.svg?style=for-the-badge)](https://coveralls.io/github/MikeHeiber/Ising_OPV?branch=development) |

## Contact

If you would like to contribute to the development of this project or would like some help in using the tool for your research, please contact me (heiber@mailaps.org) to discuss a collaboration. 
You can check out my research using this tool and other work on [Researchgate](https://www.researchgate.net/profile/Michael_Heiber).

## Using Ising_OPV

#### Building and Testing the Executable

This software package uses Message Passing Interface (MPI) to utilize parallel computing power. 
As a result, using Ising_OPV requires that an MPI library is pre-installed on your system, and the final Ising_OPV executable must be built on your specific system. 
Pre-built binaries for your system will not be supplied. 
Contact your HPC admin to determine the protocols for building MPI applications on your HPC system. 
In many cases, the HPC environment will already be configured for you, and the package comes with a default makefile that can be used with the GCC compiler. 

If you wish, you can also install MPI on your own personal workstation and then build Excimontec there as well. For development and preliminary simulation tests, sometimes it is more efficient to run on your own workstation instead of an HPC system. 
More information about common MPI packages can be found here:
- http://www.open-mpi.org/
- http://www.mpich.org/
- http://mvapich.cse.ohio-state.edu/

In addition, Ising_OPV requires a C++ 11 compatible compiler. 
The makefile that accompanies this package is setup to work with the GCC and PGI compilers, but can be easily modified for other compilers.

Once you have an MPI library installed and have an appropriate compiler, to build Ising_OPV, clone the master branch of Ising_OPV to your machine:

```git clone --recurse-submodules https://github.com/MikeHeiber/Ising_OPV```

Then, navigate to the Ising_OPV directory and run `make`. 
Once the normal build is successful, you should test Ising_OPV on your own hardware using the unit and system tests provided before you use the tool. 
Build the testing executable by running `make test`. 
Once the test build is complete, run the two test executables `./test/Ising_OPV_tests.exe` and `./test/Ising_OPV_MPI_tests.exe`.
Please report any build or testing errors in the [Issues](https://github.com/MikeHeiber/Ising_OPV/issues) section. If you do not have any build or testing errors, then you are ready to go!

#### Running Simulations

In most cases, your HPC system will use a job scheduler to manage the computing workload. 
For performing Ising_OPV simulations, it is recommended to submit batch jobs where you will request the resources needed to perform the simulation. 
An example batch script for the SLURM job scheduling system is provided with this package (slurm_script.sh). 
Similar batch scripts can also be written for TORQUE or other job schedulers.

Regardless of the job scheduler, the program execution command is essentially the same. 
Ising_OPV.exe takes one required input argument, which is the filename of the input parameter file. 
An annotated default parameter file is provided with this package (parameters_default.txt).

For example, to create 10 morphologies using 10 processors with the default parameters, the command is:

```mpiexec -n 10 Ising_OPV.exe parameters_default.txt```

Users can also import morphology sets previously generated by the Ising_OPV tool for further modification and analysis by enabling morphology import in the parameter file and running the simulation with -n set to the size of the morphology set.
This will import the morphologies (morphology_0.txt, morphology_1.txt, etc) and assign one to each processor.
The morphology files must be located in the working directory to be found and imported.

Finally, users can also import experimental 3D tomography image data, generate a morphology set from the data, and then perform further operations by enabling tomogram import in the parameter file and specifying the name of tomogram dataset. 
The tomogram metadata is imported from an XML metadata file and then that is used for interpreting a RAW binary data file that contains the image data. 
The metadata format required by Ising_OPV is defined in XML schema definition file, tomogram_metadata.xsd.
Once the tomogram data is loaded, the morphology can be segmented into a number of equally size sub-volumes to form a new morphology set, and then the rest of the analysis is performed. 
Again, the tomogram dataset files must be located in the working directory to be found and imported. 

For more detailed examples, please see the [Examples](./examples/examples.md) file.

#### Simulation Output

Ising_OPV will create several output files:
- analysis_summary.txt -- This text file will contain statistics about the set of morphologies that has been created.
- areal_composition_map_#.txt -- This text file will be created for each morphology when areal maps calculation is enabled and contains data to construct areal composition maps to check for in-plane blend compositional variations.
- areal_tortuosity_map_#.txt -- This text file will be created for each morphology when areal maps calculation is enabled and contains data to construct areal end-to-end tortuosity maps to check for in-plane tortuosity variations.
- correlation_data_#.txt -- This text file will be created for each morphology when the correlation function calculation is enabled.
- correlation_data_avg.txt -- This text file will be created when the correlation function calculation is enabled and will contain the average correlation function from all morphologies in the set.
- depth_dependent_data_#.txt -- This text file will be created for each morphology when the depth dependent calculation is enabled and will contain depth dependent blend composition and domain size data for each site type in the lattice.
- depth_dependent_data_avg.txt -- This text file will be created for each morphology when the depth dependent calculation is enabled and will contain average depth dependent blend composition and domain size data for each site type and depth dependent interfacial volume fraction data for the whole morphology set.
- interfacial_distance_histograms.txt -- This text file will be created when interfacial distance histogram calculation is enabled and will contain histogram data for each domain type averaged over all morphologies in the set.
- morphology_#_compressed.txt -- This text file will be created for each morphology generated in the set and stores the data for that morphology.
- morphology_#_cross_section.txt -- This text file will be created for each morphology generated in the set when enabled and will contain uncompressed data for a cross sectional image through the x=Length/2 plane.
- tortuosity_histograms.txt -- This text file will be created when tortuosity calculation is enabled and contains the end-to-end tortuosity histogram data for each domain type averaged over all morphologies in the set.

#### Data Analysis

For [Igor Pro](https://www.wavemetrics.com/) users, I am developing an open-source procedures package for loading, analyzing, and plotting data from Ising_OPV called [Ising_OPV_Analysis](https://github.com/MikeHeiber/Ising_OPV_Analysis). 
This is a good starting point for managing the data generated by Ising_OPV, and the Igor Pro scripting environment provides a nice playground where users can perform more advanced data analysis as needed.

#### Software API

While this tool is designed to be primarily controlled through the parameter file, API documentation for the Ising_OPV package can be viewed [here](https://mikeheiber.github.io/Ising_OPV/) to enable code developers to utilize Ising_OPV functionality as part of a larger or customized software tool. This software package is written in modern object oriented C++ and can be readily modified by other developers.

## Citing this Work

If you find Ising_OPV to be helpful for your research, please cite the original study:

[M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **2**, 014008 (2014).](http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.2.014008) [[ResearchGate]](https://www.researchgate.net/publication/264419218_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling)

If your work involves investigating the effects of morphological tortuosity, please also cite the study that introduced the tortuosity features:

[M.C. Heiber, K. Kister, A. Baumann, V. Dyakonov, C. Deibel, and T.-Q. Nguyen, Phys. Rev. Appl. **8**, 054043 (2017).](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.8.054043) [[ResearchGate]](https://www.researchgate.net/publication/321226076_Impact_of_Tortuosity_on_Charge-Carrier_Transport_in_Organic_Bulk_Heterojunction_Blends)

In addition, please also cite the DOI for the specific version that you used from [Zenodo.org](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%22592959%22&sort=-version&all_versions=True).

## Additional Reference List
Several peer-reviewed publications discuss the development and application of this software tool:

- [M. C. Heiber and A. Dhinojwala, J. Phys. Chem. C **42**, 21627 (2013).](http://pubs.acs.org/doi/abs/10.1021/jp403396v) [[ResearchGate]](https://www.researchgate.net/publication/257768674_Estimating_the_Magnitude_of_Exciton_Delocalization_in_Regioregular_P3HT)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **2**, 014008 (2014).](http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.2.014008) [[ResearchGate]](https://www.researchgate.net/publication/264419218_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling)
- [M. C. Heiber, C. Baumbach, V. Dyakonov, and C. Deibel, Phys. Rev. Lett. **114**, 136602 (2015).](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.136602) [[ResearchGate]](https://www.researchgate.net/publication/274375035_Encounter-Limited_Charge-Carrier_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)
- [M. C. Heiber, T.-Q. Nguyen, and C. Deibel, Phys. Rev. B **93**, 205204 (2016).](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.205204) [[ResearchGate]](https://www.researchgate.net/publication/302940594_Charge_Carrier_Concentration_Dependence_of_Encounter-Limited_Bimolecular_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **8**, 019902 (2017).](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.8.019902) [[ResearchGate]](https://www.researchgate.net/publication/318592832_Erratum_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling_Phys_Rev_Applied_2_014008_2014)
- [M.C. Heiber, K. Kister, A. Baumann, V. Dyakonov, C. Deibel, and T.-Q. Nguyen, Phys. Rev. Appl. **8**, 054043 (2017).](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.8.054043) [[ResearchGate]](https://www.researchgate.net/publication/321226076_Impact_of_Tortuosity_on_Charge-Carrier_Transport_in_Organic_Bulk_Heterojunction_Blends)

## Acknowledgments
Thank you to Dr. Dean M. DeLongchamp at NIST for providing access to computing resources that support the ongoing development of v4.0. 
Development of v4.0 is supported by financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).

Thank you to Klaus Kister for contributing to the development of the updated domain smoothing algorithm and the simplified morphology import procedure in v2.0 and help with testing for v3.0.

Thank you to Prof. Thuc-Quyen Nguyen at the University of California, Santa Barbara  for providing access to computing resources that supported development of v3.x. The development of v3.x used the Extreme Science and Engineering Discovery Environment (XSEDE), which is supported by National Science Foundation grant number ACI-1053575.

Thank you to Prof. Vladimir Dyakonov at the University of Würzburg and Prof. Carsten Deibel at Chemnitz University of Technology for providing access to computing resources that supported development of v2.0.

Thank you to Prof. Ali Dhinojwala and Prof. Mesfin Tsige at The University of Akron for providing access to computing resources that supported development of v1.0.
