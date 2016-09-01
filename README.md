Ising_OPV v3.0
=========

This software tool is used to generate and analyze model bulk heterojunction morphologies in a parallel computing environment.

### What's New in v3.0?
Version 3.0 fixes a significant bug present in the previous versions, adds several exciting new features, and contains several major performance improvements.

Major New Features:
- Added the ability to shrink the lattice by a specified integer factor
- Added the ability to use a 3D checkerboard starting configuration instead of a random blend
- Added the ability to modify the interaction energy in one of the directions to allow anisotropic domain growth
- Improved the analysis file output to include a list of the properties of all morphologies in the set
- Added output of the average correlation function for the morphology set to the correlation_data_avg.txt file

Minor New Features:
- Added the ability to extend the correlation function calculation out to the second correlation maximum
- Added characterization of domain anisotropy as a standard metric
- Improved the analysis file to include a header that specifies which version of the software was used to generate the morphology set
- Improved the analysis file output to also specify which morphology has the median domain size and the median tortuosity

Performance Improvements:
- Increased the speed of the phase separation process by increasing the speed of the energy calculation during the Ising swapping stage
- Improved the speed of the tortuosity calculation by increasing the speed of the pathfinding and path distance calculations
- Added the ability to enable a reduced memory usage tortuosity calculation algorithm.  This algorithm is significantly slower, but it is useful when simulating very large lattices where memory limits are reached.

Bug Fixes:
- Corrected a major bug with the random number generator used to choose which neighboring site to use for a swapping attempt (This bug caused anisotropic domain growth in the previous versions of Ising_OPV)
- Improved the correlation function to catch rare cases where the correlation function does not cross the mix fraction value, and in these cases the domain size is set to the position of the first correlation minimum

### Compiling
Compiling requires the boost library for random number generation and an MPI library for parallel processing.

More information about these packages can be found here:
- http://www.boost.org/
- http://www.mpich.org/ or http://www.open-mpi.org/

### Usage
Ising_OPV.exe takes one required input argument, which is the filename of the input parameter file.

An example parameter file is provided with parameters_default.txt

As an example, to create one morphology on a single processor, the command is:
>    Ising_OPV.exe parameters_default.txt

To run in a parallel processing environment and create 10 morphologies on 10 processors, the command is:
>    mpiexec -n 10 Ising_OPV.exe parameters_default.txt

An optional second input argument is the path to an existing morphology file for importing a previously created morphology into the program for further modifications.

As an example, to import a compressed morphology that is in current working directory on a single processor, the command is:
>    Ising_OPV.exe parameters_default.txt ./morphology_0_compressed.txt

To import an entire set of morphologies in a parallel processing environment:
>    mpiexec -n 10 Ising_OPV.exe parameters_default.txt ./morphology_#_compressed.txt

This command will import 10 morphologies (morphology_0_compressed.txt, morphology_1_compressed.txt, etc) and assign one to each processor.

These statements can be implemented into batch scripts for running Ising_OPV v3.0 in a supercomputing environment.

### Output
Ising_OPV will create several output files:
- morphology_#_compressed.txt -- This text file will be created for each morphology generated and stores the data for that morphology.
- analysis_summary.txt -- This text file will contain statistics about the set of morphologies that has been created.
- correlation_data_avg.txt -- This text file will contain the average correlation function data for the morphology set.
- interfacial_distance_histograms.txt -- This text file will be created when interfacial distance histogram calculation is enabled and will contain histogram data for each domain type.
-  tortuostiy_histograms.txt -- This text file will be created when tortuosity calculation is enabled and contain the overall tortuosity histogram data for each domain type.
-  end-to-end_path_data1.txt and end-to-end_path_data2.txt -- These text files will be created when tortuosity calculation is enabled and will contain the lengths of the shortest end-to-end paths through each domain type.

### Additional Information
A simplified web-based version of this software (based on v1.0) can be found at https://nanohub.org/resources/bhjmorphology/

Several peer-reviewed publications discuss the development and application of this software tool.  See the publications below:
- [M. C. Heiber and A. Dhinojwala, J. Phys. Chem. C **42**, 21627 (2013).](http://pubs.acs.org/doi/abs/10.1021/jp403396v) [[ResearchGate]](https://www.researchgate.net/publication/257768674_Estimating_the_Magnitude_of_Exciton_Delocalization_in_Regioregular_P3HT)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **2**, 014008 (2014).](http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.2.014008) [[ResearchGate]](https://www.researchgate.net/publication/264419218_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling)
- [M. C. Heiber, C. Baumbach, V. Dyakonov, and C. Deibel, Phys. Rev. Lett. **114**, 136602 (2015).](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.136602) [[ResearchGate]](https://www.researchgate.net/publication/274375035_Encounter-Limited_Charge-Carrier_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)
- [M. C. Heiber, T.-Q. Nguyen, and C. Deibel, Phys. Rev. B **93**, 205204 (2016).](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.205204) [[ResearchGate]](https://www.researchgate.net/publication/302940594_Charge_Carrier_Concentration_Dependence_of_Encounter-Limited_Bimolecular_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)

If you would like some assistance in customizing this software tool for your particular interest or application, please contact me to discuss collaboration options.

### Acknowledgements
Thank you to Klaus Kister for contributing to the development of the updated domain smoothing algorithm and the simplified morphology import procedure in v2.0 and help with testing in v3.0.

Thank you to Prof. Ali Dhinojwala and Prof. Mesfin Tsige at The University of Akron, Prof. Vladimir Dyakonov at the University of Würzburg, and Prof. Carsten Deibel at Chemitz University of Technology, and Prof. Thuc-Quyen Nguyen at the University of California, Santa Barbara for providing access to cluster computing resources that facilitated the testing of this software tool during various stages of development.  This work also used the Extreme Science and Engineering Discovery Environment (XSEDE), which is supported by National Science Foundation grant number ACI-1053575.

