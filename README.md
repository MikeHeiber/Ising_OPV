Ising_OPV v2.0
=========

This software tool is used to generate and analyze model bulk heterojunction morphologies in a parallel computing environment.

### What's New in v2.0?
Version 2.0 adds several exciting new features and contains several major performance improvements.

Major New Features:
- Controlled interfacial mixing following phase separation can now be performed.
- Several advanced morphology characterization procedures that calculate interfacial distance histograms, tortuosity histograms, end-to-end paths, island volume fraction, and interfacial volume fraction.
- Users can enable periodic boundary conditions in the z-direction.

Performance Improvements:
- Output morphology file can now be created using a compressed format that dramatically reduces disk space usage.
- Domain size calculation speed has been greatly increased resulting in a significantly lower overall calculation time when creating morphologies with large domain sizes.
- Memory usage reductions increase the maximum lattice size possible.

Removed Features:
- Ability to specify a target domain size has been removed do to the challenge in performing a satisfactory estimate for the wide variety of possible input parameters.

Additional Updates:
- Users can now test the effect of asymmetric interaction energies.
- Simplified morphology import.
- More complete code commenting.
- Improved smoothing algorithm that reduces blocky artifacts when using a large rescale factor.

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

These statements can be implemented into batch scripts for running Ising_OPV v2.0 in a supercomputing environment.

### Output
Ising_OPV will create several output files:
- morphology_#_compressed.txt -- This text file will be created for each morphology generated and stores the data for that morphology.
- analysis_summary.txt -- This text file will contain statistics about the set of morphologies that has been created.  
- interfacial_distance_histograms.txt -- This text file will be created when interfacial distance histogram calculation is enabled and will contain histogram data for each domain type.
-  tortuostiy_histograms.txt -- This text file will be created when tortuosity calculation is enabled and contain the overall tortuosity histogram data for each domain type.
-  end-to-end_path_data1.txt and end-to-end_path_data2.txt -- These text files will be created when tortuosity calculation is enabled and will contain the lengths of the shortest end-to-end paths through each domain type.

### Additional Information
A simplified web-based version of this software can be found at https://nanohub.org/resources/bhjmorphology/

Several peer-reviewed publications discuss the development and application of this software tool.  See the publications below:
- [M. C. Heiber and A. Dhinojwala, J. Phys. Chem. C **42**, 21627 (2013).](http://pubs.acs.org/doi/abs/10.1021/jp403396v)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **2**, 014008 (2014).](http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.2.014008)
- [M. C. Heiber, C. Baumbach, V. Dyakonov, and C. Deibel, Phys. Rev. Lett. **114**, 136602 (2015).](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.136602)

If you would like some assitance in customizing this software tool for your particular interest or application, please contact me to discuss collaboration options.

