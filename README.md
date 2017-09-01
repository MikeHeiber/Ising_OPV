Ising_OPV
=========

This highly optimized software tool uses an Ising-based model to quickly and efficiently generate three-dimensional bulk heterojunction morphologies on a cubic lattice in a parallel computing environment. Generated or imported morphologies are then rigorously analyzed to determine important morphological features such as the domain size, tortuosity, interfacial area to volume ratio, and more.  Development of this tool represents an attempt to standardize the Ising-based morphology model for the reliable creation and analysis of morphologies to be further used in kinetic Monte Carlo simulations of organic photovoltaic devices. If you would like some assistance in customizing this software tool for your particular research interest or application, please contact me to discuss collaboration options or feel free to contribute to the development of this open-source software tool.

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

These statements can be implemented into batch scripts for running Ising_OPV in a supercomputing environment.

### Output
Ising_OPV will create several output files:
- morphology_#_compressed.txt -- This text file will be created for each morphology generated and stores the data for that morphology.
- analysis_summary.txt -- This text file will contain statistics about the set of morphologies that has been created.
- correlation_data_avg.txt -- This text file will contain the average correlation function data for the morphology set.
- interfacial_distance_histograms.txt -- This text file will be created when interfacial distance histogram calculation is enabled and will contain histogram data for each domain type.
-  tortuostiy_histograms.txt -- This text file will be created when tortuosity calculation is enabled and contain the overall tortuosity histogram data for each domain type.
-  end-to-end_path_data1.txt and end-to-end_path_data2.txt -- These text files will be created when tortuosity calculation is enabled and will contain the lengths of the shortest end-to-end paths through each domain type.

### Additional Information
Several peer-reviewed publications discuss the development and application of this software tool:
- [M. C. Heiber and A. Dhinojwala, J. Phys. Chem. C **42**, 21627 (2013).](http://pubs.acs.org/doi/abs/10.1021/jp403396v) [[ResearchGate]](https://www.researchgate.net/publication/257768674_Estimating_the_Magnitude_of_Exciton_Delocalization_in_Regioregular_P3HT)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **2**, 014008 (2014).](http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.2.014008) [[ResearchGate]](https://www.researchgate.net/publication/264419218_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling)
- [M. C. Heiber, C. Baumbach, V. Dyakonov, and C. Deibel, Phys. Rev. Lett. **114**, 136602 (2015).](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.136602) [[ResearchGate]](https://www.researchgate.net/publication/274375035_Encounter-Limited_Charge-Carrier_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)
- [M. C. Heiber, T.-Q. Nguyen, and C. Deibel, Phys. Rev. B **93**, 205204 (2016).](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.205204) [[ResearchGate]](https://www.researchgate.net/publication/302940594_Charge_Carrier_Concentration_Dependence_of_Encounter-Limited_Bimolecular_Recombination_in_Phase-Separated_Organic_Semiconductor_Blends)
- [M. C. Heiber and A. Dhinojwala, Phys. Rev. Appl. **8**, 019902 (2017).](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.8.019902) [[ResearchGate]](https://www.researchgate.net/publication/318592832_Erratum_Efficient_Generation_of_Model_Bulk_Heterojunction_Morphologies_for_Organic_Photovoltaic_Device_Modeling_Phys_Rev_Applied_2_014008_2014)

### Acknowledgements
Thank you to Dr. Dean DeLongchamp at NIST for providing access to computing resources that support the ongoing development of v4.0. Development of v4.0 is supported by financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).

Thank you to Klaus Kister for contributing to the development of the updated domain smoothing algorithm and the simplified morphology import procedure in v2.0 and help with testing for v3.0.

Thank you to Prof. Ali Dhinojwala and Prof. Mesfin Tsige at The University of Akron for providing access to computing resources that supported development of v1.0.

Thank you to Prof. Vladimir Dyakonov at the University of Würzburg and Prof. Carsten Deibel at Chemitz University of Technology for providing access to computing resources that supported development of v2.0.

Thank you to Prof. Thuc-Quyen Nguyen at the University of California, Santa Barbara  for providing access to computing resources that supported development of v3.x. The development of v3.x used the Extreme Science and Engineering Discovery Environment (XSEDE), which is supported by National Science Foundation grant number ACI-1053575.



