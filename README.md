Ising_OPV
=========

This software tool is used to generate model bulk heterojunction morphologies in a parallel computing environment for organic solar cell simulations.

Compiling:
Compiling requires the boost library for random number generation and MPI for parallel processing.

Usage:
Ising_OPV.exe takes one input argument that is the filename of the input parameter file.
An example of the parameter file is provided with parameters_default.txt
As an example, to run create one morphology on a single processor, the command is:
    Ising_OPV.exe parameters.txt
To run in a parallel processing environment using MPI and create 10 morphologies on 10 processors, the command is:
    mpiexec -n 10 Ising_OPV.exe parameters.txt

Output:
Ising_OPV will create several output files.  One text file will be created for each morphology generated that stores the data for that morphology.  In addition another text file, analysis_summary.txt
