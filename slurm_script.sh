#!/bin/bash
#SBATCH -J Ising_OPV # Job name
#SBATCH -p partition_name
#SBATCH -n 8 # Number of tasks
#SBATCH -t 00:30:00 # Maximum walltime
#SBATCH --cpus-per-task=1

ParameterNum=default

# Setup job directory
mkdir $SLURM_JOB_ID
cd $SLURM_JOB_ID
cp ../parameters_$ParameterNum.txt ./parameters_$ParameterNum.txt

# Execute Simulation
mpiexec -n 8 ../Ising_OPV.exe parameters_$ParameterNum.txt > output.txt

# Save Results and Cleanup
cd ..
tar -zcf $SLURM_JOB_ID.tar.gz $SLURM_JOB_ID
mkdir -p results
cp $SLURM_JOB_ID.tar.gz ./results/$SLURM_JOB_ID.tar.gz
rm -f $SLURM_JOB_ID.tar.gz
rm -rf $SLURM_JOB_ID
