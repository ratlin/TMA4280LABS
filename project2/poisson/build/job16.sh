#!/bin/bash
#PBS -N proj2
#PBS -q training
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:15:00
#PBS -l select=2:ncpus=20:mpiprocs=16:ompthreads=4

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
#export OMP_SCHEDULE static
mpirun poisson 1
mpirun poisson 2
mpirun poisson 4
mpirun poisson 8
mpirun poisson 16
mpirun poisson 32
mpirun poisson 64
mpirun poisson 128
mpirun poisson 256
mpirun poisson 512
mpirun poisson 1024
mpirun poisson 2048
mpirun poisson 4096
mpirun poisson 8192
