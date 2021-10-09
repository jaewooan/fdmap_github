#!/bin/bash
#
#SBATCH --partition=twohour
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --output=gp229_%j.out
#SBATCH --error=gp229_%j.err
#SBATCH --time=01:00:00

module purge
ml compilers/intel16
ml openmpi/openmpi-2.0.1_intel16

date
mpirun fdmap gp229/gp229sw3ex1.in
date

