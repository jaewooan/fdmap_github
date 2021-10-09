#!/bin/bash
#
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --output=gp229_%j.out
#SBATCH --error=gp229_%j.err
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=01:00:00

module purge
source /oak/stanford/schools/ees/share/serc_env.sh
ml load gcc-cees/10.1.0
ml load openmpi-cees/4.0.5
ml load openblas-cees/0.3.14

date
mpirun fdmap gp229/gp229sw3ex1.in
date

