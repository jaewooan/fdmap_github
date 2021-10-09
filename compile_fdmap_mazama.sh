#!/bin/bash
#
#SBATCH --partition=twohour
#SBATCH --ntasks=1
#SBATCH --output=fdmap_compile_%j.out
#SBATCH --error=fdmap_compile_%j.err
#SBATCH --time=5
#
module purge
ml compilers/intel16
ml openmpi/openmpi-2.0.1_intel16
export FC=mpifort
export CC=mpicc
export CXX=mpicxx
export MPIFC=mpifort
export MPICC=mpicc
export MPICXX=mpicxx
#
export LD="${MPIFC}"
#
# debug:
#export FCFLAGS="-r8 -i4 -g -O0 -init=snan,arrays -traceback -check all -check noarg_temp_created \
#	-warn all -ftrapuv -fpe0 -fp-stack-check -std03"
#export LDFLAGS="${FCFLAGS}"
#
# production:
export FCFLAGS="-r8 -i4 -g -O2"
export LDFLAGS="${FCFLAGS}"

#
export INCL=""
export MKLPATH=""
export LIBS="-L/usr/local/intel/mkl -Wl,--start-group -lmkl_intel_lp64 \
 -lmkl_sequential -lmkl_core -Wl,--end-group"

#
make -f Makefile.noflags clean
make -f Makefile.noflags



