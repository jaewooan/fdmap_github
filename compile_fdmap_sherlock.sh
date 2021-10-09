#!/bin/bash
#
#SBATCH --partition=serc
#SBATCH --ntasks=1
#SBATCH --output=fdmap_compile_%j.out
#SBATCH --error=fdmap_compile_%j.err
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=5
#
# compile_fdmap-sherlock.sh
# principal authors: Mark R. Yoder, Stanford Research Computing
#  Randy "RC" White, Stanford Research Computing
#
module purge
#
#run on Rome node
source /oak/stanford/schools/ees/share/serc_env.sh
ml load gcc-cees/10.1.0
ml load openmpi-cees/4.0.5
ml load openblas-cees/0.3.14
#
#echo "modules: "
#module list
#
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
#export FCFLAGS="-fdefault-real-8 -fdefault-double-8 -g -Wall -Wextra -Wconversion -fbounds-check -fbacktrace -fimplicit-none -ffpe-summary=all -ffpe-trap=invalid,zero,overflow -std=f2003 -Wno-compare-reals"
#export LDFLAGS="${FCFLAGS}"
#
# production:
export FCFLAGS="-fdefault-real-8 -fdefault-double-8 -O5 -Wuninitialized"
export LDFLAGS="-O5 -Wuninitialized"
#
#
export INCL="-fall-intrinsics "
LIBS="-Wl,-rpath -Wl,--enable-new-dtags `pkg-config --libs-only-L ompi-fort` `pkg-config --libs-only-l ompi-fort` "
#
# gcc/blas:
export LIBS="${LIBS} -L/oak/stanford/schools/ees/share/cees/software/spack_/zen2/spack/opt/spack/linux-centos7-zen2/gcc-10.1.0/openblas-0.3.14-pvdl3wuymvfqnsmegbupg6rkevxf5bfg/lib -lopenblas -BStatic "

#
make -f Makefile.noflags clean
make -f Makefile.noflags



