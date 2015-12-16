# Default build type

#BUILD = debug
BUILD = production

# Default executables

Exe = fdmap
IOExe = testio

# Compiler options for different machines

HOST = $(shell hostname)
UNAME = $(shell uname)

# Mac (using gfortran and openmpi)

ifeq ($(findstring Darwin,$(UNAME)),Darwin)
 F95 = mpif90 -fdefault-real-8 -fdefault-double-8
 LD = mpif90
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -Wall -Wextra -Wconversion -fbounds-check -fbacktrace \
	-fimplicit-none -std=f2003
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O5 -Wuninitialized
 endif
 LDFLAGS = $(F95FLAGS)
 LIBS = -framework Accelerate
 INCL = -fall-intrinsics
endif

# Stanford CEES

ifeq ($(findstring cees,$(HOST)),cees)
 F95 = mpif90 -r8 -i4 -i_dynamic
 LD = mpif90
 ifeq ($(BUILD),debug)
  F95FLAGS =  -g -check all -check noarg_temp_created -warn all \
	-traceback -ftrapuv -fpe0 -fp-stack-check -fltconsistency -std03
 endif
 ifeq ($(BUILD),profile)
  F95FLAGS = -g -pg
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -g -O2
#  F95FLAGS = -fast
 endif
 LDFLAGS = $(F95FLAGS)
 MKLPATH = /usr/local/intel/mkl
 LIBS =	-L$(MKLPATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -Wl,--end-group
 INCL = 
endif

# TACC Ranger

ifeq ($(findstring ranger,$(HOST)),ranger)
 F95 = mpif90 -r8 -i4
 LD = mpif90
 ifeq ($(BUILD),debug)
  F95FLAGS =  -g -check all -check noarg_temp_created -warn all \
	-traceback -ftrapuv -fpe0 -fp-stack-check -fltconsistency -std03
 endif
 ifeq ($(BUILD),profile)
  F95FLAGS = -g -pg
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -g -O2 -xW
 endif
 LDFLAGS = -g
 LIBS =	-Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -lmkl -lguide
 INCL = 
endif

# KAUST Shaheen BlueGene/P

ifeq ($(findstring fen,$(HOST)),fen)
 F95 = mpixlf95_r -qautodbl=dbl4
 LD = mpixlf95_r
 ifeq ($(BUILD),debug)
  F95FLAGS = -g
 endif
 ifeq ($(BUILD),production)
  F95FLAGS =
 endif
 LDFLAGS = -g
 LIBS = -L/bgsys/ibm_essl/sles10/prod/opt/ibmmath/lib -lesslbg
 INCL = -I/opt/share/include/mpi/ibm
endif

# Files

Files = basal_traction.f90 boundaries.f90 checkpoint.f90 domain.f90 energy.f90 \
	interfaces.f90 fd.f90 fd_coeff.f90 friction.f90 fields.f90 \
	geometry.f90 grid.f90 hydrofrac.f90 io.f90 main.f90 material.f90 \
	mms.f90 mpi_routines.f90 mpi_routines2d.f90 output.f90 plastic.f90 \
	rates.f90 rates_heterogeneous.f90 source.f90 thermpres.f90 \
	time_step.f90 tsunami.f90 utilities.f90

IOFiles = io.f90 mpi_routines.f90 mpi_routines2d.f90 testio.f90

Obs = 	$(Files:.f90=.o)

IOObs = $(IOFiles:.f90=.o)

all: code iotest

code: $(Exe)

iotest: $(IOExe)

$(Exe): $(Obs)
	$(LD) $(LDFLAGS) -o $(Exe) \
	$(Obs) $(LIBS)

$(IOExe): $(IOObs)
	$(LD) $(LDFLAGS) -o $(IOExe) \
	$(IOObs) $(LIBS)

%.o : %.f90
	$(F95) $(F95FLAGS) -c $< -o $@ $(INCL)

clean:
	rm -f *.o *.mod $(Exe) $(IOExe)

intel:
	sed 's/!use ifport/use ifport/g' mms.f90 > temp.f90 && \
	mv temp.f90 mms.f90

notintel:
	sed 's/ use ifport/ !use ifport/g' mms.f90 > temp.f90 && \
	mv temp.f90 mms.f90

.SUFFIXES: .o .f90

# DO NOT DELETE THIS LINE - used by make depend
boundaries.o: fields.o geometry.o grid.o io.o mms.o mpi_routines.o tsunami.o

checkpoint.o: domain.o fields.o interfaces.o io.o mpi_routines.o

domain.o: boundaries.o fd_coeff.o fields.o grid.o interfaces.o io.o material.o
domain.o: mpi_routines.o mpi_routines2d.o source.o utilities.o

energy.o: domain.o fields.o geometry.o grid.o material.o mpi_routines.o
energy.o: mpi_routines2d.o

fd.o: fd_coeff.o io.o

fd_coeff.o: io.o

fields.o: fd.o geometry.o grid.o io.o material.o mms.o mpi_routines.o
fields.o: mpi_routines2d.o utilities.o

friction.o: geometry.o io.o material.o mms.o mpi_routines.o utilities.o

geometry.o: io.o

grid.o: fd.o fd_coeff.o geometry.o io.o mpi_routines.o mpi_routines2d.o
grid.o: utilities.o

hydrofrac.o: fd.o fd_coeff.o io.o mms.o mpi_routines.o mpi_routines2d.o

interfaces.o: fields.o friction.o geometry.o grid.o hydrofrac.o io.o mms.o
interfaces.o: mpi_routines.o mpi_routines2d.o thermpres.o

io.o: mpi_routines.o

main.o: checkpoint.o domain.o io.o mpi_routines.o output.o time_step.o
main.o: utilities.o

material.o: io.o mpi_routines.o mpi_routines2d.o

mms.o: io.o material.o


mpi_routines2d.o: io.o mpi_routines.o

output.o: domain.o io.o mpi_routines.o utilities.o

plastic.o: fields.o grid.o io.o material.o

rates.o: fd_coeff.o material.o

rates_heterogeneous.o: fd_coeff.o

source.o: fields.o grid.o io.o material.o mms.o mpi_routines.o utilities.o

testio.o: io.o mpi_routines.o mpi_routines2d.o

thermpres.o: io.o mpi_routines.o

time_step.o: domain.o energy.o fields.o grid.o interfaces.o io.o material.o
time_step.o: mpi_routines2d.o output.o plastic.o rates.o rates_heterogeneous.o
time_step.o: source.o


utilities.o: io.o

