# Set default build type

#BUILD = debug
BUILD = production

Exe = fdmap
IOExe = testio

# Compiler options for different machines

HOST = $(shell hostname)
UNAME = $(shell uname)
USER = $(shell users)

# Mac Pro

ifeq ($(findstring Darwin,$(UNAME)),Darwin)
# F95 = openmpif90 -r8 -i4
# LD = openmpif90
# ifeq ($(BUILD),debug)
#  F95FLAGS = -g -Wall -Wextra -Wno=165 -fbounds-check -ftrace=full \
#	-fimplicit-none -finteger=100000 -freal=nan -fpointer=invalid -std=f2003
# endif
 F95 = openmpif90 -fdefault-real-8 -fdefault-double-8
 LD = openmpif90
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -Wall -Wextra -Wconversion -fbounds-check -fbacktrace \
	-fimplicit-none -std=f2003
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O5 -Wuninitialized
 endif
 LDFLAGS = $(F95FLAGS)
 LIBS = -framework vecLib
# INCL = -fintrinsic-extensions
 INCL = -fall-intrinsics
endif

# TACC Ranger

ifeq ($(findstring ranger,$(HOST)),ranger)
 ifeq ($(findstring edunham,$(USER)),edunham)
  Exe = /work/01331/edunham/fdmap/fdmap
  IOExe = /work/01331/edunham/fdmap/testio
 else
  Exe = fdmap
  IOExe = testio
 endif
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

# Stanford CEES

ifeq ($(findstring cees,$(HOST)),cees)
 # ifeq ($(findstring edunham,$(USER)),edunham)
 #  Exe = /data/dunham1/edunham/flaunch/fdmap
 #  IOExe = /data/dunham1/edunham/flaunch/testio
 # else
  Exe = fdmap
  IOExe = testio
 # endif
 F95 = mpif90 -r8 -i4 -i_dynamic
 LD = mpif90
 ifeq ($(BUILD),debug)
  F95FLAGS =  -g -check all -check noarg_temp_created -warn all \
	-traceback -ftrapuv -fpe0 -fp-stack-check -fltconsistency -std03
#	-traceback -ftrapuv -ftz -std03
 endif
 ifeq ($(BUILD),profile)
  F95FLAGS = -g -pg
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -g -O2
#  F95FLAGS = -g -O3 -ipo -ip
 endif
 LDFLAGS = $(F95FLAGS)
 MKLPATH = /usr/local/INTEL/Compiler/11.1/069/mkl
 LIBS =	-L$(MKLPATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -Wl,--end-group -lpthread
 INCL = 
endif

# Harvard BlueGene/L

ifeq ($(findstring bgfen,$(HOST)),bgfen)
 F95 = blrts_xlf95 -qautodbl=dbl4
 LD = blrts_xlf95
 ifeq ($(BUILD),debug)
  F95FLAGS =  -O0 -qarch=440 -qmaxmem=64000 -g \
	-qflttrap=overflow:underflow:zerodivide:invalid:enable
 endif
 ifeq ($(BUILD),production)
  F95FLAGS =  -O5 -qstrict -qarch=440d -qtune=440 -qmaxmem=-1 -g
 endif
 LDFLAGS = -g
 LIBS =	-L/bghome/edunham/software/lib -lmpiP \
	-L/bgl/BlueLight/ppcfloor/bglsys/lib \
        -lchkpt.rts -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts
 INCL = -I/bgl/BlueLight/ppcfloor/bglsys/include -I/bghome/edunham/include
endif

# Files

Files = boundaries.f90 checkpoint.f90 domain.f90 energy.f90 \
	erupt.f90 fault.f90 fd.f90 fd_coeff.f90 friction.f90 fields.f90 \
	geometry.f90 gravity.f90 grid.f90 io.f90 main.f90 material.f90 \
	mpi_routines.f90 mpi_routines2d.f90 output.f90 plastic.f90 rates.f90 \
	thermpres.f90 time_step.f90 utilities.f90

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
	sed 's/!use ifport/use ifport/g' fields.f90 > temp.f90 && \
	mv temp.f90 fields.f90

notintel:
	sed 's/ use ifport/ !use ifport/g' fields.f90 > temp.f90 && \
	mv temp.f90 fields.f90

.SUFFIXES: .o .f90

# DO NOT DELETE THIS LINE - used by make depend
boundaries.o: erupt.o fields.o friction.o geometry.o grid.o io.o material.o
boundaries.o: mpi_routines.o mpi_routines2d.o thermpres.o

rates.o: fd_coeff.o material.o

checkpoint.o: boundaries.o domain.o fields.o io.o mpi_routines.o

domain.o: boundaries.o fault.o fd_coeff.o fields.o gravity.o grid.o io.o
domain.o: material.o mpi_routines.o mpi_routines2d.o utilities.o

energy.o: domain.o fields.o grid.o material.o mpi_routines.o mpi_routines2d.o

erupt.o: fd.o fd_coeff.o io.o mpi_routines.o mpi_routines2d.o utilities.o

fault-old.o: boundaries.o fields.o friction.o geometry.o io.o material.o

fault.o: boundaries.o erupt.o fields.o friction.o geometry.o io.o material.o
fault.o: mpi_routines2d.o thermpres.o

fd.o: fd_coeff.o io.o

fd_coeff.o: io.o

fields.o: fd.o geometry.o grid.o io.o material.o mpi_routines.o
fields.o: mpi_routines2d.o utilities.o

friction.o: fields.o geometry.o io.o material.o mpi_routines.o utilities.o

geometry.o: io.o

gravity.o: fields.o io.o mpi_routines.o

grid.o: fd.o fd_coeff.o geometry.o io.o mpi_routines.o mpi_routines2d.o
grid.o: utilities.o

io.o: mpi_routines.o

main.o: checkpoint.o domain.o io.o mpi_routines.o output.o time_step.o
main.o: utilities.o

material.o: io.o mpi_routines.o


mpi_routines2d.o: io.o mpi_routines.o

output.o: domain.o io.o mpi_routines.o utilities.o

plastic.o: fields.o grid.o io.o material.o

testio.o: io.o mpi_routines.o mpi_routines2d.o

thermpres.o: io.o mpi_routines.o

time_step.o: boundaries.o rates.o domain.o energy.o fault.o fd.o
time_step.o: fd_coeff.o fields.o gravity.o grid.o io.o material.o
time_step.o: mpi_routines2d.o output.o plastic.o

utilities.o: io.o

