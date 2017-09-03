#
# Setup for compiling the CUDA version of AMReX with
# CUDA C
# Assumes you have set USE_CUDA=TRUE, and have
# set the variables PGI_PATH to the root PGI
# directory and CUDA_PATH to the root CUDA directory.
#
CXX = nvcc
CC  = nvcc
# FC  = pgfortran
# F90 = pgfortran
FC  = gfortran
F90 = gfortran


ifeq ($(USE_MPI),TRUE)
    CXXFLAGS = -Wno-deprecated-gpu-targets -x cu --std=c++11 -ccbin=mpic++ -O3
    CFLAGS   = -Wno-deprecated-gpu-targets -x c -ccbin=mpicc -c99 -O3
else
    CXXFLAGS = -Wno-deprecated-gpu-targets -ccbin=g++ -dc
    CFLAGS   = -Wno-deprecated-gpu-targets -ccbin=gcc -dc
endif 

# other options 
# verbose: -Xptxas=-v
# CXXFLAGS += -Xptxas -dlcm=ca
# CFLAGS += -Xptxas -dlcm=ca
FFLAGS   =
F90FLAGS =

########################################################################

pgi_version := $(shell $(CXX) -V 2>&1 | grep 'target')

COMP_VERSION := $(pgi_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  CXXFLAGS += -G -Xptxas=-v -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'
  CFLAGS   += -G -Xptxas=-v -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'
  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  CXXFLAGS += -Xcompiler='-g -O3 -std=c++11'
  CFLAGS   += -Xcompiler='-g -O3 -std=gnu99'
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3
  # FFLAGS   += -gopt -fast
  # F90FLAGS += -gopt -fast

endif

########################################################################


# F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
# FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend
FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir)
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir) -fimplicit-none

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  # for g++ and gcc
  CXXFLAGS += -Xcompiler='-fopenmp'
  CFLAGS += -Xcompiler='-fopenmp'
  # # for pgfortran
  # FFLAGS   += -Minfo=mp
  # F90FLAGS += -Minfo=mp
  # GENERIC_COMP_FLAGS +=
    FFLAGS   += -fopenmp 
    F90FLAGS += -fopenmp 
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_COMP_FLAGS += -acc -Minfo=acc -ta=nvidia -lcudart -mcmodel=medium
else
  GENERIC_COMP_FLAGS += 
endif

# TODO:
# actually invoking this Makefile already indicates that USE_CUDA=TRUE
ifeq ($(USE_CUDA),TRUE)
  # CXXFLAGS +=
  # CFLAGS   +=
  # FFLAGS   += -Mnomain
  # F90FLAGS += -Mnomain

  # override XTRALIBS += -lstdc++
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################
# ask gfortran the name of the library to link in.  First check for the
# static version.  If it returns only the name w/o a path, then it
# was not found.  In that case, ask for the shared-object version.
gfortran_liba  := $(shell $(F90) -print-file-name=libgfortran.a)
gfortran_libso := $(shell $(F90) -print-file-name=libgfortran.so)

ifneq ($(gfortran_liba),libgfortran.a)  # if found the full path is printed, thus `neq`.
  LIBRARY_LOCATIONS += $(dir $(gfortran_liba))
else
  LIBRARY_LOCATIONS += $(dir $(gfortran_libso))
endif

# Because we do not have a Fortran main

ifeq ($(which_computer),$(filter $(which_computer),summit))
override XTRALIBS += -pgf90libs -L /sw/summitdev/gcc/5.4.0new/lib64/ -latomic
else
override XTRALIBS += -lgfortran -lquadmath -lstdc++
# override XTRALIBS += -pgf90libs -latomic -lquadmath
endif
