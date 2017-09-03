#
# Setup for compiling the CUDA version of AMReX with
# CUDA C
# Assumes you have set USE_CUDA=TRUE, and have
# set the variables PGI_PATH to the root PGI
# directory and CUDA_PATH to the root CUDA directory.
#
CXX = nvcc
CC  = nvcc
FC  = pgfortran
F90 = pgfortran
# FC  = gfortran
# F90 = gfortran


ifeq ($(USE_MPI),TRUE)
    CXXFLAGS = -Wno-deprecated-gpu-targets -x cu --std=c++11 -ccbin=mpic++ -O3
    CFLAGS   = -Wno-deprecated-gpu-targets -x c -ccbin=mpicc -c99 -O3
else
    CXXFLAGS = -Wno-deprecated-gpu-targets -ccbin=pgc++ -dc
    CFLAGS   = -Wno-deprecated-gpu-targets -ccbin=pgcc -dc
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

  #TODO: change these to pgc++ flags
  CXXFLAGS += -G -Xptxas=-v -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'
  CFLAGS   += -G -Xptxas=-v -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'
  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  CXXFLAGS += -Xcompiler='-gopt -fast --c++11 -Mcuda=cuda8.0'
  CFLAGS   += -Xcompiler='-gopt -fast -c99 -Mcuda=cuda8.0'
  # FFLAGS   += -g -O3
  # F90FLAGS += -g -O3
  FFLAGS   += -gopt -fast -Mcuda=cuda8.0 -Mnomain -Mcuda=lineinfo
  F90FLAGS += -gopt -fast -Mcuda=cuda8.0 -Mnomain -Mcuda=lineinfo

endif

########################################################################


F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend
# FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir)
# F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir) -fimplicit-none

########################################################################

GENERIC_COMP_FLAGS =
ifeq ($(USE_OMP),TRUE)
  CXXFLAGS += -Xcompiler='-mp=nonuma -Minfo=mp -noacc'
  CFLAGS   += -Xcompiler='-mp=nonuma -Minfo=mp -noacc'
  FFLAGS   += -mp=nonuma -Minfo=mp -noacc
  F90FLAGS   += -mp=nonuma -Minfo=mp -noacc
endif


ifeq ($(USE_ACC),TRUE)
  GENERIC_COMP_FLAGS += -acc -Minfo=acc -ta=nvidia -lcudart -mcmodel=medium
else
  GENERIC_COMP_FLAGS += 
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)


ifeq ($(which_computer),$(filter $(which_computer),summit))
override XTRALIBS += -pgf90libs -L /sw/summitdev/gcc/5.4.0new/lib64/ -latomic -lstdc++
else
override XTRALIBS += -pgf90libs -latomic -lquadmath  -lstdc++
endif
