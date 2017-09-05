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


ifeq ($(USE_MPI),TRUE)
    CXXFLAGS = -Wno-deprecated-gpu-targets -x cu --std=c++11 -ccbin=mpic++ -O3
    CFLAGS   = -Wno-deprecated-gpu-targets -x c -ccbin=mpicc -c99 -O3
else
    CXXFLAGS = -Wno-deprecated-gpu-targets -x cu -std=c++11 -ccbin=g++ -O3 -dc
    CFLAGS   = -Wno-deprecated-gpu-targets -x c -c99 -ccbin=gcc -O3 -dc
    # CXXFLAGS = -Wno-deprecated-gpu-targets -x cu -ccbin=pgc++ -O3 -dc
    # CFLAGS   = -Wno-deprecated-gpu-targets -x c -ccbin=pgcc -O3 -dc
endif 

# other options 
# verbose: -Xptxas=-v
# CXXFLAGS += -Xptxas -dlcm=ca
# CFLAGS += -Xptxas -dlcm=ca

########################################################################

pgi_version := $(shell $(CXX) -V 2>&1 | grep 'target')

COMP_VERSION := $(pgi_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  #TODO: change these to pgc++ flags
  CXXFLAGS += -G -Xptxas=-v -Xcompiler='-g -O0 -Mbounds'
  CFLAGS   += -G -Xptxas=-v -Xcompiler='-g -O0 -Mbounds'
  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  CXXFLAGS += -lineinfo -rdc=true -Xcompiler='-g -O3'
  CFLAGS   += -lineinfo -rdc=true -Xcompiler='-g -O3'
  FFLAGS   += -gopt -fast -Mcuda=cuda8.0 -Mnomain -Mcuda=lineinfo -Mcuda=rdc
  F90FLAGS += -gopt -fast -Mcuda=cuda8.0 -Mnomain -Mcuda=lineinfo -Mcuda=rdc

  # CXXFLAGS += -lineinfo -rdc=true -Xcompiler='-gopt -fast --c++11 -Mcuda=cuda8.0'
  # CFLAGS   += -lineinfo -rdc=true -Xcompiler='-gopt -fast -c99 -Mcuda=cuda8.0'

endif

########################################################################


F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  # CXXFLAGS += -Xcompiler='-mp=nonuma -Minfo=mp -noacc'
  # CFLAGS   += -Xcompiler='-mp=nonuma -Minfo=mp -noacc'
  CXXFLAGS += -Xcompiler='-fopenmp'
  CFLAGS   += -Xcompiler='-fopenmp'
  FFLAGS   += -mp=nonuma -Minfo=mp -noacc
  F90FLAGS += -mp=nonuma -Minfo=mp -noacc
  override XTRALIBS += -lgomp
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
