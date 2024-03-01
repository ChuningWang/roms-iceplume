#
# Include file for OASIS3 Makefile for ADA 
# using mpiifort
#
###############################################################################
#
# CHAN	: communication technique used in OASIS3 (MPI1/MPI2/NONE)
CHAN            = MPI1
#
# Paths for libraries, object files and binaries
#
# COUPLE	: path for oasis3-mct main directory
COUPLE          = $(HOME)/oasis/OASIS3-MCT_3.0_branch/oasis3-mct
#
# ARCHDIR       : directory created when compiling
ARCHDIR         = $(HOME)/oasis/compile_oasis3-mct_v3_netcdf4.1.3
#ARCHDIR         = $(HOME)/oasis/compile_oasis3-mct_netcdf_serial
#
# MPI library
MPIDIR      = 
MPIBIN      = 
MPI_INCLUDE = 
MPILIB      = 
#
# NETCDF library
NETCDF_INCLUDE  = $(shell $(NETCDF_CONFIG) --fflags)
NETCDF_LIBRARY  = $(shell $(NETCDF_CONFIG) --flibs)
#
# Compiling and other commands
MAKE        = gmake
F90         = mpiifort
F           = $(F90)
f90         = $(F90)
f           = $(F90)
CC          = mpiicc 
LD          = mpiifort
AR          = ar
ARFLAGS     = r
#
# CPP keys and compiler options
# 
#CPPDEF    = -Duse_libMPI -Duse_comm_$(CHAN)
CPPDEF     = -Duse_libMPI -Duse_comm_$(CHAN)
#
# -g is necessary in F90FLAGS and LDFLAGS for pgf90 versions lower than 6.1
# 
F90FLAGS_1  = -O2 -I. 
#F90FLAGS_1  = -O3 -xAVX -I. -fp-model precise
f90FLAGS_1  = $(F90FLAGS_1)
FFLAGS_1    = $(F90FLAGS_1)
fFLAGS_1    = $(F90FLAGS_1)
CCFLAGS_1   = -O2 -xAVX
LDFLAGS     = 
#
###################
#
# Additional definitions that should not be changed
#
FLIBS		= $(NETCDF_LIBRARY)
# BINDIR        : directory for executables
BINDIR          = $(ARCHDIR)/bin
# LIBBUILD      : contains a directory for each library
LIBBUILD        = $(ARCHDIR)/build/lib
# INCPSMILE     : includes all *o and *mod for each library
INCPSMILE       = -I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/mct -I$(LIBBUILD)/scrip
#
F90FLAGS  = $(F90FLAGS_1) $(INCPSMILE) $(CPPDEF) -I$(NETCDF_INCLUDE)
f90FLAGS  = $(f90FLAGS_1) $(INCPSMILE) $(CPPDEF) -I$(NETCDF_INCLUDE)
FFLAGS    = $(FFLAGS_1) $(INCPSMILE) $(CPPDEF) -I$(NETCDF_INCLUDE)
fFLAGS    = $(fFLAGS_1) $(INCPSMILE) $(CPPDEF) -I$(NETCDF_INCLUDE)
CCFLAGS   = $(CCFLAGS_1) $(INCPSMILE) $(CPPDEF) -I$(NETCDF_INCLUDE)
#
#############################################################################
