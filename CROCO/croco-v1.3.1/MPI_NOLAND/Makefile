NETCDF_ROOT=/usr/local
NETCDF_INC=$(NETCDF_ROOT)/include
NETCDF_LIB=$(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
#
F90=gfortran
F90FLAGS= -O3 -fdefault-real-8
LDFLAGS= 
#
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90))
EXEC=mpp_optimiz

$(EXEC): $(SRCS)
	$(F90) $(LDFLAGS) -o $@ $< -L$(NETCDF_LIB)

%.o: %.f90
	$(F90) $(F90FLAGS) -c  $< -I$(NETCDF_INC) 

clean:
	rm -f *.o *.mod $(EXEC)
