

FC = gfortran
NETCDF = /usr/local
debug = -fcheck=bounds
OPTS = $(debug) -O3 --default-real-8 -g -I$(NETCDF)/include -L$(NETCDF)/lib -lnetcdf -lnetcdff


#FC = ifort
#NETCDF = /gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.1.3
#OPTS = -O3 -r8 -g -traceback

OBJS = long_surf_gravity_waves.o netcdf_write.o

shallow.exe: $(OBJS)
	$(FC) $(OBJS) $(OPTS) -o $@ 

netcdf_write.o: ../netcdf_write/netcdf_write.f90
	$(FC) $(OPTS) -c $<

long_surf_gravity_waves.o: long_surf_gravity_waves.f90 netcdf_write.o
	$(FC) $(OPTS) -c $< 

clean:
	rm -f *.o *.mod *.exe out*.nc out*.dat
