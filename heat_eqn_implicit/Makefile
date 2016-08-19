
FC = gfortran
NETCDF = /usr/local
OPTS = -O3 --default-real-8 -g

#FC = ifort
#NETCDF = /gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.1.3
#OPTS = -O3 -r8 -g -traceback

heat_eqn_implicit.exe: heat_eqn_implicit.o netcdf_write.o
	$(FC) $(OPTS) -o heat_eqn_implicit.exe netcdf_write.o heat_eqn_implicit.o -llapack -lblas -L$(NETCDF)/lib -lnetcdf -lnetcdff

netcdf_write.o: ../netcdf_write/netcdf_write.f90
	$(FC) $(OPTS) -I$(NETCDF)/include -c ../netcdf_write/netcdf_write.f90

heat_eqn_implicit.o: heat_eqn_implicit.F90 netcdf_write.o
	$(FC) $(OPTS) -c heat_eqn_implicit.F90

clean:
	rm -f *.o *.mod heat_eqn_implicit.exe out*.nc out*.dat
