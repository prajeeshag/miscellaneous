
FC = gfortran

#NETCDF = /gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.1.3
NETCDF = /usr/local

OPTS = -O3 --default-real-8

#OPTS = -03 -r8

laplace.exe: laplace_soln_implicit.o netcdf_write.o
	$(FC) $(OPTS) -o laplace.exe netcdf_write.o  laplace_soln_implicit.o -llapack -lblas -L$(NETCDF)/lib -lnetcdf -lnetcdff

netcdf_write.o: netcdf_write.f90
	$(FC) $(OPTS) -I$(NETCDF)/include -c netcdf_write.f90

laplace_soln_implicit.o: laplace_soln_implicit.f90 netcdf_write.o
	$(FC) $(OPTS) -c laplace_soln_implicit.f90

clean:
	rm -f *.o *.mod laplace.exe
