F90=gfortran
FFLAGS=-ffpe-trap=zero,overflow,invalid -Wall -W -ffree-form -std=gnu -fopenmp -fimplicit-none -g -fbounds-check -fbacktrace -fdefault-real-8
FFLAGSCOV=-O2 -fdefault-real-8 -ffree-form -fopenmp -pg --coverage
#FFLAGS=-O2 -fdefault-real-8 -ffree-form -fopenmp
%.o: %.f90

	$(F90) $(FFLAGS) -c $<

2DMC: main.o timestep.o output.o collisions.o advection.o grid.o parallel_sort.o mrgrnk.o initproblem.o discstruct.o parameters.o types.o constants.o
	$(F90) $(FFLAGS) $(LDFLAGS) $? -o $@

constants.o:   constants.f90
types.o:       types.f90
parameters.o:  parameters.f90 constants.o
discstruct.o:  discstruct.f90 parameters.o constants.o
advection.o:   advection.f90 discstruct.o types.o parameters.o constants.o
initproblem.o: initproblem.f90 advection.o constants.o types.o discstruct.o parameters.o
mrgrnk.o:      mrgrnk.f90 types.o
parallel_sort.o: parallel_sort.f90 mrgrnk.o types.o
grid.o:        grid.f90 initproblem.o types.o constants.o parallel_sort.o mrgrnk.o
collisions.o:  collisions.f90 initproblem.o grid.o types.o discstruct.o parameters.o constants.o
output.o:      output.f90 grid.o initproblem.o discstruct.o parameters.o types.o constants.o
timestep.o:    timestep.f90 advection.o grid.o discstruct.o parameters.o types.o constants.o
main.o:        main.f90 timestep.o output.o collisions.o advection.o grid.o parallel_sort.o mrgrnk.o initproblem.o discstruct.o parameters.o types.o constants.o

clean:
	rm -rf *.o *.mod *.dat *.out *.err 2DMC 

cleancov:
	rm -rf *.gcno *.gcda *.gcov *.o *.mod *.dat *.out *.err 2DMC
