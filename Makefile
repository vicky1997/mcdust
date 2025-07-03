F90=gfortran
H5F=h5fc
FFLAGS=-ffpe-trap=zero,overflow,invalid -fopenmp -Wall -W -ffree-form -g -fbounds-check -fbacktrace -fdefault-real-8 -fcheck=all -O2
#FFLAGS=-O2 -fdefault-real-8 -ffree-form -fopenmp

#add preprocessors:
CFLAGS=
# setup:
ifndef SETUP_FILE
	$(error Specify SETUP_FILE. Example: make SETUP_FILE=default)
 endif
# setuptest:
# 	SETUP_FILE=tests

OPT_FILE ?= setups/$(SETUP_FILE)/prepocs.opt
-include $(OPT_FILE)

SRC_DIR=src
OBJ_DIR=obj/$(SETUP_FILE)
TEST_DIR=unit_tests

EXECUTABLE = $(SETUP_FILE)
TEST=test1
vpath %.F90 $(SRC_DIR) $(TEST_DIR)


all: $(EXECUTABLE)
test: $(TEST)

$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)
	$(H5F) $(FFLAGS) $(CFLAGS) -J$(OBJ_DIR) -c $< -o $@

$(EXECUTABLE): $(OBJ_DIR)/main.o $(OBJ_DIR)/timestep.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o 
	$(H5F) $(FFLAGS) $(CFLAGS) $(LDFLAGS) $? -o $@

# setup:
# 	ifndef SETUP_FILE
# 		$(error Specify SETUP_FILE. Example: make SETUP_FILE=default)
# 	 endif
# setuptest:
# 	SETUP_FILE=tests

#test: 
#	$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)
#		$(H5F) $(FFLAGSTEST) $(CFLAGS) -J$(OBJ_DIR) -c $< -o $@

$(TEST): $(OBJ_DIR)/testsuite.o $(OBJ_DIR)/timestep.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o 
	$(H5F) $(FFLAGS) $(CFLAGS) $(LDFLAGS) $? -o $@


$(OBJ_DIR)/constants.o: $(SRC_DIR)/constants.F90
$(OBJ_DIR)/types.o: $(SRC_DIR)/types.F90
$(OBJ_DIR)/parameters.o: $(SRC_DIR)/parameters.F90 $(OBJ_DIR)/constants.o
$(OBJ_DIR)/discstruct.o: $(SRC_DIR)/discstruct.F90 $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/advection.o: $(SRC_DIR)/advection.F90 $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/types.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/initproblem.o: $(SRC_DIR)/initproblem.F90 $(OBJ_DIR)/advection.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/types.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o
$(OBJ_DIR)/mrgrnk.o: $(SRC_DIR)/mrgrnk.F90 $(OBJ_DIR)/types.o
$(OBJ_DIR)/parallel_sort.o: $(SRC_DIR)/parallel_sort.F90 $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/types.o
$(OBJ_DIR)/grid.o: $(SRC_DIR)/grid.F90 $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem.o
$(OBJ_DIR)/collisions.o: $(SRC_DIR)/collisions.F90 $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/types.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/timestep.o: $(SRC_DIR)/timestep.F90 $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/hdf5output.o: $(SRC_DIR)/hdf5output.F90 $(OBJ_DIR)/grid.o $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.F90 $(OBJ_DIR)/timestep.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/testsuite.o: $(TEST_DIR)/testsuite.F90 $(OBJ_DIR)/timestep.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem.o $(OBJ_DIR)/discstruct.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)


clean:
	rm -rf $(OBJ_DIR)/ $(EXECUTABLE)
cleantest:
	rm -rf $(OBJ_DIR)/ $(TEST)
cleancov:
	rm -rf *.gcno *.gcda *.gcov *.o *.mod *.dat *.out *.err $(EXECUTABLE)
