## Installation

`git clone git@github.com:vicky1997/mcdust.git`

## Prerequisites

`gfortran` `hdf5-serial` `python`

To install the required software in Ubuntu(this requires root permissions):
`sudo apt-get install gfortran`
`sudo apt-get install libhdf5-serial-dev`

Python is not required to run the code. But if you want to use the routines to read/write data from the simulation you will need a python installation.

## To compile the code: 

`make SETUP_FILE=default`

This uses the setup files from the default run in `/setups/default/` 
The setup files can be modified to suit your needs and can be added to the `setups` directory. An executable in the name of the setup file will be created in the root directory which is `default` in this case.

parameter file: setups/default/setup.par

To run the unit tests:

`make test SETUP_FILE=tests`
## To run the code: 

`export OMP_NUM_THREADS=number_of_CPUs`

`./default /setups/default/setup.par`

To clear the run

`make clean SETUP_FILE=default`

## Unit tests
To run the unit tests in `unit_tests/testsuite.F90`

`make test SETUP_FILE=tests`

To clear the run

`make cleantest SETUP_FILE=tests`

For more details: [mcdust Documentation](https://mcdust.readthedocs.io/en/latest/)

## Acknowledgments
The authors acknowledge funding from the European Union under the European Union’s Horizon Europe Research \& Innovation Programme 101040037 (PLANETOIDS). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.