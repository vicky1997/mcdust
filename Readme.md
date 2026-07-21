Markdown:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.10140/status.svg)](https://doi.org/10.21105/joss.10140)

HTML:
<a style="border-width:0" href="https://doi.org/10.21105/joss.10140">
  <img src="https://joss.theoj.org/papers/10.21105/joss.10140/status.svg" alt="DOI badge" >
</a>

reStructuredText:
.. image:: https://joss.theoj.org/papers/10.21105/joss.10140/status.svg
   :target: https://doi.org/10.21105/joss.10140
   
## Installation

`git clone git@github.com:vicky1997/mcdust.git`

## Prerequisites
 - A Linux/Unix based OS
 - A Fortran compiler (preferably `gfortran` with `OpenMP`)
 - A Fortran-OpenMP installation
 - A HDF5 serial installation ()
 - A Python installation (if you want to use the pre-written post-processing scripts)

`gfortran` `hdf5-serial` `python` `OpenMP`

To install the required software in Ubuntu(this requires root permissions):
`sudo apt-get install gfortran`
`sudo apt-get install libhdf5-serial-dev`

Python is not required to run the code. But if you want to use the routines to read/write data from the simulation you will need a python installation.
v 
## To compile the code: 

`make SETUP_FILE=default`

This uses the setup files from the default run in `setups/default/` 
The setup files can be modified to suit your needs and can be added to the `setups` directory. An executable in the name of the setup file will be created in the root directory which is `default` in this case.

parameter file: setups/default/setup.par

To run the unit tests:

`make test SETUP_FILE=tests`
## To run the code: 

`export OMP_NUM_THREADS=number_of_CPUs`

`./default setups/default/setup.par`

To clear the run

`make clean SETUP_FILE=default`

## Unit tests
To run the unit tests in `unit_tests/testsuite.F90`

`make test SETUP_FILE=tests`

`./test1 setups/tests/setup.par`

To clear the run

`make cleantest SETUP_FILE=tests`

For more details: [mcdust Documentation](https://mcdust.readthedocs.io/en/latest/)

## Acknowledgments
The authors acknowledge funding from the European Union under the European Union’s Horizon Europe Research \& Innovation Programme 101040037 (PLANETOIDS). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.