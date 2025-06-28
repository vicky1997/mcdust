## Installation

`git clone git@github.com:vicky1997/mcdust.git`

## Prerequisites

`gfortran` `hdf5-serial` `python`

## To compile the code: 

`make SETUP_FILE=default`

This uses the setup files from the default run in `/setups/default/` 
The setup files can be modified to suit your needs and can be added to the `setups` directory.

parameter file: setups/default/setup.par

## To run the code: 

export OMP_NUM_THREADS=number_of_CPUs

./default setup.par

For more details: [mcdust Documentation](https://mcdust.readthedocs.io/en/latest/)
