Quickstart
=====

Pre-requisites:

    * A Linux/Unix based OS

    * A Fortran compiler (preferably :code:`gfortran`), 
    
    * A Fortran-OpenMP installation

    * A HDF5 serial installation
    
    * A Python installation (if you want to use the pre-written post-processing scripts)

Now one can clone the GitHub repository to use the code.

:code:`git clone git@github.com:vicky1997/mcdust.git`

First run
+++++++++
#. Make the :code:`/setups/default` setup 
    :code:`make SETUP_FILE=default`
#. Set the number of threads for OpenMP parallelisation.
    :code:`export OMP_NUM_THREADS=$OMP_NUM_THREADS`
#. Run the code
    | :code:`./default /setups/default/setup.par`

Restart
+++++++

To (re)start the simulation from a specific snapshot, follow the steps.

#. Copy the snapshot of your choice to the setup directory from where you want to restart
    :code:`cp your/data/directory/foo.h5 setups/restartrun/`
#. Rename the file to :code:`restart.h5`
    :code:`mv foo.h5 restart.h5`
#. Set the :code:`restart` condition in the :code:`setup.par` file to :code:`.true.`

