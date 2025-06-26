Quickstart
=====

Pre-requisites:

    * A Linux/Unix based OS

    * a fortran compiler (preferably gfortran), 
    
    * a fortran-OpenMP installation

    * a hdf5 serial installation
    
    * a Python installation (if you want to use the pre-written post-processing scripts)

Now one can clone the github repository to use the code.

:code:`git clone git@github.com:astrojoanna/mcdust.git`

First run
+++++++++
#. Make the :code:`/setups/default` setup 
    :code:`make SETUP_FILE=default`
#. Set the number of threads for OpenMP parallelisation.
    :code:`export OMP_NUM_THREADS=$OMP_NUM_THREADS`
#. Run the code
    | :code:`cd setups/Setup1`
    | :code:`./default setup.par`

Restart
+++++++

To (re)start the simulation from a specific snapshot. Here are the following steps.

#. Copy the snapshot of your choice to the setup directory where you want to restart
    :code:`cp your/data/directory/foo.h5 setups/restartrun/`
#. Rename the file to :code:`restart.h5`
    :code:`mv foo.h5 restart.h5`
#. Set the :code:`restart` condition in the :code:`setup.par` file to :code:`.true.`

