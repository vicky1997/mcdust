Example: Post-processing data from hydro simulations
====================

:code:`mcdust` can be used to post-process results from hydrodynamic simulations. Currently, :code:`mcdust` can post-process steady-state 2D gas velocity data and perform sust growth simulations.
The post-processing version of mcdust is in the github branch :code:`postprocess` and to access it

.. code-block:: bash

    git checkout postprocess

Preparing the data to post-process
++++++++++++++++++++++++++++++++++
To be able to post-process velocity data the following steps are required

- Convert the data to cylindrical coordinates if it is originally not in cylindrical coordinate.
- Make your data scale-free: vertical height in the units of scale height :math:`H_g`, velocities in the unit of sound speed :math:`c_s`
- structure the data in a regular 2D grid v(r,z) (interpolation might be required)
- save the data in the following format separately for radial and vertical velocities (or depending on which one you want to implement)

.. code-block:: bash

    rsize = np.size(r)
    zsize = np.size(z)
    rzsize = [rsize, zsize]
    filename = 'radialvelocity.inp'
    with open(filename,"ab") as f:
        np.savetxt(f,np.c_[int(rsize)])
    with open(filename,"ab") as f:
        f.write(b"#r,z(AU)\n")
    with open(filename,"ab") as f:
        np.savetxt(f,np.c_[r,z])
    with open(filename,"ab") as f:
        f.write(b"#rr(AU),v(cs)\n")
    with open(filename,"ab") as f:
        np.savetxt(f,np.c_[rr.flatten(),vr.flatten()]) 

and update the file name in the setup.par and turn on the respective flag in the preprocs.opt file. It is :code:`#READGASVELOCITIES_R` for reading radial velocities and :code:`#READGASVELOCITIES_Z` for vertical velocities. An example of the setup can be found in :code:`/setups/gasvel/`