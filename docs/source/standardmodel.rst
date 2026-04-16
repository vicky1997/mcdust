Standard run
============

Here we describe the standard disk model present in :code:`discstruct.f90` and the parameter choices in the default run setup in :code:`setups/default`.

Disk model
++++++++++

The default gas structure in :code:`discstruct.f90` is static and does not evolve and we use a power law disk model. The gas surface density is given by

.. math:: 

    \Sigma_\mathrm{g}(r) = \Sigma_\mathrm{g0} \left(\frac{r}{AU}\right)^{-1}

where :math:`\Sigma_\mathrm{g0}` is the gas surface density at 1 AU. This is defined in the function :code:`sigmag(x,time)` which can be modified to suit the user's needs.

The 2D distribution of gas is given by,

.. math::

    \rho_\mathrm{g}(r,z) = \frac{\Sigma_g(r)}{\sqrt{2\pi}H_g} \mathrm{exp}\left(\frac{-z^2}{2H_g^2}\right)

where :math:`H-g` is the gas scale height. This is defined in the function :code:`densg(x,z, time)`

We use a vertically isothermal disk temperature model where the sound speed given by,

.. math::

    cs(r) = cs_0\left(\frac{r}{AU}\right)^{-0.25}

where :math:`cs_0 = \sqrt{k_BT_0/m_{H_2}}` is the sound speed at 1 AU. This is defined in the function :code:`cs(x)`.

The gas pressure is then given by,

.. math::

    P(r,z) = \rho_g c_s^2.

This is defined in the function :code:`Pg(x,z,time)`

The turbulent strength :math:`alpha` is defined in the function :code:`alpha(x)` and is constant over time and space in the default model. The gas velocity is set to zero in the default model and can be changed in the function :code:`vgas(x,time)`. 

Default run
+++++++++++

Listed below are the parameters for the default simulation. We encourage the user to run the default setup first as a test run to ensure proper working of the code.

| :code:`number_of_particles_per_cell`  - 256 (recommended to have atleast 200 particles a cell)          
| :code:`number_of_radial_zones`  -  32           
| :code:`number_of_vertical_zones` -  16                
| :code:`steps_between_outputs`  - 10000                   
| :code:`time_between_outputs_[yrs]` - 250            
| :code:`maximum_time_of_simulation_[yrs]`  - 10000
| :code:`minimum_radius_[AU]` - 1                    
| :code:`maximum_radius_[AU]` - 100                     
| :code:`monomer_radius_[cm]` - 1e-4                  
| :code:`material_density_[g/cm3]` - 1.2                 
| :code:`dmmax`  -  0.001                                 
| :code:`evaporation_radius_[AU]` - 0.99                  
| :code:`dust_to_gas_ratio` - 0.01                     
| :code:`fragmentation_velocity_[cm/s]` - 1000            
| :code:`alpha` - 1.e-3                                   
| :code:`sigma_gas_[g/cm2]` - 1000                       
| :code:`temperature_[K]`  - 280                      
| :code:`erosion_mass_ratio` - 10                   
| :code:`data_directory` - 'data'

