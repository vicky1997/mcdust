Input and Output files
============
Input
++++++
There are two required input files to run a simulation with mcdust, setup.par and prepocs.opt. The details are each file are explained below.

Input file 1: :code:`setup.par`
-----------------------

This is the file to specify the simulation parameters. The parameters to be entered in the file are written below

| :code:`number_of_particles_per_cell`  - Number of particles in a cell          
| :code:`number_of_radial_zones`  -  Number of radial grids            
| :code:`number_of_vertical_zones` -  Number of vertical grids                 
| :code:`steps_between_outputs`  - Number of iterations between outputs                   
| :code:`time_between_outputs_[yrs]` - time duration between outputs               
| :code:`maximum_time_of_simulation_[yrs]`  - simulation time in yrs
| :code:`minimum_radius_[AU]` - inner radial boundary of the simulation in AU                      
| :code:`maximum_radius_[AU]` - outer radial boundary of the simulation in AU                      
| :code:`monomer_radius_[cm]` - radius of the particle monomers in cm                      
| :code:`material_density_[g/cm3]` - internal density of the particles in g/cm3                 
| :code:`dmmax`  -  mass optimization value for grouped collisions                                 
| :code:`evaporation_radius_[AU]` - radius inside which particles will not be tracked anymore                  
| :code:`dust_to_gas_ratio` - the dust to gas ratio in the disk                        
| :code:`fragmentation_velocity_[cm/s]` - fragmentation velocity of the particles            
| :code:`alpha` - the alpha turbulence value to specifiy the strength of the turbulence in the disk                                   
| :code:`sigma_gas_[g/cm2]` - the gas surface density at 1 AU                        
| :code:`temperature_[K]`  - the temperature in K at 1 AU                         
| :code:`erosion_mass_ratio` - mass ratio condition to trigger erosion                     
| :code:`data_directory` - directory name to write the data  
| :code:`restart` - to be set to true if we are starting a simulation from a snapshot


If the values for the parameters are not entered in the file, the default values will be taken as specified in the :code:`parameters.f90` file. New parameters can also be added to the simulation via the :code:`parameters.f90` file.


Input file 2: :code:`preprocs.opt`
--------------------------
This file specifies the regions of the code to be compiled based on compiler options. The current compiler options are listed

| :code:`VERTICALMODEL`: activate this to setup a 1D vertical model
| :code:`TRANSPORT`: activate dust transport
| :code:`VERTSETTLING`: activate vertical settling of particles
| :code:`RADRIFT`: activate radial drift of particles
| :code:`COLLISIONS`: activate dust coagulation
| :code:`EROSION`: activate erosion as a collision outcome
| :code:`AUXDATA`: write timestep data

Comment on uncomment the options in the file based on the setup of the run.


Output
++++++

The output files are of the hdf5 file format with the extension :code:`.h5`. Each snapshot has its own file and can accessed in :code:`outputs/` directory.

Each data file contains the properties of each representative particle (called as a swarm). The properties that are stored in the default version are given below.

| ID number of swarm
| mass of the representative particle in grams
| cylindrical radius in AU (distance from the central star)
| height above midplane in AU
| Stokes Number of the representative particle
| Radial velocity v_r of the representative particle in cm/s
| Vertical velocity v_z of the representative particle in cm/s

Each data file also contains the mass of a swarm (in g) and the output number of the file as metadata.

The metadata of the code can be accessed using the :code:`h5dump` command as shown in the example below.

.. code-block:: bash
    
    h5dump -H data/swarms-00000.h5

The output of the command shows us the metadata of the file :code:`swarms-00000.h5`

.. code-block:: bash

    HDF5 "swarms-00000.h5" {
    GROUP "/" {
    ATTRIBUTE "author" {
        DATATYPE  H5T_STRING {
            STRSIZE 20;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
    }
    ATTRIBUTE "code" {
        DATATYPE  H5T_STRING {
            STRSIZE 9;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
    }
    ATTRIBUTE "mass_of_swarm" {
        DATATYPE  H5T_IEEE_F32LE
        DATASPACE  SCALAR
    }
    ATTRIBUTE "output_number" {
        DATATYPE  H5T_IEEE_F32LE
        DATASPACE  SCALAR
    }
    GROUP "swarms" {
        DATASET "swarmsout" {
            DATATYPE  H5T_COMPOUND {
                H5T_STD_I32LE "id_number";
                H5T_IEEE_F64LE "mass_of_a_particle_[g]";
                H5T_IEEE_F64LE "cylindrical_radius_[AU]";
                H5T_IEEE_F64LE "height_above_midplane_[AU]";
                H5T_IEEE_F64LE "Stokes_number";
                H5T_IEEE_F64LE "Radial_velocity_v_r_[cm/s]";
                H5T_IEEE_F64LE "Vertical_velocity_v_z_[cm/s]";
            }
            DATASPACE  SIMPLE { ( 1, 1048576 ) / ( H5S_UNLIMITED, 1048576 ) }
        }
    }
    GROUP "times" {
        DATASET "timesout" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
        }
    }
    }

To add new properties to be written, one can add the same in the :code:`hdf5output.f90` file. 
