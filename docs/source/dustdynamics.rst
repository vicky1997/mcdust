Dust dynamics
=============

The transport of dust particles is determined by the disk properties (radial drift, settling, turbulence).
We briefly explain the different components of the dust velocites in this section. For further detail we refer the reader to `Drążkowska et al 2013 <https://www.aanda.org/articles/aa/abs/2013/08/aa21566-13/aa21566-13.html>`_ .

Radial velocity
+++++++++++++++

The radial advection velocity is given by,

.. math:: 

    v^r_d = \frac{v_g}{1 + \mathrm{St}^2} + \frac{2v_\eta}{\mathrm{St} + \mathrm{St}^{-1}}

The above equation takes into account radial drift of dust particles due to the differences in dust and gas radial and azimuthal velocities.

he radial diffusion is implemented as a random kick on the particle distributions as in `Ciesla 2010 <https://iopscience.iop.org/article/10.1088/0004-637X/723/1/514/meta>`_ and `Zsom et al 2011 <https://www.aanda.org/articles/aa/abs/2011/10/aa16515-11/aa16515-11.html>`_ .
The total diffusion velocity is given by,

.. math:: 

    v_{\mathrm{diff}} = \frac{\Delta r}{\Delta t} + D_d\frac{1}{\rho_g}\frac{\partial \rho_g}{\partial r}

where :math:`\Delta r` is the turbulent kick obtained from a Gaussian distribution. The second term accounts for the non-homogeniety of the gas distribution.

Vertical velocity
+++++++++++++++++
In the vertical direction, settling of the dust particles due to the gravity of the central star is of importance. This is given by,

.. math:: 

    v_d^z = -z\Omega_K \mathrm{St}

The diffusion in the vertical direction is implemented in the same way as in the radial direction as turbulent kicks (replacing :math:`\Delta r` with :math:`\Delta z`).

Advection timestep
++++++++++++++++++

The timestep is limited by both advection and collisions in order to resolve the processes properly. 

To properly resolve advection we limit the timestep according to the Courant condition given by,

.. math:: 

    \Delta t^x = \frac{\Delta x_{\mathrm{min}}}{v^x_{\mathrm{max}}},

where 

This is calculated for both the radial and vertical directions and the minimum is chosen :math:`\Delta t^d = \mathrm{min}(\Delta t^r, \Delta t^z)`.

We do not want more collisions than the number particles in a cell per timestep, this limit is implemented by

.. math:: 

    \Delta t^c = \frac{n_{\mathrm{cell}}}{n_{\mathrm{max}}^{\mathrm{coll}}} \Delta t_{\mathrm{old}},

where :math:`n_{\mathrm{cell}}` is the number of particles in a cell and :math:`n_{\mathrm{max}}^{\mathrm{coll}}` is the maximum number of collisions from the previous timestep. The final timestep is chosen as a minimum between the advection and collision limits :math:`\Delta t = \mathrm{min}(\Delta t^d, \Delta t^c)`