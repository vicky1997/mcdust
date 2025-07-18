``mcdust`` Documentation
=======================

| ``mcdust`` is a code to model the growth and evolution of dust in protoplanetary disks.

| ``mcdust`` is a 2D (r-z) code that simulates the vertical and radial evolution and growth of dust. The code makes use of the representative particle approach to simulate the dust particles and uses Monte Carlo methods to simulate dust growth.

| The details are ``mcdust`` are presented in  `Drążkowska et al 2013 <https://www.aanda.org/articles/aa/abs/2013/08/aa21566-13/aa21566-13.html>`_ & Vaikundaraman et al. (in prep). The code units are in cgs and we use a cylindrical coordinate system (r-z) in the code.

| Here we outline the basic setup and usage of the code and some ways in which one can make use of the features of the code.

.. toctree::
    :maxdepth: 2
    :caption: Contents
    
    setup
    inputoutput
    codesetup
    dustdynamics
    dustcoagulation
    notebooks/mcdust_example.ipynb
    example_1
    example_2
    citation
    
