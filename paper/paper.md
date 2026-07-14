
---
title: 'mcdust: A 2D Monte Carlo code for dust coagulation in protoplanetary disks'
tags:
  - dust evolution
  - protoplanetary disks
  - FORTRAN
  - planet formation
authors:
  - name: Vignesh Vaikundaraman
    orcid: 0000-0002-2451-9574
    affiliation: 1
  - name: Nerea Gurrutxaga
    orcid: 0009-0008-3256-9564
    affiliation: 1
  - name: Joanna Drążkowska
    orcid: 0000-0002-9128-0305
    affiliation: 1
affiliations:
 - name: Max Planck Institute for Solar System Research, Justus-von-Liebig-Weg 3
         37077 Gottingen
   index: 1
date: 11 July 2025
bibliography: paper.bib
---

# Summary

The Solar System and other planetary systems are formed in a disk of gas and dust around the young, forming central star, called a protoplanetary disk. The micron-sized particles collide and grow through several processes to form the planetary architectures as we see them today. The overall gas dynamics, the interaction between the gas and dust, and the collisions between the dust particles leave an imprint on the formation pathways of the planets and, therefore, on life. Recent observations with telescopes such as ALMA [@Andrews_2018] show a diversity of dust structures in protoplanetary disks, and it is imperative that we have models that are able to describe the observations and, more generally, to help move forward towards the goal of an end-to-end planet formation theory. We present `mcdust`, a parallel Monte Carlo code aimed at modelling the evolution of dust in protoplanetary disks to better understand the dust size distributions and dust collisional evolution.

# Statement of Need

Modelling dust coagulation and evolution is an essential part of understanding protoplanetary disks and planet formation. The evolution of micrometer-sized dust grains to millimeter-sized aggregates sets the tone for planetesimal formation. Dust evolution is inextricably linked with the formation of substructures in disks [@Pinilla2012] and is of prime importance when understanding disk chemistry [@Woitke2016;@Houge2025]. Therefore, it is important to understand the processes involved in dust growth and dynamics and how they influence different aspects of disk dynamics and planet growth in detail.

Dust coagulation is generally modelled with the Smoluchowski equation [@Smoluchowski1916], an integro-differential equation. The drawback of the methods is that when we add a property to track (e.g, charge, composition, porosity), it adds a dimension to solving the Smoluchowski equation and necessitates a parameterized approach when considering dust properties to reduce the complexity in these methods.[@Stammler2017;@Blum2006;@Okuzumi2009].

One can also model dust coagulation by performing Monte Carlo simulations. Monte Carlo methods are well-suited for stochastic processes, making them suitable for dust coagulation [@Gillespie1975]. For such a method, the representative particle approach [@Zsom2008;@Ormel2007] is very suitable, allowing us to track the particles and hence their histories. We can add properties to the particles that can be tracked without much computational complexity. This advantage is very useful when we want to combine dust evolution and chemistry, allowing us to examine the thermal histories of particles and their compositions when chemical models are included [@Floresrivera2025]. 

Although there have been several publications that have made use of the Monte Carlo method to model dust coagulation in protoplanetary disks [@Krijt2018;@Houge2023], there is no available open-source code that models dust coagulation using Monte Carlo methods, and such a code would be very helpful.

# State of the field

Several codes exist to model dust growth in protoplanetary disks, and we focus on discussing the open-source codes in the field. There are semi-analytical methods like `twopoppy` [@Birnstiel2012] and `TriPoD` [@Pfeil2024] that model the overall dust size distributions without modelling full dust coagulation. Some examples of open-source codes that model full dust coagulation are `dustpy` [@Stammler2022], a 1D code to simulate gas and dust evolution, and `cuDisc` [@Robinson2024], a 2D (radial-vertical) code modelling dust coagulation and disk evolution. Both of these codes solving the Smoluchowski equation are Eulerian in nature, i.e., they follow volume rather than mass. `PHANTOM` [@Price] is a Smoothed Particle Hydrodynamics code that models dust growth using a monodisperse growth approximation and does not solve the Smoluchowski equation for dust growth [@Vericel2021]. We refer the reader to @Eriksson2026 for a comparison of the different dust evolution codes, including `mcdust`.

`mcdust` models dust coagulation and transport in the vertical and radial directions. The currently included collisional outcomes are dust growth by sticking, fragmentation of dust particles, and erosion, where a small particle chips a portion of the large particle. We employ a representative particle approach detailed in @Zsom2008 to track a limited number of particles instead of tracking every particle, saving computational time. We have a static power-law gas disk in the default version where the gas surface density is $\Sigma_g=\Sigma_{AU}(r/\mathrm{AU})^{-1}$ with $\Sigma_{AU}$ being the gas surface density at 1 AU. The temperature is assumed to be vertically isothermal, given by $T=T_{AU}(r/\mathrm{AU})^{-0.5}$ with $T_{AU}$ being the temperature at 1 AU. The static power law gas disk is a good assumption for simulations less than Myr timescales. This can be modified as per the user's needs. Dust coagulation depends on the local gas properties, and therefore, we bin particles into grids in order to perform collisions. We make use of an adaptive grid approach where we make sure that each cell has an equal number of representative particles. This guarantees that there are always sufficient particles to resolve the physics of collisions. Figure \ref{fig:grid} shows a sketch of our adaptive grid model and the different physical processes simulated by `mcdust`. `mcdust` in its first iteration was introduced in @Drazkowska2013. In its current version, the code has been optimized and modified to make it faster and adaptable to suit the needs of the user. But the physics of the code has largely remained the same, and we refer the reader to @Drazkowska2013 for the details.

The particle-based approach makes `mcdust` useful for performing dust coagulation simulations with different properties (e.g., porosity) and compositions for dust with little computational overhead. Including dust growth/dynamics in hydrodynamic simulations of protoplanetary disks can be very expensive and computationally complex [@Drazkowska2019]. This is where `mcdust` can also be used to post-process data from hydrodynamic simulations to understand the dynamics and evolution of dust in different conditions. The current version of `mcdust` can take 2D (r-z) steady-state gas velocities as input and perform dust growth simulations for different gas velocity structures seen in hydrodynamic simulations.  

![(Left) An overview of the different physical processes in `mcdust` and (right) a representative sketch  of the adaptive grid method used in `mcdust` to group particles and perform collisions. \label{fig:grid}](Figures/gridnew.png)

# Software Design

The representative particle approach can be computationally intensive and, therefore, time-consuming when running a simulation with enough resolution [@Drazkowska2014]. `mcdust` has been designed with the intention to overcome this aspect and is written in `Fortran` and parallelised with `OpenMP`, which enables the user to utilize high-performance computing systems on the node level. The modules of `mcdust` are written in such a way that the user can modify the specific module as per their requirements. For example, the major processes in the simulations can be turned off/on with the preprocessor directives before running the simulation. Large parts of the code, like the protoplanetary disk setup (`discstruct.F90`) and the initialisation of the dust particles across the disk (`initproblem.F90`), are written in a way that they can be adapted by the user for their specific needs. The data I/O is done with the HDF5 data format to ensure faster writing and reading of data for further analysis. We also provide a Python script that aids the user with the data visualisation, which can be modified to the user's needs. 

# Research Impact Statement

`mcdust` has already been used in its previous iterations to model the evolution of dust in protoplanetary disks. @Drazkowska2013 used it to model the pile-up of dust in the inner edge of the dead zones in protoplanetary disks. We refer the reader to @Drazkowska2014 and @Drazkowska2018 for other instances where `mcdust` has been used to model dust growth in protoplanetary disks. @Vaikundaraman2025 investigated the depletion of refractory carbon in the inner Solar System, and @Gurrutxaga2026a modelled the formation of carbonaceous chondrites in a pressure bump using `mcdust`, demonstrating its capabilities of adding properties to particles and tracking them. @Gurrutxaga2026b showed a more robust way to track particle properties using a modified version of `mcdust`, which will soon be implemented in `mcdust`. @Vaikundaraman2026 used `mcdust` to post-process magnetohydrodynamic simulations to understand dust evolution in complex gas flows.

## Benchmark
We compare a run of our code with a run from the open-source 1D dust coagulation code `dustpy` [@Stammler2022]. The parameters used for the simulations are listed in Table \ref{tab:params}. The gas evolution in both codes was switched off, and the codes were run with a static gas background to exclude any differences in gas treatment that might influence the outcome of dust evolution.

| Parameter                                   | Value                 |
|:------------------------------------------- | :--------------------:|
| gas surface density at 1 AU                 | 1000 g/$\mathrm{cm}^2$|
| temperature at 1 AU                         | 280 K                 |
| turbulence strength ($\alpha$)              | $10^{-3}$             |
| fragmentation velocity ($v_{\mathrm{frag}}$)| 10 m/s                |
| erosion mass ratio                          | 10                    |
| gas surface density power law $p$           | 1                     |
| temperature power law $q$                   | 0.5                   |

: Parameters used for the benchmark simulations to compare `mcdust` and `dustpy`. \label{tab:params}

We ran a simulation for 10,000 years. Figure \ref{fig:benchmark} shows the dust surface density $\sigma_d$ as a function of particle mass and distance from the star at the end of the simulation for both `mcdust` and `dustpy`. 

![A comparison of the dust surface density $\sigma_d$ as a function of particle mass and distance from the star between `mcdust` (left) and `dustpy` (right).\label{fig:benchmark}](Figures/mcdustvsdustpy.png){width=100%}


It is evident from Figure \ref{fig:benchmark} that both codes have similar overall outcomes, but there are certain differences. The most striking one is that `mcdust` does not provide coverage of the regions of parameter space that do not include a sufficiently high fraction of the dust mass. This holds for all Monte Carlo-based codes. But with higher resolution simulations, this issue can be overcome. The other important factor is that `mcdust` does not face the issue of artificially sped-up growth [@Drazkowska2014] that `dustpy` and other Smoluchowski equation-based codes tend to encounter [@Stammler2022]. And the 2D r-z structure of `mcdust` also helps us to investigate processes like sedimentation-driven coagulation [@Drazkowska2013] that are not usually seen in 1D simulations like `dustpy`. This can be seen in the image at around 50 AU, where `mcdust` has larger surface densities and higher masses when compared to `dustpy`. For a more detailed discussion of the differences between the Monte Carlo method and the Smoluchowski equation approach, we point the reader to @Drazkowska2014. Benchmarking the simulations against analytical coagulation kernels is an important way to validate the coagulation models, and we refer the reader to @Drazkowska2013 and the `mcdust` documentation for the different coagulation kernel tests.

# AI usage disclosure
No AI tools, including generative AI tools, have been used in the development of the software, creation of the documentation, or writing of the paper. 

# Acknowledgments
The authors acknowledge funding from the European Union under the European Union’s Horizon Europe Research \& Innovation Programme 101040037 (PLANETOIDS). Views and opinions expressed are, however, those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.

# References
