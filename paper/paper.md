
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

The Solar System and other planetary systems formed in a disk of gas and dust around the young forming central star called a protoplanetary disk. The overall gas dynamics, the interaction between the gas and dust and the collisions between the dust particles leave an imprint on the formation pathways of the planets and therefore, on life. Recent observations with telescopes such as ALMA [@Andrews_2018] show a diversity of dust structures in protoplanetary disks and it is imperative that we have models that are able to describe the observations and more generally, to help move forward towards the goal of an end-to-end planet formation theory. We present `mcdust`, a parallel simulation code aimed at modelling the evolution of dust in protoplanetary disks to better understand the dust size distributions and dust collisional evolution.

# Statement of Need

Modelling dust coagulation and evolution is an essential part of understanding protoplanetary disks and planet formation. The evolution of micrometer-sized dust grains to millimeter-sized aggregates sets the tone for planetesimal formation. Dust evolution is inextricably linked with the formation of substructures in disks. Dust dynamics is also of prime importance when understanding protoplanetary disk chemistry. Therefore it is very important to understand the processes involved in dust growth and dynamics and how they influence different aspects of disk dynamics and planet growth in detail.

Dust coagulation is generally modelled with the Smoluchowski equation [@Smoluchowski1916], an integro-differential equation. The drawback with the method is that it is difficult to track histories of dust particles. Furthermore, adding a property to track adds a dimension to solving the Smoluchowski Equation [@Stammler2017].

Alternatively one can model dust coagulation by performing Monte Carlo simulations. Monte Carlo methods are well suited for stochastic processes and that makes it a suitable method for dust coagulation [@Gillespie1975]. For such a method, the representative particle approach [@Zsom2008] is very suitable and the representative particle approach is Lagrangian in nature, meaning we track the particles and hence their histories. We can add properties to the particles that can be tracked without much computational complexity. This advantage is very useful when we want to combine dust evolution and chemistry, where we want to look at the thermal histories of particles and their compositions when chemical models are included. 

Although there have been several publications that have made use of the Monte Carlo method to model dust coagulation in protoplanetary disks [@Krijt2018;@Houge2023], there is no available open source code that models dust coagulation using Monte Carlo methods and such a code would be very helpful as an alternative way for modelling dust coagulation.

# State of the field

Several codes exist to model dust growth in protoplanetary disks and we focus on discussing the open sources codes in the field. There are semi-analytical methods like `twopoppy` [@Birnstiel2012] and `TriPoD` [@Pfeil2024] that model the overall dust size distributions without modelling full dust coagulation. Some examples of open source codes that model full dust coagulation are `dustpy` [@Stammler2022], a 1D code to simulate gas and dust evolution and `cuDisc` [@Robinson2024], a 2D (radial-vertical) code modelling dust coagulation and disk evolution. Both of these codes solve the Smoluchowski equation are Eulerian in nature, i.e. they follow volume rather than mass. `PHANTOM` [@Price], is an open source Smoothed Particle Hydrodynamics (SPH) code that models models dust growth using a monodisperse growth or a 'single-size' approximation and do not solve the Smoluchowski equation for dust growth [@Vericel2021]. 

`mcdust` models dust coagulation and transport in the vertical and radial directions. The currently included collisional outcomes are dust growth by sticking, fragmentation of dust particles and erosion, where a small particle chips a portion of the large particle. We employ a representative particle approach detailed in @Zsom2008 to track a limited number of particles instead of tracking every particle, saving computational time. We have a static power-law gas disk with temperature assumed to be vertically isothermal. Dust coagulation depends on the local gas properties and therefore we bin particles into grids in order to perform collisions. We make use of an adaptive grid approach where we make sure that each cell has equal number of representative particles. This guarantees that there are always sufficient particles to resolve the physics of collisions. Figure \ref{fig:grid} shows a sketch of our adaptive grid model and the different physical processes simulated by `mcdust`. `mcdust` in its first iteration was introduced in @Drazkowska2013. In its current version, the code has been optimized and modified to make it faster and adaptable to suit the needs of the user. But the physics of the code has largely remains the same and we refer the reader to @Drazkowska2013 for the details.

The particle based approach to dust growth has the advantage of being able to track particle histories and add properties/composition that can be tracked with little computational overhead making `mcdust` useful for performing dust coagulation simulations with different properties (e.g.porosity) and compositions for the dust. Including dust growth/dynamics in hydrodynamic simulations of protoplanetary disks can be very expensive and computationally complex. This is where `mcdust` can also be used to post-process data from hydrodynamic simulations to understand the dynamics and evolution of dust in different conditions without much computational complexity.  

![(Left) An overview of the different physical processes in `mcdust` and (right) a representative sketch  of the adaptive grid method used in `mcdust` to group particles and perform collisions. \label{fig:grid}](Figures/gridnew.png)

# Software Design

The representative particle approach can be computationally intensive and therefore time consuming to run a simulation with enough resolution [@Drazkowska2014]. `mcdust` has been designed with the intention to overcome this aspect and is written in `Fortran` and parallelised with `OpenMP` which enables the user to utilize high performance computing systems on the node level. The modules of `mcdust` are written in such a way that the user can modify the specific module as per their requirements. For e.g, the major processes in the simulations can be turned off/on with the preprocessor directives before running the simulation . Large parts of the code, like the protoplanetary disk setup (`discstruct.F90`) and the initialisation of the dust particles across the disk (`initproblem.F90`) are written in a way that it can be adapted by the user for their specific needs. The data I/O is done with the HDF5 data format to ensure faster writing and reading of data for further analysis. We also provide a python script that aids the user with the data visualisation which can be modified to the user's needs. 

# Research Impact Statement

`mcdust` has already been used in its previous iterations to model the evolution of dust in protoplanetary disks. @Drazkowska2013 used to model the pile-up of dust in the inner-edge of the dead zones in protoplanetary disks. We refer the reader to @Drazkowska2014 and @Drazkowska2018 for other instances where `mcdust` has been used to model dust growth in protoplantary disks. Recently, @Vaikundaraman2025 investigated the depletion of refractory carbon in the inner Solar System using `mcdust` demonstrating its capabilities of adding properties to particles and tracking them with little computational overhead.

## Benchmark
We compare a run of our code with a run from the open source 1D dust coagulation code `dustpy` [@Stammler2022]. The parameters used for the simulations are listed in the Table \ref{tab:params}. The gas evolution in both the codes was switched off and the codes were run with a static gas background to exclude any differences in gas treatment that might influence the outcome of dust evolution.

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

We ran a simulation for 10000 years. Figure \ref{fig:benchmark} shows the dust surface density $\sigma_d$ as a function of particle mass and distance from star at the end of the simulation for both `mcdust` and `dustpy`. 

![A comparison of the dust surface density $\sigma_d$ as a function of particle mass and distance from star between `mcdust` (left) and `dustpy` (right).\label{fig:benchmark}](Figures/dustpyvsmcdust.png){width=75%}


It is evident from Figure \ref{fig:benchmark} that both the codes have similar overall outcomes but they do have certain differences. The most striking one is that `mcdust` does not provide coverage of the regions of parameter space that do not include sufficiently high fraction of the dust mass. This holds true for all Monte Carlo based codes. But with higher resolution simulation this issue can be overcome. The other important factor is that `mcdust` does not face the issue of artificially sped-up growth that `dustpy` and other Smoluchowski equation based codes tend to encounter. And the 2D r-z structure of `mcdust` also helps us to investigate processes like sedimentation driven coagulation  [@Drazkowska2013] that are not usually seen in 1D simulations like `dustpy`. This can be seen in image at around 50 AU where `mcdust` has larger surface densities higher masses when compared to `dustpy`. For a more detailed discussion of the differences between the two approaches to dust growth, i.e., the Monte Carlo method and the Smoluchowski equation approach we point the reader to @Drazkowska2014.

# AI usage disclosure
No AI tools including generative AI tools have been used in the development of the software, creation of the documentation or in the writing of the paper. 

# Acknowledgments
The authors acknowledge funding from the European Union under the European Union’s Horizon Europe Research \& Innovation Programme 101040037 (PLANETOIDS). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.

# References
