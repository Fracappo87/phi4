# phi4
Simulation program for generating field configuration for a 3D scalar lattice theory with quartic interaction, using a Markov Chain Monte Carlo algorithms.
It also includes routines for implementing a gradient flow evolution, using a 4th order Runge Kutta.
The flow evolution is a simple Gaussian heat equation, according to the work published in https://inspirehep.net/record/1505605.
It has been created for numerically studying a specific set of lattice Ward Identities built with probes evolved along the flow. 
From these identities, the renormalization constants of the lattice renormalized Energy Momentum Tensor can be numerically extracted.

The program is divided in several modules, all linked together to several main files, to be used for different purposes.
Here is the list of the following modules

----------------------------------------------------------------------------------------------------------------------------

1) Geometry

This module is necessary for setting up the geometry of the system (D-dimensional lattice) and to generate next-neighbour     vectors to be used for navigate through the lattice.

----------------------------------------------------------------------------------------------------------------------------

2) Io
  
This modules provides functions for 
    a) reading data from input files and check their correctness.
    b) reading a given field configuration to be used as a starting point, and check its correctness.
    c) saving field configurations to be used later by the user.
    d) writing useful informations about the execution of the program, as well as the value of computed observables, to a specific log file.
    
 ----------------------------------------------------------------------------------------------------------------------------
 
3) Observables
 
This module provides different subroutines used to compute observables built from fundamental scalar fields and fields evolved according to the gradient flow.
 ----------------------------------------------------------------------------------------------------------------------------
 
4) Random
 
Module for the generation of pseudo random numbers: it uses the  RANLUX random number generator, written by Martin Luscher http://luscher.web.cern.ch/luscher/ranlux/ 
 
 ----------------------------------------------------------------------------------------------------------------------------

5) Update
  
This is the core module: it contains routines for the generation of field configuration using a mixture of Metropolis and cluster algorithm (Swenden-Wang formulation)
 It contains routines for the implementation of the simple Metropolis update and the creation of clusters using the Ising model embedded into the scalar theory.
 
----------------------------------------------------------------------------------------------------------------------------

6) WilsonFlow

Module for the evolution of field configurations along the flow: it uses a simple gaussian smearing as prescription for evolging lattice scalar fields.
It also implement the construction of the Jacobian map from boundary fields to field configurations at positive flow time.


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TESTING

A directory Test has also been included: it contains a set of codes to be used for testing different modules of the program:

--------------------------------------------------------------------------------------------------------------------------

1) main_corr.c    

This program computes field correlators.

--------------------------------------------------------------------------------------------------------------------------

2) main_testJACOB.c, testJACOB.c

This program check the correcteness of 4th order Runge Kutta integrator for the evolution of the Jacobian map. It computes the accumulated error due to time discretization by computing the difference between numerical and exact solution of the heat equation, as a function of time steps.

--------------------------------------------------------------------------------------------------------------------------

3) main_testRK.c, testRK.c

This program check the correcteness of 4th order Runge Kutta integrator for the evolution of the scalar field. It computes the accumulated error due to time discretization by computing the difference between numerical and exact solution of the heat equation, as a function of time steps.

--------------------------------------------------------------------------------------------------------------------------

4) main_testMASS.c   

This program is used for computing the renormalized mass of the theory (useful for check code correcteness making comparison with available measurements)


