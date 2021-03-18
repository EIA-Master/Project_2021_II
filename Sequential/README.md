
# Project_2021_II
---


# Group members:

-Ignasi Puch (admin: ignasipuch @ github)

-Edgar Alvarez

-Alba Fischer

-Jaume Garcia

-Arnau Prat


## Tasks:

- Initial state: Arnau
- Boundy conditions: Jaume
- Forces: Alba
- Statistics: Ignasi
- Integration: Edgar
- Visualitzation of results: Arnau

---

# Repository information

This repository contains nine files consisting on:

- Seven modules
- One main program
- One file to visualize the results with gnuplot
- One makefile to join all the parts

## Modules

### Parameters.f90:

Module with all the needed parameters known of the problem.

### Initialize.f90:

Contains one subroutine that creates a simple cubic lattice of particles with a random velocity associated to a Maxwell-Boltzmann distribution at the given temperature. 

### boundary.f90:

Contains one subroutine that computes the position of a particle considering the minimum image convention of the periodic boundary conditions.

### Forces.f90:

Contains one subroutine that with a given position of the particles calculates the Lennard-Jones forces associated to them and its potential, considering only peair-wise interactions.

### statistics.f90:

Contains three subroutines and one function:
- kinetic: Calculates the kinetic energy of a given set of velocities.
- insttemp: Given a kinetic energy calculates the instant temperature of the system.
- pressure: With the position of the particles and the forces applied to each calculates the instant pressure of the system.

The totalenergy function calculates the sumation of the potential energy and the kinetic energy.

### Integration.f90:

This module contains four subroutines:
- Velocity_Verlet: Which is a single step of the velocity Verlet integrator.
- Integrate: Main subroutine of the project which makes the integration and calculation of all the interesting variables.
- Andersen: Subroutine used to implement the Andersen thermostat in the system.
- Boxmuller: Soubroutine to extract a random number from a gaussian normal distribution.

### radial_distribution.f90:
This module contains two subroutines:
- radial_dist: Calculates the radial distribution function
- radial_dist_norm: Normalize radial distribution function

## Main

The main code uses all the previously mentioned modules, directly or indirectly, to calculate all the magnitudes of interest of a Helium gas at 300K.

## Visualization

This part only contains one file that plots the evolution of the magnitudes through time.

## Makefile

Makefile used to compile all the modules together along with the main program.

---

# How to compile?

1. First of all download all the files in this repository into your local computer, in the same directory.
2. Go to this directory via command line and write:

        :> make

This will create and executable named programa.exe.

3. To run the program execute programa.exe.

  3.1. If you are using LINUX:

    :> ./programa.exe

  3.2. If you are using Windows:

    :> .\programa.exe

   This will generate two files:

   3.2.1. positions.xyz: To visualize the dynamics of the system with a molecular dynamics visualizer like VMD.
  
   3.2.2. Thermodynamics.dat: Contains the evolution of all the relevant magnitudes throughout the simulation.

4. To visualize the plots you must type in your command line

        :> gnuplot grafics.gnu

this will generate three plots:

  4.1. Energy plot containing the evolution of the potential, the kinetic and the total energy. 
  
  4.2. Evolution of the intantaneous temperature.
   
  4.3. Evolution of the instantaneous pressure.

# Files from the simulation

In this repository we have also included thermodynamics.dat file with all the calculations performed.




>>>>>>> 8f94f0d17e3c1cc55fdcee4fcb3ab385ccfe5ba7
