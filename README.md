
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
- Forces/Radial distribution: Alba
- Statistics: Ignasi
- Integration: Edgar
- Visualitzation: Jaume

---

# 1. Repository information

This repository contains four main folders, all of the containing the corresponding code, data and images:

- Sequential

- Parallel_v1

- Parallel_v2

- Parallel_v3

## Modules shared by Sequential and Parallel:

### reader.f90:

Module with subroutine that reads the input of the program.

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

Makefile used to:
1. Compile all the modules together along with the main program.

       :> make 

2. Plot all the figures.

       :> make plot

3. Clean the .0 and .mod files.

       :> make clean
       
4. Get help.

       :> make help
       
       
## data.txt

This file contains all the information that the program needs to simulate the system. It is very important not to change the format of it and to put the numbers matching Fortran 90's syntax (i.e. 2.d0).

## Modules only in Parallel

### parallel.f90

This module has the subroutine that assigns the number of particles to simulate each processor.

### openmpi.sub

This is the script used to send works to the cluster's queue.

---

# 2. How to compile?

1. First of all download all the files in this repository into your local computer, in the same directory.
2. Go to this directory via command line and write:

        :> make

This will create and executable named programa.exe.

# 3. Running the program:
## Sequential

3. To run the program execute programa.exe.

  3.1. If you are using LINUX:

    :> ./programa.exe

  3.2. If you are using Windows:

    :> .\programa.exe
    
 ## Parallel
 
 3. To run the program execute programa.exe.

        :> mpirun -np 4 programa.exe
    
    where this 4 is the number of processors you want to execute the program.
    
 # 4. Results
 
   This will generate two files:

   3.2.1. positions.xyz: To visualize the dynamics of the system with a molecular dynamics visualizer like VMD.
  
   3.2.2. Thermodynamics.dat: Contains the evolution of all the relevant magnitudes throughout the simulation.

   3.2.3. gdr.dat: Contains the information to plot the radial distribution.


4. To visualize the plots you must type in your command line

        :> make plots

this will generate three plots:

  4.1. Energy plot containing the evolution of the potential, the kinetic and the total energy. 
  
  4.2. Evolution of the intantaneous temperature.
   
  4.3. Evolution of the instantaneous pressure.
