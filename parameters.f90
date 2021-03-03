! Author: Jaume Garcia
! Mòdul que conté tots els paràmetres necessaris per fer la simulació i
! la conversió. 
module parameters
contains
implicit none
double precision T,rho,dt,sigma,eps,mass
integer Natoms,Nsteps
logical true 


rho = 0.8d0 ! Density in reduced units
Natoms = 125 ! Number of atoms of the system

Nsteps = 50000 ! Number of steps of our simulation
dt = 0.00001d0 ! Time step

sigma = 2.64d0 ! Sigma of the gas (Angstroms)
eps = 10.9d0     ! Epsilon of the gas (K)
mass = 4.0026d0  ! Mass (g/mol)

T = 27.52d0        ! Heat bath temperature (reduced units)

end module parameters
