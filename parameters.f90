! Author: Jaume Garcia
! Mòdul que conté tots els paràmetres necessaris per fer la simulació i
! la conversió. 
module parameters
implicit none

double precision, parameter :: T = 27.52d0 ! Heat bath temperature (reduced units)
double precision, parameter :: rho = 0.8d0 ! Density in reduced units
double precision, parameter :: dt = 0.00001d0 ! Time step
double precision, parameter :: sigma = 2.64d0 ! Sigma of the gas (Angstroms)
double precision, parameter :: eps = 10.9d0     ! Epsilon of the gas (K)
double precision, parameter :: mass = 4.0026d0  ! Mass (g/mol)
integer, parameter :: Natoms = 125 ! Number of atoms of the system
integer, parameter :: Nsteps = 50000 ! Number of steps of our simulation

end module parameters
