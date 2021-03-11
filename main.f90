! Author: Jaume Garcia
! Programa principal que crida les subrutines per inicialitzar el sistema
! i realitza la simulació amb un bany tèrmic. 
program main
use Initialize
use Integration
use parameters 

implicit none 
real*8 L,rcut

LOGICAL true ! Paràmetre lògic que afegeix el termostat

REAL*8, allocatable :: pos0(:,:)
REAL*8, allocatable :: vel0(:,:)
REAL*8, allocatable :: pf(:,:)
REAL*8, allocatable :: vf(:,:)
REAL*8, allocatable :: ff(:,:)


allocate(pos0(3,Natoms))
allocate(vel0(3,Natoms))
allocate(pf(3,Natoms))
allocate(vf(3,Natoms))
allocate(ff(3,Natoms))

! Inicilitzem el sistema calculant les posicions i les velocitats inicials
call initial(Natoms,rho,T,pos0,vel0,L)

rcut= L/2.d0 ! Cut-off

! Integrem les equacions del moviment i escrivim els resultats en fitxers
call Integrate(Nsteps,Natoms,Nradial,T,dt,rho,rcut,L,sigma,true,pos0, &
vel0,pf,vf,ff)

deallocate(pos0)
deallocate(vel0)
deallocate(pf)
deallocate(vf)
deallocate(ff)
return

end program main
