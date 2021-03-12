! Author: Jaume Garcia
! Programa principal que crida les subrutines per inicialitzar el sistema
! i realitza la simulació amb un bany tèrmic. 
program main
use Initialize
use Integration
use parameters 

implicit none 
<<<<<<< HEAD
real*8 L,rcut
=======
INTEGER Nsteps, Natoms
REAL*8 T, dt, rho, L, rcut
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c

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
<<<<<<< HEAD
call Integrate(Nsteps,Natoms,Nradial,T,dt,rho,rcut,L,sigma,true,pos0, &
vel0,pf,vf,ff)
=======
call Integrate(Nsteps,Natoms,T,dt,rho,rcut,L,true,pos0,vel0,pf,vf,ff)
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c

deallocate(pos0)
deallocate(vel0)
deallocate(pf)
deallocate(vf)
deallocate(ff)
return

end program main
