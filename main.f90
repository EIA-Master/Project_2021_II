program main
use Initialize
use Integration

implicit none 
INTEGER Nsteps, Natoms
REAL*8 T, dt, rho, L, rcut
LOGICAL true
REAL*8, allocatable :: pos0(:,:)
REAL*8, allocatable :: vel0(:,:)
REAL*8, allocatable :: pf(:,:)
REAL*8, allocatable :: vf(:,:)
REAL*8, allocatable :: ff(:,:)

T=300.d0
rho=0.8d0
Natoms=125
Nsteps=50000
dt=0.00001

allocate(pos0(3,Natoms))
allocate(vel0(3,Natoms))
allocate(pf(3,Natoms))
allocate(vf(3,Natoms))
allocate(ff(3,Natoms))

call initial(Natoms,rho,T,pos0,vel0,L)

rcut= L/2.d0

call Integrate(Nsteps,Natoms,T,dt,rho,rcut,L,true,pos0,vel0,pf,vf,ff)

deallocate(pos0)
deallocate(vel0)
deallocate(pf)
deallocate(vf)
deallocate(ff)
return
end program main
