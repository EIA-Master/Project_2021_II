! Author: Jaume Garcia
! Programa principal que crida les subrutines per inicialitzar el sistema
! i realitza la simulació amb un bany tèrmic. 
program main
use Initialize
use Integration
use reader 
use parallel

implicit none 

real*8 L,rcut
real*8 t_ini,t_final
integer master

LOGICAL true ! Paràmetre lògic que afegeix el termostat

REAL*8, allocatable :: pos0(:,:)
REAL*8, allocatable :: vel0(:,:)
REAL*8, allocatable :: pf(:,:)
REAL*8, allocatable :: vf(:,:)
REAL*8, allocatable :: ff(:,:)
REAL*8 temp,density,timestep,sigma,epsilon,mass
integer Natoms,Nsteps,Nradial
integer numproc


call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
t_ini = MPI_WTIME()
master = 0

! Lectura del fitxer amb les dades necessàries per la simulació
call read_files(temp,density,timestep,sigma,epsilon,mass,Natoms,Nsteps,Nradial)

call assignacio(Natoms,numproc)

allocate(pos0(3,Natoms))
allocate(vel0(3,Natoms))
allocate(pf(3,Natoms))
allocate(vf(3,Natoms))
allocate(ff(3,Natoms))

! Inicilitzem el sistema calculant les posicions i les velocitats inicials
call initial(Natoms,density,temp,pos0,vel0,L)

rcut= L/2.d0 ! Cut-off

! Compartim les posicions i velocitats inicials amb tots els processadors
call MPI_BCAST(pos0, Natoms*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
call MPI_BCAST(vel0, Natoms*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

! Integrem les equacions del moviment i escrivim els resultats en fitxers
call Integrate(Nsteps,Natoms,Nradial,temp,timestep,density,rcut,L,sigma,true,pos0, &
vel0,pf,vf,ff)

deallocate(pos0)
deallocate(vel0)
deallocate(pf)
deallocate(vf)
deallocate(ff)

t_final = MPI_WTIME()
IF (taskid==0) then
    print*, "Número de processadors = ",numproc
    print*, "Temps de simulació = ",endtime-starttime
END IF
call MPI_FINALIZE(IERROR)
return

end program main
