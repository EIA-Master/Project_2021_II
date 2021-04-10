! Author: Jaume Garcia
! Programa principal que crida les subrutines per inicialitzar el sistema
! i realitza la simulació amb un bany tèrmic. 
program main
use Initialize
use Integration
use reader 
use parallel

implicit none
include 'mpif.h' 

real*8 L,rcut
real*8 t_ini,t_final
integer master
integer:: taskid,comm,ierror,reslen,message,partner,request
integer:: stat(MPI_STATUS_SIZE)

LOGICAL true ! Paràmetre lògic que afegeix el termostat

REAL*8 temp,density,timestep,sigma,epsilon,mass  
integer Natoms,Nsteps,Nradial
REAL*8, allocatable :: pos0(:,:)
REAL*8, allocatable :: vel0(:,:)
!Parallel variables
integer numproc
integer, allocatable  :: index_particles(:,:)
integer, allocatable ::  allgather_pointer(:)
integer, allocatable :: num_send(:)
integer ii
! Lectura del fitxer amb les dades necessàries per la simulació
call read_files(temp,density,timestep,sigma,epsilon,mass,Natoms,Nsteps,Nradial)

call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
if (taskid.eq.0) t_ini = MPI_WTIME()
master = 0

allocate(index_particles(numproc,2))
allocate(allgather_pointer(numproc))
allocate(num_send(numproc))

call assignacio(Natoms,numproc,index_particles,allgather_pointer,num_send)

allocate(pos0(3,Natoms))
allocate(vel0(3,Natoms))

! Inicilitzem el sistema calculant les posicions i les velocitats inicials
call initial(Natoms,density,temp,pos0,vel0,L)
rcut= L/2.d0 ! Cut-off

! Compartim les posicions i velocitats inicials amb tots els processadors
call MPI_BCAST(pos0, Natoms*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
call MPI_BCAST(vel0, Natoms*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

! Integrem les equacions del moviment i escrivim els resultats en fitxers
call Integrate(Nsteps,Natoms,Nradial,temp,timestep,density,rcut,L,sigma,true,pos0, &
vel0,numproc,index_particles,num_send,allgather_pointer,taskid)

deallocate(pos0)
deallocate(vel0)

deallocate(index_particles)
deallocate(allgather_pointer)
deallocate(num_send)

IF (taskid==0) then
    t_final = MPI_WTIME()
    print*, "Número de processadors = ",numproc
    print*, "Temps de simulació = ",t_final-t_ini
END IF
call MPI_FINALIZE(ierror)

end program main
