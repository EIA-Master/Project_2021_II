! Author: Alba Fischer

module forces
! Module that contains the subroutine to implement forces and 
! compute Lenard Jones potential. Double work parallel version
use boundary
use initialize
implicit none
contains

! Calculates the energy and force of a system of N identical 
! particles interacting through a Lennard-Jones potential, 
! considering periodic boundary conditions.
! INPUT: r (array with x, y, z coordinates of particles),Natoms
! (number of particles), cutoff (threshold),L (box length)
! OUTPUT: F (force) and Upot(potential energy)

subroutine force_LJ(Natoms,L,cutoff,r,numproc,index_particles,F,Upot)
implicit none
double precision, dimension(:,:) :: r,F
double precision L,Far,Uar
double precision dr(3),cutoff, cf6,cf12,Upot
double precision d2,d4,d6,d8,d12,d14,modf     
integer i,j,k,Natoms

! Variables paralel      
integer numproc,taskid
integer :: index_particles(numproc,2)
integer :: ierror
      
include 'mpif.h'

call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)

F = 0.d0
Upot = 0.d0

do i = index_particles(taskid+1,1), index_particles(taskid+1,2)
      do j = 1, Natoms
        if (i .ne. j) then
            dr(:) = r(:,i)-r(:,j)

            do k = 1,3
                  call PBC1(L,dr(k))
            enddo

            d2 = dr(1)**2 + dr(2)**2 + dr(3)**2

            if (d2.lt.(cutoff*cutoff)) then
                  d4 = d2*d2
                  d6 = d4*d2
                  d8 = d6*d2
                  d12 = d6*d6
                  d14 = d12*d2
                  cf6 = cutoff**6
                  cf12 = cf6*cf6
                  modf = (48.d0/d14 - 24.d0/d8)
                  F(:,i) = F(:,i) + modf*dr(:)
!                  F(:,j) = F(:,j) - modf*dr(:)
                  Upot = Upot+ 4.d0*(1.d0/d12 - 1.d0/d6) - &
                  4.d0*(1.d0/cf12 - 1.d0/cf6)
            endif
        endif
      enddo
enddo

Upot = Upot / 2d0 !Descomptem el doble comptatge.

! Compartim la força de cada partícula per tots els processadors     
do i=1,Natoms
    do j=1,3
        call MPI_ALLREDUCE(F(j,i),Far,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        F(j,i)=Far
    enddo
enddo

! Sumem totes les contribucions de l'energia potencial de tots els processadors      
call MPI_ALLREDUCE(Upot,Uar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
Upot=Uar
return
end subroutine force_LJ
      
end module forces
