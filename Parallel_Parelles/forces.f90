module forces
      use boundary
      use initialize
      implicit none
      contains

!     Calculates the energy and force of a system of N identical 
!     particles interacting through a Lennard-Jones potential, 
!     considering periodic boundary conditions.
!     INPUT: r (array with x, y, z coordinates of particles),Natoms
!     (number of particles), cutoff (threshold),L (box length)
!     OUTPUT: F (force) and Upot(potential energy)

      subroutine force_LJ(Natoms,L,cutoff,r,numproc,pair_limits,taskid,F,Upot)
      implicit none
      double precision, dimension(:,:):: r,F
      double precision L
      double precision dr(3),cutoff, cf6,cf12,Upot
      double precision d2,d4,d6,d8,d12,d14,modf     
      integer i,j,k,p,Natoms
! Variables paralel      
      integer numproc,taskid
!      integer :: index_particles(numproc,2)
      integer :: MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror
!      integer :: npairs
      integer :: pairs(Natoms*(Natoms-1)/2, 2)
!      integer :: pair_per_proc
      integer :: pair_limits(numproc, 2)
!	COMMON/pairs/npairs,pairs,pair_per_proc,pair_limits

      
      F = 0.d0
      Upot = 0.d0
      
      do p = pair_limits(taskid+1,1), pair_limits(taskid+1,2)
      	i = pairs(p,1) ! Particle 1 of pair p.
      	j = pairs(p,2) ! Particle 2 of pair p.
!      do i = index_particles(taskid+1,1), index_particles(taskid+1,2)
!            do j = i+1,Natoms
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
                   F(:,j) = F(:,j) - modf*dr(:)
                   Upot = Upot + 4.d0*(1.d0/d12 - 1.d0/d6) - 4.d0*(1.d0/cf12 - 1.d0/cf6)
                  endif
!            enddo
      enddo
!	  Compartim la força de cada partícula per tots els processadors     
      call MPI_ALLREDUCE(F,F,Natoms*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
! 	  Sumem totes les contribucions de l'energia potencial de tots els processadors      
      call MPI_REDUCE(Upot,Upot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      return
      end subroutine force_LJ
      
end module forces
