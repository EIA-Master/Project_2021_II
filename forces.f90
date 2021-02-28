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

      subroutine force_LJ(Natoms,L,cutoff,r,F,Upot)
      implicit none
      double precision r(3,Natoms),F(Natoms,3),L
      double precision dr(3),cutoff, cf6,cf12,Upot
      double precision d2,d4,d6,d8,d12,d14,modf     
      integer i,j,k,Natoms
      
      F = 0.d0
      Upot = 0.d0
      do i = 1,Natoms
            do j = 1,Natoms
                  dr(:) = r(i,:)-r(j,:)
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
                        F(i,:) = F(i,:) + modf*dr(:)
                        F(j,:) = F(j,:) - modf*dr(:)
                        Upot = Upot+ 4.d0*(1.d0/D12 - 1.d0/D6) - &
                        4.d0*(1.d0/cf12 - 1.d0/cf6)
                  endif
            enddo
      enddo
      return
      end 
      
end module forces
