! Author: Alba Fischer Carles
module radial_distribution      
      use boundary
      implicit none
      contains
      
!     Calculates the radial distribution function of a system of N  
!     identical particles interacting through a Lennard-Jones potential, 
!     considering periodic boundary conditions.
!     INPUT: r (array with x, y, z coordinates of particles),Natoms
!     (number of particles), Nradial(number of bins),L (box length)
!     OUTPUT: g (radial distribution normalized)
      subroutine radial_dist(Natoms,Nradial,L,r,g)
      implicit none
      double precision g(Nradial),r(3,Natoms),dr(3),dist,bin_size,L
      integer Natoms,i,j,k,coef,Nradial
      
      bin_size =  L/(2.d0*dble(Nradial))
      do i = 1,Natoms-1
            do j = i+1,Natoms
                  dr(:) = r(:,i)-r(:,j)
                  do k = 1,3
                        call PBC1(L,dr(k))
                  enddo
                  dist = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                  coef = int(dist/bin_size) + 1
                  if (coef .lt. (Nradial + 1)) then
                        g(coef) = g(coef) + 1.d0
                  endif
            enddo
      enddo

      return
      end subroutine radial_dist

      subroutine radial_dist_norm(Nradial,Natoms,Nsteps,L,rho,g)
      implicit none
      double precision g(Nradial),V(Nradial),bin_size,rho,pi,L
      integer Nradial,i,Natoms,Nsteps
      pi = atan(1.d0)*4 
      bin_size =  L/(2.d0*dble(Nradial))
      
      ! Volume of every bin

      do i  = 1,Nradial
            V(i) = (4.d0/3.d0)*pi*(dble(i)*bin_size)**3-(4.d0/3.d0)*pi &
            *(dble(i-1)*bin_size)**3
      enddo

      ! Normalize radial distribution g
      g = g*(2.d0/(V*rho*Natoms*Nsteps))

      return
      end subroutine radial_dist_norm
end module radial_distribution
