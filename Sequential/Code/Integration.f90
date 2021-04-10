! Edgar Alvarez Galera
! Eines Informàtiques Avançades: Project
! Last modification: 02/03/2021

	MODULE integration
	IMPLICIT NONE
	CONTAINS 
	

	SUBROUTINE Velocity_Verlet(N,dt,L,rcut,r,v,F,rnew,vnew,Fnew,pot)
	! This subroutine implements one step of the velocity Verlet algorithm.
	! INPUT:
	!	N  --> Number of particles.
	! 	dt --> Time-step.
	!	L  --> Length of the lattice.
	!	r  --> Position of the particles.
	!	v  --> Velocity of the particles.
	!	F --> Forces of the particles (previus step)

	! OUTPUT:
	!	rnew --> New positions of the particles (after implementing the step).
	!	vnew --> New velocities of the particles.
	!	Fnew --> New forces between particles.
	use boundary
	use forces
	IMPLICIT NONE
	INTEGER, intent(in) :: N
	REAL*8, intent(in) :: dt, L, rcut
	REAL*8, intent(inout) :: r(3,N), rnew(3,N)
	REAL*8, intent(inout) :: v(3,N), vnew(3,N)
	REAL*8, intent(inout) :: F(3,N), Fnew(3,N)
	INTEGER i, j
	REAL*8 pot

	rnew(:,:) = r(:,:) + v(:,:)*dt + .5d0*F(:,:)*dt*dt ! New coordinates.
	
	do i=1,N
	  do j=1,3
		call pbc1(L,rnew(j,i)) ! Set periodic boundary conditions (put back particles that escape from the box).
	  enddo
	enddo
	
	call force_LJ(N,L,rcut,rnew,Fnew,pot) ! New forces.	
	vnew(:,:) = v(:,:) + (F(:,:)+Fnew(:,:))*.5d0*dt ! New velocities.

	return	
	END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	SUBROUTINE Integrate(Nsteps,Npart,Nradial,T,dt,rho,rcut,L, &
     	sigma,thermostat,r0,v0,pf,vf,ff)
	! This subroutine implements the integration of the equations of motion for all the particles of the system.
  	! INPUT:
  	!	Nsteps --> Total number of steps in which the integration is implemented.
  	!	Npart --> Number of particles of the system.
  	!	Nradial --> Number of points for the radial distribution function.
  	!	T --> Temperature.
  	!	dt --> Time-step.
  	!	rho --> Density of particles.
  	!	rcut --> Cutoff radius.
  	!	L --> Linear length of the system.
  	!	thermostat --> Logical variable. If is true, an Andersen thermostat will be implemented to keep temperature roughly constant (around the temperature of the bath T). 
  	!	r0 --> Initial positions of the particles.
  	!	v0 --> Initial velocities of the particles.
  	! OUTPUT:
  	!	pf --> Final positions.
  	!	vf --> Final velocities.
  	!	ff --> Final forces.
  	
  	! Notice that the final arrays are returned as output so that the final configuration can be used as the initial one for a consecutive run.
  	
	use statistics
	use forces
	use radial_distribution
	IMPLICIT NONE
	INTEGER, intent(in) :: Nsteps, Npart, Nradial
	REAL*8, intent(in) :: T, dt, rho, L, rcut, sigma
	LOGICAL, intent(in) :: thermostat
	REAL*8, intent(in) :: r0(3,Npart), v0(3,Npart)
	REAL*8, intent(out) :: pf(3,Npart), vf(3,Npart), ff(3,Npart)
	REAL*8 pos(3,Npart), vel(3,Npart), forc(3,Npart)
	REAL*8 np(3,Npart), nv(3,Npart), nf(3,Npart)
	INTEGER i, j
	REAL*8 time, KE, PE, totalE, Tinst, pressio, g(Nradial)
	
	
1  	FORMAT(A1,2X,3(F14.8,2X))
2  	FORMAT(6(F14.8,2X))	
	
	open(14,file="Positions.xyz")
	open(15,file="Thermodynamics.dat")
	open(16,file="gdr.dat")
	
	! Set initial state:
	time = 0d0	
	pos(:,:) = r0(:,:)
	vel(:,:) = v0(:,:)
	call force_LJ(Npart,L,rcut,pos,forc,PE)	

	write(14,*) Npart
	write(14,*) ""
	do j=1,Npart  
		write(14,*) "A", pos(:,j)
	enddo

	! Compute the initial thermodynamic quantities of the system:	
	call kinetic(Npart,vel,KE)
	call insttemp(Npart,KE,Tinst)
	totalE = totalenergy(PE,KE)
	call pressure(Npart,L,rho,pos,forc,Tinst,pressio)
	
	write(15,2) time, KE, PE, totalE, Tinst, pressio !Write the values in "Thermodynamics.dat"		
	
	
	do i=1,Nsteps
		time = dble(i)*dt

		call Velocity_Verlet(Npart,dt,L,rcut,pos,vel,forc,np,nv,nf,PE)

		if (thermostat .eqv. .true.) call Andersen(Npart,T,nv) ! If indicated, couple the system to an Andersen thermostat. Otherwise temperature can evolve with time.
		
		! Write the new positions to in XYZ format (trajectories):	
		if (mod(i,500) .eq. 0) then	
			write(14,*) Npart  
			write(14,*) "" 
			do j=1,Npart
				write(14,1) "A", pos(:,j)
			enddo
		endif
		! For the next iteration, update positions, velocities and forces:
		pos = np
		vel = nv
		forc = nf
		
		! Compute the statistics:
		call kinetic(Npart,nv,KE)
		call insttemp(Npart,KE,Tinst)
		totalE = totalenergy(PE,KE)
		call pressure(Npart,L,rho,pos,forc,T,pressio)
		call radial_dist(Npart,Nradial,L,pos,g)

		if (mod(i,500).eq.0) write(15,2) time, KE, PE, totalE, Tinst, pressio 		
		
	enddo	
	
	call radial_dist_norm(Nradial,Npart,Nsteps,L,rho,g)
	
	
	do i=1,Nradial
		write(16,*) sigma*(dble(i)-.5d0)*L/(2d0*dble(Nradial)),g(i)
	enddo
	
	close(14)
	close(15)
	close(16)

	! Set final arrays (outputs).	
	pf = pos
	vf = vel
	ff = forc
	return
	END SUBROUTINE Integrate
	
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	SUBROUTINE ANDERSEN(N,T,v)
	! Andersen thermostat with probability 10% of interacting with the bath.
	! INPUT:
	!	N --> Number of particles.
	!	T --> Temperature of the bath.
	! OUTPUT:
	!	v --> Velocities.
      IMPLICIT NONE
      INTEGER N, i, k
      REAL*8 T, sigma, NU, v(3,N), v1, v2
      PARAMETER(nu = 0.1d0)
      
      SIGMA = DSQRT(T)

	call srand(12345678)

      do i=1, N
      if (rand().lt.NU) then
      do k=1,3      
            CALL BOXMULLER(SIGMA, dble(rand()), dble(rand()), v1, v2)
            v(i,k) = v1
      enddo
      endif
      enddo
      
      return
      END SUBROUTINE
  
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	SUBROUTINE BOXMULLER(SIGMA,X1,X2,XOUT1,XOUT2)
	! This subroutine converts two random numbers (x1 and x2) from a uniform distribution U(0,1) into two others (xout1 and xout2) compatible with a Gaussian one of standard deviation "sigma".
      IMPLICIT NONE
      REAL*8 PI, sigma, x1, x2, xout1, xout2
      
      PI = 4d0*datan(1d0)
      
      XOUT1=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*PI*x2)
      XOUT2=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dsin(2d0*PI*x2)
       
      END SUBROUTINE BOXMULLER
     
     END MODULE integration
