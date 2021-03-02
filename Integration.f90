! Edgar Alvarez Galera
! Eines Informàtiques Avançades: Project
! Last modification: 02/03/2021

	MODULE integration
	IMPLICIT NONE
	CONTAINS 
	
	SUBROUTINE Velocity_Verlet(N,dt,L,rcut,r,v,F,rnew,vnew,Fnew)
! This subroutine implements one step of the velocity Verlet algorithm.
! INPUT
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
	USE boundary, forces
	IMPLICIT NONE
	INTEGER, intent(in) :: N
	REAL*8, intent(in) dt, L, rcut
	REAL*8, intent(inout) :: r(3,N), rnew(3,N)
	REAL*8, intent(inout) :: v(3,N), vnew(3,N)
	REAL*8, intent(inout) :: F(3,N), Fnew(3,N)
	INTEGER i, j
	REAL*8, intent(out) :: pot
	


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

C------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	SUBROUTINE Integrate(Nsteps,Npart,T,dt,rho,rcut,L,thermostat
     &	p0,v0,pf,vf,ff)
	use statistics, forces
	IMPLICIT NONE
	INTEGER, intent(in) :: Nsteps, Npart
	REAL*8, intent(in) :: T, dt, rho, L, rcut
	LOGICAL, intent(in) :: thermostat
	REAL*8, intent(in) :: r0(3,Npart), v0(3,Npart)
	REAL*8, intent(out) :: pf(3,Npart), vf(3,Npart), ff(3,Npart)
	REAL*8 pos(3,Npart), v(3,Npart), forc(3,Npart)
	REAL*8 newp(3,Npart), newv(3,Npart), newf(3,Npart)
	INTEGER i, j
	REAL*8 time, KE, PE, totalE, Tinst, pressure
	
	
1  	FORMAT(A1,2X,3(F14.8,2X))
2  	FORMAT(5(F14.8,2X))	
	
	open(14,file="Positions.xyz")
	open(15,file="Thermodynamics.dat")
	
! Set initial state:
	time = 0d0
	
	pos(:,:) = r0(:,:)
	vel(:,:) = v0(:,:)
	call force_LJ(Npart,L,rcut,pos,forc,PE)	

	write(14,*) Npart
	write(14,*) ""
	do j=1,Npart
		write(14,*) "A", pos(j,:)
	enddo
	
	call kinetic(Npart,nv,KE)
	call insttemp(Npart,KE,T)
	totalE = totalenergy(PE,KE)
	call pressure(Npart,L,rho,pos,forc,T,pressure)
	
	write(15,2) time, KE, PE, totalE, Tinst, pressure !Write the values in "thermodynamics.dat"		
	
	
	do i=1,Nsteps
		time = dble(i)*dt
		call Velocity_Verlet(Npart,dt,L,pos,vel,forc,np,nv,nf)
		if (thermostat .eqv. .true.) call Andersen(Npart,T,nv)
		
! Write the new positions to in XYZ format (trajectories):		
		write(14,*) Npart  ! Number of particles in simulation
		write(14,*) "" ! Blank line
		do j=1,Npart
			write(14,1) "A", pos(j,:)
		enddo
! For the next iteration:
		pos = np
		vel = nv
		forc = nf
		
! Compute the statistics:
		call kinetic(Npart,nv,KE)
		call insttemp(Npart,KE,T)
		totalE = totalenergy(PE,KE)
		call pressure(Npart,L,rho,pos,forc,T,pressure)

		write(15,2) time, KE, PE, totalE, Tinst, pressure !Write the values in "thermodynamics.dat"		
		
	enddo	
	
	close(14)
	close(15)

! Set final arrays (outputs).	
	pf = pos
	vf = vel
	ff = forc
	return
	END SUBROUTINE Integrate
	
C------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	SUBROUTINE ANDERSEN(N,T,v)
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
