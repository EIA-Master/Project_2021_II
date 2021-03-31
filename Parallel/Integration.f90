! Edgar Alvarez Galera
! Eines Informàtiques Avançades: Project
! Last modification: 02/03/2021

	MODULE integration
	IMPLICIT NONE
	CONTAINS 
	

	SUBROUTINE Velocity_Verlet(N,dt,L,rcut,r,v,F,numproc,index_part,numsend,allgather,taskid,rnew,vnew,Fnew,pot)
	! ************************************************************************************************ !
	! This subroutine implements one step of the velocity Verlet algorithm.
	! INPUT:
	! - Sequential variables:
	!	N  --> Number of particles.
	! 	dt --> Time-step.
	!	L  --> Length of the lattice.
	!	r  --> Position of the particles.
	!	v  --> Velocity of the particles.
	!	F --> Forces of the particles (previous step)

	! - Parallelization variables:
	! 	index_part --> vector with beginning and ending particle to simulate
	! 	numsend --> vector containing the number of particles each processor has to take care of.
	! 	allgather --> vector with necessary information for MPI_ALLGATHERV
	! 	numproc --> Number of processors

	! OUTPUT:
	!	rnew --> New positions of the particles (after implementing the step).
	!	vnew --> New velocities of the particles.
	!	Fnew --> New forces between particles.

	! --------- !
	!  MODULES  !
	! --------- !
	use boundary
	use forces
	use parallel 
	! --------- !
	! ************************************************************************************************ !
	IMPLICIT NONE
	! Input:
	! - Sequential
	INTEGER, intent(in) :: N
	REAL*8, intent(in) :: dt, L, rcut
	! - Parallel
	INTEGER, intent(in) :: numproc,taskid
	INTEGER, intent(in) :: index_part(numproc,2),numsend(numproc),allgather(numproc)
	! Other variables:
	INTEGER i,j,k
	! Output:
	REAL*8, dimension(:,:) :: r,rnew,v,vnew,F,Fnew
	REAL*8 pot
	! ************************************************************************************************ !

	! Calculating new coordinates.
	rnew(:,index_part(taskid+1,1):index_part(taskid+1,2)) = & 
	r(:,index_part(taskid+1,1):index_part(taskid+1,2)) + &
	v(:,index_part(taskid+1,1):index_part(taskid+1,2))*dt + &
	.5d0*F(:,index_part(taskid+1,1):index_part(taskid+1,2))*dt*dt 

	! Sharing the information all processors have calculated.
	do k=1,3
    	
		call MPI_ALLGATHERV(rnew(k,index_part(taskid+1,1):index_part(taskid+1,2)),&
		numsend,MPI_DOUBLE_PRECISION, &
		rnew(k,index_part(taskid+1,1):index_part(taskid+1,2)), &
		numsend,allgather,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)

    enddo
	
	! Set periodic boundary conditions (put back particles that escape from the box).
	do i=1,N
	  do j=1,3
		call pbc1(L,rnew(j,i)) 
	  enddo
	enddo

	call force_LJ(N,L,rcut,rnew,numproc,index_part,taskid,Fnew,pot) ! New forces.	

	! New velocities.
	vnew(index_part(:,taskid+1,1):index_part(taskid+1,2)) = &
	v(index_part(:,taskid+1,1):index_part(taskid+1,2)) +&
	(F(index_part(:,taskid+1,1):index_part(taskid+1,2))+Fnew(index_part(:,taskid+1,1):index_part(taskid+1,2)))*.5d0*dt 

	return	
	END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	SUBROUTINE Integrate(Nsteps,Npart,Nradial,T,dt,rho,rcut,L, &
     	sigma,thermostat,r0,v0,numproc,index_part,numsend,allgather,taskid)
	! ************************************************************************************************ !
	! This subroutine implements the integration of the equations of motion for all the particles of the system.
  	! INPUT:
	! - Sequential variables:
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

	! - Parallelization variables:
	! 	index_part --> vector with beginning and ending particle to simulate
	! 	numsend --> vector containing the number of particles each processor has to take care of.
	! 	allgather --> vector with necessary information for MPI_ALLGATHERV
	! 	numproc --> Number of processors

  	! OUTPUT:
  	!	pf --> Final positions.
  	!	vf --> Final velocities.
  	!	ff --> Final forces.
  	
  	! Notice that the final arrays are returned as output so that the final configuration can be used as the initial one for a consecutive run.

	! --------- !
	!  MODULES  !
	use statistics
	use forces
	use radial distribution
	use parallel 
	! --------- !

	! ************************************************************************************************ !
	IMPLICIT NONE
	INTEGER, intent(in) :: Nsteps, Npart, Nradial
	REAL*8, intent(in) :: T, dt, rho, L, rcut, sigma
	LOGICAL, intent(in) :: thermostat
	REAL*8, dimension(:,:) :: r0, v0
	REAL*8, dimension(:,:) :: pf,vf,ff,pos,vel,forc,np,nv,nf
	INTEGER i,j,k
	REAL*8 time, KE, PE, totalE, Tinst, pressio, g(Nradial)
	! - Parallel
	INTEGER, intent(in) :: numproc,taskid
	INTEGER, intent(in) :: index_part(numproc,2),numsend(numproc),allgather(numproc)
	! ************************************************************************************************ !	
	
1  	FORMAT(A1,2X,3(F14.8,2X))
2  	FORMAT(6(F14.8,2X))	
	
	open(14,file="Positions.xyz")
	open(15,file="Thermodynamics.dat")
	open(16,file="gdr.dat")
	
	! Set initial state:
	time = 0d0	
	pos(:,index_part(taskid+1,1):index_part(taskid+1,2)) = r0(:,index_part(taskid+1,1):index_part(taskid+1,2))
	vel(:,index_part(taskid+1,1):index_part(taskid+1,2)) = v0(:,index_part(taskid+1,1):index_part(taskid+1,2))
	
	! Sharing the information all processors have calculated.
	do k=1,3
    	
		call MPI_ALLGATHERV(pos(k,index_part(taskid+1,1):index_part(taskid+1,2)),&
		numsend,MPI_DOUBLE_PRECISION, &
		pos(k,index_part(taskid+1,1):index_part(taskid+1,2)), &
		numsend,allgather,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)

    enddo

	! Calculating forces
	call force_LJ(Npart,L,rcut,pos,numproc,index_part,taskid,forc,PE) 	

	! Calculating pressure and kinetic energy
	call pressure(Npart,L,rho,pos,forc,Tinst,numproc,index_particles,taskid,pressio)
	call kinetic(Npart,vel,numproc,index_particles,taskid,KE)

	! Master processor tasks
	if (taskid == 0) then

		! Write initial configuration.
		write(14,*) Npart
		write(14,*) ""
		do j=1,Npart  
			write(14,*) "A", pos(:,j)
		enddo

		! Compute the Temperature and total energy
		call insttemp(Npart,KE,Tinst)
		totalE = totalenergy(PE,KE)	

		write(15,2) time, KE, PE, totalE, Tinst, pressio ! Write the values in "Thermodynamics.dat"	

	endif

	! Main loop.
	do i=1,Nsteps

		time = dble(i)*dt

		call Velocity_Verlet(Npart,dt,L,rcut,pos,vel,forc,numproc,index_part,numsend,allgather,taskid,np,nv,nf,PE)

		! If indicated, couple the system to an Andersen thermostat. Otherwise temperature can evolve with time.										  
		if (thermostat .eqv. .true.) call Andersen(T,index_particles,numproc,taskid,nv) 

		! Pressure and kinetic energy per step.
		call pressure(Npart,L,rho,np,nf,T,numproc,index_particles,taskid+1,pressio)
		call kinetic(Npart,nv,numproc,index_particles,taskid,KE)
		
		if (taskid == 0) then

			! Write the new positions to in XYZ format (trajectories):	
			if (mod(i,1000) .eq. 0) then	
				write(14,*) Npart  
				write(14,*) "" 
				do j=1,Npart
					write(14,1) "A", pos(:,j)
				enddo
			endif

			! Compute the statistics:
			call insttemp(Npart,KE,Tinst)
			totalE = totalenergy(PE,KE)
			call radial_dist(Npart,Nradial,L,np,g)

			if (mod(i,500).eq.0) write(15,2) time, KE, PE, totalE, Tinst, pressio 

		endif		

		! For the next iteration, update positions, velocities and forces:
		pos = np
		vel = nv
		forc = nf
		
	enddo	
	
	! Normalize all the radial distribution with master processor.
	if (taskid == 0) then
	
		call radial_dist_norm(Nradial,Npart,Nsteps,L,rho,g)

		do i=1,Nradial
			write(16,*) sigma*(dble(i)-.5d0)*L/(2d0*dble(Nradial)),g(i)
		enddo
	
	endif
	
	close(14)
	close(15)
	close(16)

	return
	END SUBROUTINE Integrate
	
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	SUBROUTINE ANDERSEN(T,index_particles,numproc,taskid,v)
	! Andersen thermostat with probability 10% of interacting with the bath.
	! INPUT:
	! - Sequential
	!	T --> Temperature of the bath.

	! - Parallelization variables:
	! 	index_part --> vector with beginning and ending particle to simulate
	! 	numproc --> Number of processors

	! OUTPUT:
	!	v --> Velocities.
      IMPLICIT NONE
      INTEGER N, i, k
      REAL*8 T,sigma,NU,v1,
	  REAL*8, dimension(:,:):: v
      PARAMETER(nu = 0.1d0)
	  INTEGER, intent(in) :: numproc,taskid
	  INTEGER, intent(in) :: index_part(numproc,2)
      
      SIGMA = DSQRT(T)

	call srand(12345678)

      do i=index_part(taskid+1,1):index_part(taskid+1,2)
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
