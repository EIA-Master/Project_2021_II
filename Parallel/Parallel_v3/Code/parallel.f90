! Author: Edgar Alvarez 
module parallel
!     Module that contains the subroutine to assign the corresponding
!     particles to each processor.
    IMPLICIT NONE

    contains
    SUBROUTINE assignacio(Natoms,numproc,index_particles,allgather_pointer,num_send,pairs,pair_limits)

!     Assigns the corresponding particles to each processor.
!     INPUT: Natoms (number of particles), numproc (Number of processors)
!     OUTPUT: index_particles (array that contains in (i,1) component, 
!     lower limit sent to processor and in (i,2) upper limit), 
!     allgather_pointer(punter d'assignació), num_send (number of particles
!     sent to each processor), pairs (number of particle pairs),
!     pair_limits (number of particle pairs in each processor)

        implicit none
        integer :: i,j,k
        integer :: particle_per_proc,residu
        integer :: Natoms, numproc
        integer index_particles(numproc,2)
        integer allgather_pointer(numproc)
        integer num_send(numproc)
        
        integer :: npairs
        integer :: pairs(Natoms*(Natoms-1)/2 ,2)
        integer :: pair_per_proc
        integer :: pair_limits(numproc,2)
       
        particle_per_proc = int(real(Natoms)/real(numproc))
        residu = mod(Natoms,numproc)

!	Nombre de partícules que li corresponen a cada processador
            DO i=1,numproc
              IF((i-1)<residu)then
                index_particles(i,1)=(particle_per_proc+1)*(i-1)+1
                index_particles(i,2)=index_particles(i,1)+particle_per_proc
              else
              index_particles(i,1)=(i-1)*particle_per_proc+residu+1
              index_particles(i,2)=index_particles(i,1)+particle_per_proc-1
              end if
            END DO
!	Punter d'assignació per la funció MPI_ALLGATHERV
            DO i=1,numproc
              if(i.eq.1)then
                allgather_pointer(i)=0
                else
                allgather_pointer(i)=index_particles(i-1,2)
                end if
                num_send(i)=index_particles(i,2)-index_particles(i,1)+1
            end do  

!	Parelles de partícules            
            k = 0
		DO i=1,Natoms-1
	  	  DO j=i+1,Natoms
			k = k + 1
			pairs(k,1) = i 
			pairs(k,2) = j
                        write(1964,*) pairs(k,1), pairs(k,2)
	  	  ENDDO
		ENDDO  
              write(1964,*)
            write(1964,*)  
!	Nombre de parelles a cada processador
                pair_per_proc = Natoms*(Natoms-1) / (2*numproc)
		residu = mod(Natoms*(Natoms-1)/2,numproc)
		DO i=1,numproc
		!pair_limits(i,1) = pair_limits(i-1,2) + 1
	  	!pair_limits(i,2) = pair_limits(i,1) + nint(real(i*npairs)/real(numproc)) - 1
		  if ((i-1) .lt. residu) then
	  		pair_limits(i,1) = (pair_per_proc+1)*(i-1) + 1
	  		pair_limits(i,2) = pair_limits(i,1) + pair_per_proc
	  	  else
	  		pair_limits(i,1) = (pair_per_proc)*(i-1) + 1 + residu
	  		pair_limits(i,2) = pair_limits(i,1) + pair_per_proc - 1	
	  	  endif
                write(1964,*) pair_limits(i,1), pair_limits(i,2)        
	  	ENDDO 
END Subroutine
END module parallel
