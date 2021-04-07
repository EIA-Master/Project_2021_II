module parallel
    IMPLICIT NONE

    contains
    SUBROUTINE assignacio(Natoms,numproc,index_particles,allgather_pointer,num_send)
        implicit none
        integer :: i, j, k
        integer :: particle_per_proc,residu
        integer :: Natoms, numproc
        integer :: index_particles(numproc,2)
        integer :: allgather_pointer(numproc)
        integer :: num_send(numproc)
        
        integer :: npairs
        integer :: pairs(Natoms*(Natoms-1)/2 ,2)
        integer :: pair_per_proc
        integer :: pair_limits(numproc,2)
        
!        COMMON/pairs/npairs,pairs,pair_per_proc,pair_limits
        
        particle_per_proc = int(real(Natoms)/real(numproc))
        residu = mod(Natoms,numproc)
        
!        allocate(index_particles(numproc,2),allgather_pointer(numproc),num_send(numproc))
        

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
      
!      npairs = Natoms*(Natoms-1)/2  
!	allocate(pairs(npairs,2), pair_limits(numproc,2))
!	Parelles de partícules
      k = 0
	DO i=1,Natoms-1
	  DO j=i+1,Natoms
		k = k + 1
		pairs(k,1) = i 
		pairs(k,2) = j
	  ENDDO
	ENDDO

!	Nombre de parelles per a cada processador
	pair_per_proc = npairs / numproc
	residu = mod(npairs, numproc)
	
!	pair_limits(1,1) = 1
	DO i=1,numproc
		!pair_limits(i,1) = pair_limits(i-1,2) + 1
	  	!pair_limits(i,2) = pair_limits(i,1) + nint(real(i*npairs)/real(numproc)) - 1
	  if ((i-1) .lt. residu) then
	  	pair_limits(i,1) = (pair_per_proc+1)*(i-1) + 1
	  	pair_limits(i,2) = pair_limits(i,1) + pair_per_proc
	  else
	  	pair_limits(i,1) = (pair_per_proc+1)*(i-1) + 1 + residu
	  	pair_limits(i,2) = pair_limits(i,1) + pair_per_proc - 1	
	  endif	
	  
	  	
	ENDDO   
	
	
	
    END SUBROUTINE
END module
