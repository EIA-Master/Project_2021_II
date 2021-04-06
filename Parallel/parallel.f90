module parallel
    IMPLICIT NONE

    contains
    SUBROUTINE assignacio(Natoms,numproc)
        implicit none
        integer :: i,particle_per_proc,residu
        integer :: Natoms, numproc
        integer :: index_particles(:,:)
        integer :: allgather_pointer(:)
        integer :: num_send(:)
        
        particle_per_proc = int(real(Natoms)/real(numproc))
        residu = mod(Natoms,numproc)
        
        allocate(index_particles(numproc,2),allgather_pointer(numproc),num_send(numproc))

!	Nombre de partícules que li corresponen a cada processador
            DO i=1,numproc
              IF((i-1)<residu)then
                index_particles(i,1)=(particle_per_proc+1)*(i-1)+1
                index_particles(i,2)=index_particles(i,1)+particle_per_proc
              else
              index_particles(i,1)=(i-1)*particle_per_proc+res+1
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
    END SUBROUTINE
END module
