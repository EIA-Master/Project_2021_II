module initialize
implicit none
contains


subroutine initial(Npart,ro,temp,posi,vel,L)
! given a number of particles of a system of density ro and temperature temp,
! this subroutine returns the position and velocities of those particles.
implicit none
double precision posi(3,Npart),vel(3,Npart),aver(3),aver2,Ro,L,A,temp,x1,x2,phi,pi
double precision p,v,a1,a2
integer Npart,ii,jj,kk,ll,atom,irank
integer comm,rank,Nproc,Nppp,ierror,icomm

include 'mpif.h'

! Each processor will compute the position and velocities of a certain number of
! particles, Nppp (Number of Particles Per Processor)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierror)

Nppp=Npart/Nproc

L=(Npart/ro)**(1./3.)
pi=4.d0*datan(1.d0)
aver=0.d0
aver2=0.d0

call RANDOM_SEED()
ii=-1
jj=0
posi=0.d0
vel=0.d0
! Each processor wil start from a different particle and will compute the
! position and velocity of Nppp particles
atom=Nppp*rank
do ll=1,Nppp
    atom=atom+1
    ii=ii+1
    if (ii*L/dsqrt(dble(Nppp+1)).gt.L) then
        ii=0
        jj=jj+1
    endif

    posi(1,atom)=L*ii/dsqrt(dble(Nppp+1))-L/2.d0
    posi(2,atom)=L*jj/dsqrt(dble(Nppp+1))-L/2.d0
    posi(3,atom)=L*dble(rank)/dble(Nproc+1)-L/2.d0

    
    ! The velocities are initialized in a MB distribution
    do kk=1,3
        call RANDOM_NUMBER(x1)
        call RANDOM_NUMBER(x2)
        x1=1.d0-x1
        x2=1.d0-x2
        phi=2.d0*pi*x2
        vel(kk,atom)=dsqrt(-2.d0*temp*dlog(x1))*dsin(phi)

        ! We compute the averages so we can rescale the velocities later
        aver(kk)=aver(kk)+vel(kk,atom)
        aver2=aver2+vel(kk,atom)**2.d0
    enddo
    
enddo

! Since normally Nppp*Nproc < Npart, we will have some missing particles. Those
! ones will be computed by the main processos (rank = 0) using the same
! algorithm as before
if (rank.eq.0)then
    atom=Nppp*Nproc
    ii=-1
    jj=0
    do while (atom.lt.Npart)
        atom=atom+1

        ii=ii+1
        if (ii.gt.(dsqrt(dble(Nppp))+1))then
            ii=0
            jj=jj+1
        endif

        posi(1,atom)=L*ii/dsqrt(dble(Nppp+1))-L/2.d0
        posi(2,atom)=L*jj/dsqrt(dble(Nppp+1))-L/2.d0
        posi(3,atom)=L*dble(Nproc)/dble(Nproc+1)-L/2.d0

        
        do kk=1,3
            call RANDOM_NUMBER(x1)
            call RANDOM_NUMBER(x2)               
            x1=1.d0-x1
            x2=1.d0-x2
            phi=2.d0*pi*x2

            vel(kK,atom)=dsqrt(-2.d0*temp*dlog(x1))*dsin(phi)
            aver(kk)=aver(kk)+vel(kk,atom)
            aver2=aver2+vel(kk,atom)**2.d0
        enddo                                                 
        
    enddo
    
endif

! We wait until all positions and velocities are computed 
call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! We add all the positions and velocities, so all the processors have all the
! information
do ii=1,Npart
    do kk=1,3
        call MPI_ALLREDUCE(posi(kk,ii),p,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        call MPI_ALLREDUCE(vel(kk,ii),v,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        
        posi(kk,ii)=p
        vel(kk,ii)=v

    enddo
enddo

! We also add the averages together, so we can re-scale the velocities
do kk=1,3
    call MPI_ALLREDUCE(aver(kk),a1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    aver(kk)=a1
enddo
    
call MPI_ALLREDUCE(aver2,a2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
aver2=a2

! Re-scaling
do ii=1,Npart
    do kk=1,3
        vel(kk,ii)=vel(kk,ii)-aver(kk)/Npart
    enddo
enddo

vel=vel*dsqrt(3.d0*Npart*temp/aver2)


call MPI_BARRIER(MPI_COMM_WORLD,ierror)

return
end subroutine initial



end module initialize
