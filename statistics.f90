module statistics
implicit none
contains

subroutine statistics(Natoms,posi,vel,Kin)
implicit none

double precision posi(3,Natoms),vel(3,Natoms),Kin
integer Natoms,II,JJ,KK


Kin=0.d0
do II=1,Natoms
    do JJ=1,3
        Kin=Kin+Vel(JJ,II)**2.d0
    enddo
enddo
Kin=Kin*0.5d0



return
end