module statistics
! Module with functions and subroutines to calculate averages.
implicit none
contains

! ------------------------------------------------------------------ !
!                           KINETIC ENERGY                           !
! ------------------------------------------------------------------ !
subroutine kinetic(Natoms,vel,kin)
implicit none
! Input
integer Natoms
double precision vel(3,Natoms)
! Output
double precision kin
! Other variables
integer II,JJ,KK
! ------------------------------------------------------------------ !

kin = 0.d0
do II = 1,Natoms
    do JJ = 1,3

        kin = kin + vel(JJ,II)**2.d0

    enddo
enddo
kin = kin*0.5d0


return
end subroutine kinetic

! ------------------------------------------------------------------ !
!                       INSTANT TEMPERATURE                          !
! ------------------------------------------------------------------ !

<<<<<<< HEAD
double precision function insttemp(Natoms,kine)
=======
subroutine insttemp(Natoms,kine,temp)
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c
implicit none
! Input
integer Natoms
double precision kine
<<<<<<< HEAD
=======
! Output
double precision temp
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c
! Other variables
integer Nf
! ------------------------------------------------------------------ !

Nf = 3*Natoms - 3
temp = 2.d0*kine/dble(Nf)

<<<<<<< HEAD
end function insttemp
=======
end subroutine insttemp
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c

! ------------------------------------------------------------------ !
!                           TOTAL ENERGY                             !
! ------------------------------------------------------------------ !

double precision function totalenergy(potene,kinene)
implicit none
! Input
double precision potene,kinene
! ------------------------------------------------------------------ !

totalenergy = potene + kinene

end function totalenergy


! ------------------------------------------------------------------ !
!                              PRESSURE                              !
! ------------------------------------------------------------------ !

subroutine pressure(Natoms,L,rho,posis,force,temp,pres)
implicit none
! Input
integer Natoms
double precision L,rho,posis(3,Natoms),force(3,Natoms),temp
! Output
double precision pres
! Other variables
integer i,j
! ------------------------------------------------------------------ !

pres = 0.d0

do i = 1, Natoms
    do j = 1, 3

        pres = pres + posis(j,i)*force(j,i)

    enddo
enddo

pres = rho*temp + pres/(3.d0*(L**3)*Natoms)

end subroutine pressure

end module statistics
