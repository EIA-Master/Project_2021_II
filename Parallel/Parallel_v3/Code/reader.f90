! Author: Ignasi Puch
module reader 
! Module that contains subroutine to read all the data in data.txt file
contains 

! ------------------------------------------------------------------ !
!                             FILE READER                            !
! ------------------------------------------------------------------ !

subroutine read_files(temp,density,timestep,sigma,epsilon,mass,Natoms,Nsteps,Nradial)
implicit none
! MPI variables: 
integer:: taskid,comm,ierror,reslen,message,partner,request
! Input:
character*8 inputfile 
! Other variables
character:: dummie_variable
! Output
double precision:: temp,density,timestep,sigma,epsilon,mass
integer:: Natoms,Nsteps,Nradial
! ****************************************************************** !
! This subroutine reads all the data in the data.txt file necessary
! in order to perform the simulation with the experimental values
! wanted.
! ****************************************************************** !

inputfile = 'data.txt'

open(unit = 50, file = inputfile)

read(50,*) dummie_variable,temp
read(50,*) dummie_variable,density
read(50,*) dummie_variable,timestep
read(50,*) dummie_variable,sigma
read(50,*) dummie_variable,epsilon
read(50,*) dummie_variable,mass
read(50,*) dummie_variable,Natoms
read(50,*) dummie_variable,Nsteps
read(50,*) dummie_variable,Nradial

close(50)

end subroutine read_files

end module reader

