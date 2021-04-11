# Author: Jaume Garcia

! ----------------------------------------------------------------------------------- !
!                               ENERGIES                                              !
! ----------------------------------------------------------------------------------- !
set term eps
set output "energies.eps"
set title 'Energies' temporal evolution'
set xlabel "time (dimensionless)"
set ylabel "energy (dimensionless)"
set key outside

plot "Thermodynamics.dat" i 0 u 1:2 w l t 'Kinetic',"Thermodynamics.dat"\
 i 0 u 1:3 w l t 'Potential',"Thermodynamics.dat" i 0 u 1:4 w l t 'Total'

! ----------------------------------------------------------------------------------- !
!                               TEMPERATURE                                           !
! ----------------------------------------------------------------------------------- !
set output "temperature.eps"
set title 'Temperature's temporal evolution'
set xlabel "time (dimensionless)"
set ylabel "temperature (dimensionless)"
set key outside

plot "Thermodynamics.dat" i 0 u 1:5 w l notitle

! ----------------------------------------------------------------------------------- !
!                                PRESSURE                                             !
! ----------------------------------------------------------------------------------- !
set output "pressure.eps"
set title 'Pressure's temporal evolution'
set xlabel "time (dimensionless)"
set ylabel "pressure (dimensionless)"
set key outside

plot "Thermodynamics.dat" i 0 u 1:6 w l notitle

! ----------------------------------------------------------------------------------- !
!                            RADIAL DISTRIBUTION                                      !
! ----------------------------------------------------------------------------------- !
set output "radial_distribution.eps"
set title 'Radial distribution'
set xlabel "distance (dimensionless)"
set ylabel "g(r) (dimensionless)"
set key outside

plot "gdr.dat" u 1:2 w l notitle

