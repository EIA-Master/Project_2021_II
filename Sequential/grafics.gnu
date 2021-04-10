# Author: Jaume Garcia

set term png
set output "energies.png"
set xlabel "time (dimensionless)"
set ylabel "energy (dimensionless)"
plot "Thermodynamics.dat" i 0 u 1:2 w l t 'Kinetic',"Thermodynamics.dat" i 0 u 1:3 w l t 'Potential',"Thermodynamics.dat" i 0 u 1:4 w l t 'Total'
pause -1

set output "temperature.png"
set xlabel "time (dimensionless)"
set ylabel "temperature (dimensionless)"
plot "Thermodynamics.dat" i 0 u 1:5 w l notitle
pause -1

set output "pressure.png"
set xlabel "time (dimensionless)"
set ylabel "pressure (dimensionless)"
plot "Thermodynamics.dat" i 0 u 1:6 w l notitle
pause -1

set output "radial_distribution.png"
set xlabel "distance (dimensionless)"
set ylabel "g(r) (dimensionless)"
plot "gdr.dat" u 1:2 w l notitle
pause -1
