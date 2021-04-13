# Author: Jaume Garcia
# Helper: Ignasi Puch

# Variables to dimesionalize:
mass = 4.0026
epsilon = 10.9
sigma = 2.64

# ----------------------------------------------------------------------------------- #
#                               ENERGIES                                              #
# ----------------------------------------------------------------------------------- #
set term pngcairo
set output "energies.png"

set title 'Energies'
set xlabel "time (ns)"
set ylabel "energy (kJ/mol)"

set xrange[*:12]
set yrange[*:*]

set key inside right center

plot "Thermodynamics.dat" i 0 u ($1)*sqrt(mass*sigma*sigma/(epsilon*0.8314)):($2)*epsilon*0.008314\
 w l t 'Kinetic',"Thermodynamics.dat" i 0 u ($1)*sqrt(mass*sigma*sigma/(epsilon*0.8314)):($3)*epsilon*0.008314\
 w l t 'Potential',"Thermodynamics.dat" i 0 u ($1)*sqrt(mass*sigma*sigma/(epsilon*0.8314)):($4)*epsilon*0.008314\
 w l t 'Total'

# ----------------------------------------------------------------------------------- #
#                               TEMPERATURE                                           #
# ----------------------------------------------------------------------------------- #
set output "temperature.png"

set title 'Temperature'
set xlabel "time (ns)"
set ylabel "temperature (K)"

set xrange[*:12]
set yrange[*:*]

plot "Thermodynamics.dat" i 0 u ($1)*sqrt(mass*sigma*sigma/(epsilon*0.8314)):($5)*epsilon\
 w l notitle lc rgb 'blue'

# ----------------------------------------------------------------------------------- #
#                                PRESSURE                                             #
# ----------------------------------------------------------------------------------- #
set output "pressure.png"

set title 'Pressure'
set xlabel "time (ns)"
set ylabel "pressure (MPa)"

set xrange[*:12]
set yrange[*:*]

plot "Thermodynamics.dat" i 0 u \
 ($1)*sqrt(mass*sigma*sigma/(epsilon*0.8314)):($6)*epsilon*1.38*10/(sigma*sigma*sigma)\
 w l notitle lc rgb 'red'


# ----------------------------------------------------------------------------------- #
#                            RADIAL DISTRIBUTION                                      #
# ----------------------------------------------------------------------------------- #
set output "radial_distribution.png"

set title 'Radial distribution'
set xlabel "distance (Angstroms)"
set ylabel "g(r) (dimensionless)"

set xrange[*:18]
set yrange[*:*]

plot "gdr.dat" u ($1)*sigma:($2) w l notitle lc rgb 'black'

