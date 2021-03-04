# Start of the makefile
# Defining variables

objects = forces.o main.o Initialize.o Integration.o boundary.o statistics.o radial_distribution.o parameters.o
f90comp = gfortran

OPT= -O3
# Makefile
programa.exe: $(objects)
	$(f90comp) -o programa.exe $(OPT) $(objects)

boundary.o: boundary.f90
	$(f90comp) -c $(OPT) boundary.f90
parameters.o: parameters.f90
	$(f90comp) -c $(OPT) parameters.f90
Initialize.o: Initialize.f90
	$(f90comp) -c $(OPT) Initialize.f90
forces.o: boundary.o Initialize.o forces.f90
	$(f90comp) -c $(OPT) forces.f90	
statistics.o: statistics.f90
	$(f90comp) -c $(OPT) statistics.f90	
Integration.o: boundary.o forces.o statistics.o radial_distribution.o Integration.f90
	$(f90comp) -c $(OPT) Integration.f90
radial_distribution.o: boundary.o radial_distribution.f90
	$(f90comp) -c $(OPT) radial_distribution.f90	
main.o: Initialize.o Integration.o parameters.o main.f90
	$(f90comp) -c $(OPT) main.f90
  
# Cleaning everything
clean:
	rm -f programa.exe
	rm -f $(objects)
# End of the makefile
