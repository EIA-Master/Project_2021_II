# Start of the makefile
# Defining variables
objects = forces.o main.o Initialize.o Integration.o boundary.o statistics.o
f90comp = gfortran
OPT= -O3
# Makefile
programa.x: $(objects)
	$(f90comp) -o programa.x $(OPT) $(objects)
#boundary.mod: boundary.f90
#	$(f90comp) -c $(OPT) boundary.f90
boundary.o: boundary.f90
	$(f90comp) -c $(OPT) boundary.f90
#forces.mod: forces.f90
#	$(f90comp) -c $(OPT) forces.f90
#Initialize.mod: Initialize.f90
#	$(f90comp) -c $(OPT) Initialize.f90
Initialize.o: Initialize.f90
	$(f90comp) -c $(OPT) Initialize.f90
forces.o: boundary.o Initialize.o forces.f90
	$(f90comp) -c $(OPT) forces.f90	
#statistics.mod: statistics.f90
#	$(f90comp) -c $(OPT) statistics.f90
statistics.o: statistics.f90
	$(f90comp) -c $(OPT) statistics.f90	
#Integration.mod: Integration.f90
#	$(f90comp) -c $(OPT) Integration.f90
Integration.o: boundary.o forces.o statistics.o Integration.f90
	$(f90comp) -c $(OPT) Integration.f90	
main.o: Initialize.o Integration.o main.f90
	$(f90comp) -c $(OPT) main.f90

.PHONY:plots
plots: gnuplot grafics.gnu
# Cleaning everything
clean:
	rm -f programa.x
	rm -f $(objects)
# End of the makefile
