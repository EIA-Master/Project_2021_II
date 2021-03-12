# Start of the makefile
# Defining variables

<<<<<<< HEAD
objects = forces.o main.o Initialize.o Integration.o boundary.o statistics.o radial_distribution.o parameters.o
=======
objects = forces.o main.o Initialize.o Integration.o boundary.o statistics.o parameters.o
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c
f90comp = gfortran

OPT= -O3
# Makefile
programa.exe: $(objects)
	$(f90comp) -o programa.exe $(OPT) $(objects)
<<<<<<< HEAD
=======

>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c
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
<<<<<<< HEAD
Integration.o: boundary.o forces.o statistics.o radial_distribution.o Integration.f90
	$(f90comp) -c $(OPT) Integration.f90
radial_distribution.o: boundary.o radial_distribution.f90
	$(f90comp) -c $(OPT) radial_distribution.f90	
main.o: Initialize.o Integration.o parameters.o main.f90
	$(f90comp) -c $(OPT) main.f90
  
## clean: Remove *.o , *.mod and *.exe files
.PHONY: clean
clean:
	rm -f programa.exe
	rm -f $(objects)
	rm -f *.mod

## plots: Generate gnuplots
.PHONY: plots
plots:
	gnuplot grafics.gnu
	
#Help
.PHONY: help
help:
	@sed -n 's/^##//p' Makefile
=======
Integration.o: boundary.o forces.o statistics.o Integration.f90
	$(f90comp) -c $(OPT) Integration.f90	
main.o: Initialize.o Integration.o parameters.o main.f90
	$(f90comp) -c $(OPT) main.f90
  
# Cleaning everything
clean:
	rm -f programa.exe
	rm -f $(objects)
>>>>>>> f7f33d3c2708777b005facdce400c724a0909c0c
# End of the makefile
