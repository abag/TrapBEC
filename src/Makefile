F90=openmpif90
FFLAGS=-O3   # SP/ 3rd level of optimisation
   #
   #  -pg   =  profiler (gprof)
   #  -g    = debugger  (gdb)
   #

## Generic rules:
.SUFFIXES:          # get rid of that annoying Modula rule
.SUFFIXES: .f .f90 .c .o .obj .x .exe .inc .h  

#
# Rules (how to make a Blergh from a Blah
#
.f90.o:
	    $(F90) $(FFLAGS) -o $*.o -c $*.f90

#
# Special targets
#
default: run.x
#
# Targets
#
run.x: run.f90 cdata.o mpi_var.o statistics.o initial.o ghost.o deriv.o finite_diff.o output.o diagnostic.o remesh.o Makefile
	$(F90) $(FFLAGS) -o ./run.x run.f90 statistics.o cdata.o mpi_var.o initial.o ghost.o deriv.o finite_diff.o output.o diagnostic.o remesh.o
#
#
clean: 
	rm -f *.o *.x *.exe *.obj *.mod
cleann: clean
	rm -f *~
pristine: cleann
	rm -f *.dat *.log *.png
#
# Dependencies
#
cdata.o: cdata.f90 cparam.local
# 
mpi_var.o: mpi_var.f90 cdata.o 
#
ghost.o: ghost.f90 cdata.o mpi_var.o 
#
output.o: output.f90 cdata.o
#
deriv.o: deriv.f90 cdata.o
#
statistics.o: statistics.f90 cdata.o
#
finite_diff.o: finite_diff.f90 mpi_var.o deriv.o cdata.o
#
initial.o: initial.f90 statistics.o cdata.o
#
diagnostic.o: diagnostic.f90 cdata.o deriv.o
#
remesh.o: remesh.f90 cdata.o mpi_var.o 
