#variables
FC=mpifort
CFLAGS=-c -g -Og -Wall

## Turn on for debugging
DFLAGS=-fbacktrace -fmax-errors=5

#linking
a.out: mobcal_jay.o
	$(FC)  commons.o xtb_mctc_accuracy.o MobCal_MPI.o

#compiling
mobcal_jay.o: MobCal_MPI.F90 #commons.f90 xtb_mctc_accuracy.f90
	$(FC) $(CFLAGS) $(DFLAGS) commons.f90 xtb_mctc_accuracy.f90 MobCal_MPI.F90

#cleanup
#clean:
#	rm *.o a.out

#run
run:
	make
	./a.out
