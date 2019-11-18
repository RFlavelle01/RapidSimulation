# Set up the compiler
FC=gfortran

# Set the compiler flags
#FCFLAGS= -O3 -Wall -fcheck=all -Wno-unused-dummy-argument
FCFLAGS= -O3 -Wall -Wno-unused-dummy-argument

# Set up any libraries needed
LFLAGS=

# Set up the hooks
all: mesher execute

mesher: precision.o mesh_generator.o
	$(FC) $(FCFLAGS) precision.o mesh_generator.o -o $@ $(LFLAGS)
	mv $@ ../bin/.

execute: precision.o mesh.o equations.o flow.o work.o io.o main.o
	$(FC) $(FCFLAGS) precision.o mesh.o equations.o flow.o work.o io.o main.o -o $@ $(LFLAGS)
	mv $@ ../bin/.

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod rm *~ 

