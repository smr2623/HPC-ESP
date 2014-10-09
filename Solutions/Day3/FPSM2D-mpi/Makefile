
FC = mpif90
#FC = gfortran

#SunFire libraries:
FFTWlibdir = -L/opt/lib        
FFTWlib = -lfftw3f

FFTWincdir= -I/opt/fftw-3.2.2/include

#SunFire options:
FOPTS = -fastsse -O3 -Mipa=fast -mp -mcmodel=medium
FOPTS = -O0 -g -C -Ktrap=divz -mp 
FOPTS = -g –Ktrap=fp -C

#Fermi gnu options
#FOPTS = -fopenmp -O3

#######################################

MYHOME= /home/pklin
FEXEC = $(MYHOME)/bin/fpsm2D.mpi
F = FPSM2D_

SRCDIR= $(MYHOME)/src/Fury/f90
SF13=$(SRCDIR)/$(F13)

.SUFFIXES: .f90 .o

#OBJ=$(SRC:.f90=.o)
OBJ = FPSM2D_MAIN.o $(F)Pinput.o $(F)Simulation.o $(F)Wave.o \
$(F)Output.o $(F)Operator.o $(F)Deriv.o DFourier.o \
$(F)Geninput.o $(F)Ddata.o $(F)Fileunits.o $(F)Parallelism.o
 
all : $(FEXEC)

$(FEXEC): $(OBJ)
	$(FC) $(FOPTS) -o $@ $(OBJ) $(FFTWlibdir) $(FFTWlib)
	
FPSM2D_MAIN.o: FPSM2D_MAIN.f90 $(F)Simulation.o
	$(FC) $(FOPTS) -c FPSM2D_MAIN.f90
$(F)Simulation.o: $(F)Simulation.f90 $(F)Pinput.o $(F)Wave.o 
	$(FC) $(FOPTS) -c $(F)Simulation.f90
$(F)Wave.o: $(F)Wave.f90 $(F)Pinput.o $(F)Output.o $(F)Operator.o
	$(FC) $(FOPTS) -c $(F)Wave.f90
$(F)Output.o: $(F)Output.f90 $(F)Pinput.o $(F)Fileunits.o $(F)Parallelism.o
	$(FC) $(FOPTS) -c $(F)Output.f90
$(F)Operator.o: $(F)Operator.f90 $(F)Pinput.o $(F)Deriv.o
	$(FC) $(FOPTS) -c $(F)Operator.f90
$(F)Deriv.o: $(F)Deriv.f90 DFourier.o $(F)Pinput.o $(F)Parallelism.o
	$(FC) $(FOPTS) -c $(F)Deriv.f90 
$(F)Pinput.o: $(F)Pinput.f90 $(F)Geninput.o $(F)Ddata.o $(F)Parallelism.o
	$(FC) $(FOPTS) -c $(F)Pinput.f90
$(F)Parallelism.o: $(F)Parallelism.f90
	$(FC) $(FOPTS) -c $(F)Parallelism.f90
$(F)Geninput.o: $(F)Geninput.f90 $(F)Fileunits.o
	$(FC) $(FOPTS) -c $(F)Geninput.f90
$(F)Ddata.o: $(F)Ddata.f90 $(F)Fileunits.o
	$(FC) $(FOPTS) -c $(F)Ddata.f90 
$(F)Fileunits.o: $(F)Fileunits.f90
	$(FC) $(FOPTS) -c $(F)Fileunits.f90
DFourier.o: DFourier.f90
	$(FC) $(FOPTS) -c DFourier.f90 $(FFTWincdir)

clean : 
	rm -rf *.o
