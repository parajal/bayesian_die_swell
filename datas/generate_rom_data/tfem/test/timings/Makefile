TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS =

MAINP = stokes_ma41 stokes_ma57

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

stokes_ma41: $(OBJS) stokes_ma41.o
	$(FC) $(FFLAGS) -o $@ stokes_ma41.o $(OBJS) $(LIBS)

stokes_ma57: $(OBJS) stokes_ma57.o
	$(FC) $(FFLAGS) -o $@ stokes_ma57.o $(OBJS) $(LIBS)

clean:
	$(RM)  stokes_ma41 stokes_ma57 *.o *.mod *.bin tags
