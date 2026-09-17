TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS =

MAINP = stokes stokes_3D

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

stokes: $(OBJS) stokes.o
	$(FC) $(FFLAGS) -o $@ stokes.o $(OBJS) $(LIBS)

stokes_3D: $(OBJS) stokes_3D.o
	$(FC) $(FFLAGS) -o $@ stokes_3D.o $(OBJS) $(LIBS)

clean:
	$(RM)  stokes stokes_3D *.o *.mod *.bin tags

outclean:
	$(RM) $(OUTCLEAN)

allclean: clean outclean
