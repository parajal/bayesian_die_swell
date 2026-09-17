TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS =

MAINP = stokes11par

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

stokes11par: $(OBJS) stokes11par.o
	$(FC) $(FFLAGS) -o $@ stokes11par.o $(OBJS) $(LIBS)

clean:
	$(RM)  stokes11par *.o *.mod *.bin tags

outclean:
	$(RM) $(OUTCLEAN)

allclean: clean outclean
