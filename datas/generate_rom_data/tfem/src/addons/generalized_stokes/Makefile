TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = generalized_newtonian_elements.o generalized_newtonian_set_globals.o generalized_stokes_elements.o generalized_stokes_globals.o generalized_stokes_set_globals.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

generalized_newtonian_elements.o: generalized_stokes_globals.o

generalized_newtonian_set_globals.o: generalized_newtonian_elements.o

generalized_stokes_elements.o: generalized_newtonian_elements.o generalized_stokes_globals.o

generalized_stokes_set_globals.o: generalized_stokes_elements.o

clean:
	$(RM) *.o tags
