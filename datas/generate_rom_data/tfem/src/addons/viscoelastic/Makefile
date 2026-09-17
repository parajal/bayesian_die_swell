TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = viscoelastic_elements.o viscoelastic_elements_DG_2D.o viscoelastic_elements_SUPG.o viscoelastic_elements_embedded.o viscoelastic_elements_gammap.o viscoelastic_elements_generic.o viscoelastic_elements_implicit.o viscoelastic_elements_integrate.o viscoelastic_globals.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

viscoelastic_elements.o: viscoelastic_elements_DG_2D.o viscoelastic_elements_SUPG.o viscoelastic_elements_embedded.o viscoelastic_elements_gammap.o viscoelastic_elements_implicit.o viscoelastic_elements_integrate.o

viscoelastic_elements_DG_2D.o: viscoelastic_elements_generic.o viscoelastic_globals.o

viscoelastic_elements_SUPG.o: viscoelastic_elements_generic.o viscoelastic_globals.o

viscoelastic_elements_embedded.o: viscoelastic_elements_generic.o viscoelastic_globals.o

viscoelastic_elements_gammap.o: viscoelastic_elements_generic.o viscoelastic_globals.o

viscoelastic_elements_generic.o: viscoelastic_globals.o

viscoelastic_elements_implicit.o: viscoelastic_elements_SUPG.o viscoelastic_elements_generic.o viscoelastic_globals.o

viscoelastic_elements_integrate.o: viscoelastic_elements_generic.o viscoelastic_globals.o

clean:
	$(RM) *.o tags
